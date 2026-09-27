#!/usr/bin/env python3
"""Generate random nonoverlapping oblate spheroids in a fully periodic box.

Python standard library only. Positions and orientations are proposed randomly.
Periodic cell lists limit collision queries to nearby particles. Successful
insertion never moves existing particles; consecutive failures trigger bounded
geometric Monte Carlo moves of the blockers and their nearest neighbours only.
These are overlap-rejecting coordinate proposals, with no physical time,
dynamics or aging. Optional final MC sweeps retain their existing meaning.
The CSV uses OpenLB's body-to-world Rz(angle_z) Ry(angle_y) Rx(angle_x), degrees.
Contact/clearance uses Perram--Wertheim plus a conservative separating plane:
J. Comput. Phys. 58, 409--416 (1985), doi:10.1016/0021-9991(85)90171-8.
"""

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import random
import sys
import time


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def unit_quaternion(rng):
    """Uniform SO(3), stored as (w,x,y,z); not uniform Euler angles."""
    u, v, w = rng.random(), rng.random(), rng.random()
    return (math.sqrt(u) * math.cos(2 * math.pi * w),
            math.sqrt(1 - u) * math.sin(2 * math.pi * v),
            math.sqrt(1 - u) * math.cos(2 * math.pi * v),
            math.sqrt(u) * math.sin(2 * math.pi * w))


def quat_product(p, q):
    w, x, y, z = p
    a, b, c, d = q
    return (w*a-x*b-y*c-z*d, w*b+x*a+y*d-z*c,
            w*c-x*d+y*a+z*b, w*d+x*c-y*b+z*a)


def quat_matrix(q):
    w, x, y, z = q
    return ((1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)),
            (2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)),
            (2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)))


def euler_degrees(q):
    r = quat_matrix(q)
    # Random orientations almost surely avoid gimbal lock; handle it anyway.
    ay = math.asin(max(-1.0, min(1.0, -r[2][0])))
    if abs(math.cos(ay)) > 1.e-12:
        ax, az = math.atan2(r[2][1], r[2][2]), math.atan2(r[1][0], r[0][0])
    else:
        ax, az = math.atan2(-r[1][2], r[1][1]), 0.0
    return tuple(math.degrees(t) for t in (ax, ay, az))


def particle(position, quaternion):
    qnorm = math.sqrt(dot(quaternion, quaternion))
    q = tuple(t / qnorm for t in quaternion)
    matrix = quat_matrix(q)
    return (tuple(position), q, tuple(row[2] for row in matrix))


def support(normal, direction, a, c):
    """Support function, valid for non-unit direction too."""
    return math.sqrt(max(0.0, a*a*dot(direction, direction)
                         + (c*c-a*a)*dot(normal, direction)**2))


def matrix_coefficients(n, a, c):
    b = c*c - a*a
    x, y, z = n
    return (a*a+b*x*x, a*a+b*y*y, a*a+b*z*z, b*x*y, b*x*z, b*y*z)


def solve_symmetric(m, r):
    xx, yy, zz, xy, xz, yz = m
    c00, c11, c22 = yy*zz-yz*yz, xx*zz-xz*xz, xx*yy-xy*xy
    c01, c02, c12 = xz*yz-xy*zz, xy*yz-xz*yy, xy*xz-xx*yz
    det = xx*c00 + xy*c01 + xz*c02
    x, y, z = r
    return ((c00*x+c01*y+c02*z)/det,
            (c01*x+c11*y+c12*z)/det,
            (c02*x+c12*y+c22*z)/det)


def pw_contact(r, n1, n2, a, c, iterations=22):
    """Return max_lambda F and its separating direction.

    A_i=R_i diag(a^2,a^2,c^2) R_i^T,
    F=lambda(1-lambda) r^T[(1-lambda)A_1+lambda A_2]^{-1}r.
    Fmax < 1 means overlap, =1 tangency, >1 separated.  Golden-section
    maximization is conservative for rejection: an underestimated F rejects a
    valid placement rather than accepting an overlapping one.
    """
    a1, a2 = matrix_coefficients(n1, a, c), matrix_coefficients(n2, a, c)

    def evaluate(t):
        m = tuple((1-t)*v + t*w for v, w in zip(a1, a2))
        direction = solve_symmetric(m, r)
        return t*(1-t)*dot(r, direction), direction

    lo, hi = 0.0, 1.0
    ratio = (math.sqrt(5.0)-1.0)/2.0
    left, right = hi-ratio*(hi-lo), lo+ratio*(hi-lo)
    fl, fr = evaluate(left), evaluate(right)
    for _ in range(iterations):
        if fl[0] < fr[0]:
            lo, left, fl = left, right, fr
            right = lo+ratio*(hi-lo)
            fr = evaluate(right)
        else:
            hi, right, fr = right, left, fl
            left = hi-ratio*(hi-lo)
            fl = evaluate(left)
    return fl if fl[0] >= fr[0] else fr


def separated(r, n1, n2, a, c, gap):
    distance2 = dot(r, r)
    if distance2 >= (2*a+gap)**2:
        return True
    if distance2 < (2*c+gap)**2:
        return False
    distance = math.sqrt(distance2)
    # Any separating plane with the requested margin proves clearance.
    if distance2 - support(n1, r, a, c) - support(n2, r, a, c) >= gap*distance:
        return True
    # Thin-axis planes often separate nearby platelets without a PW solve.
    # A successful plane test is a geometric proof, not a distance surrogate.
    for axis in (n1, n2):
        if abs(dot(r, axis))-support(n1, axis, a, c)-support(n2, axis, a, c) >= gap:
            return True
    fmax, direction = pw_contact(r, n1, n2, a, c)
    if fmax <= 1.0 + 1.e-12:
        return False
    norm = math.sqrt(dot(direction, direction))
    margin = (dot(r, direction)-support(n1, direction, a, c)
              -support(n2, direction, a, c))/norm
    return margin >= gap - 1.e-12*a


def image_displacements(p, q, length, cutoff):
    """Relevant periodic images in all three axes, at initial shear phase zero."""
    dr = tuple(q[i]-p[i] for i in range(3))
    if length > 2*cutoff:
        # At most one image per pair can fall within the bounding-sphere range.
        dr = tuple(d-length*math.floor(d/length+0.5) for d in dr)
        if all(abs(d) <= cutoff for d in dr):
            yield dr
        return
    ranges = [range(math.ceil((-cutoff-d)/length),
                    math.floor((cutoff-d)/length)+1) for d in dr]
    for image in itertools.product(*ranges):
        yield tuple(dr[i]+image[i]*length for i in range(3))


class PeriodicCellList:
    """Mutable centre index; each cell is at least one collision cutoff wide.

    All 27 neighbouring cells are sufficient, including wrapped boundary
    cells. Deduplication matters when there are only one or two cells per axis.
    The index is only a broad phase: the existing spheroid predicate and every
    relevant periodic image still decide overlap/clearance.
    """

    def __init__(self, length, cutoff, bodies=()):
        if not (math.isfinite(length) and math.isfinite(cutoff)
                and length > 0 and cutoff > 0):
            raise ValueError('Cell-list length and cutoff must be positive and finite')
        self.length = length
        self.count = max(1, int(length/cutoff))
        self.width = length/self.count
        self.bins = {}
        self.locations = {}
        self.neighbour_keys = {}
        for index, body in enumerate(bodies):
            self.add(index, body[0])

    def key(self, position):
        return tuple(min(self.count-1, int((x % self.length)/self.width))
                     for x in position)

    def add(self, index, position):
        if index in self.locations:
            raise ValueError('Particle already in cell list: %d' % index)
        key = self.key(position)
        self.locations[index] = key
        self.bins.setdefault(key, set()).add(index)

    def remove(self, index):
        key = self.locations.pop(index)
        self.bins[key].remove(index)
        if not self.bins[key]:
            del self.bins[key]

    def move(self, index, position):
        if self.key(position) != self.locations[index]:
            self.remove(index)
            self.add(index, position)

    def neighbors(self, position):
        key = self.key(position)
        if key not in self.neighbour_keys:
            self.neighbour_keys[key] = sorted(set(
                tuple((key[d]+offset[d]) % self.count for d in range(3))
                for offset in itertools.product((-1, 0, 1), repeat=3)))
        # Stable ordering keeps seeded runs independent of set iteration order.
        return sorted(index for cell in self.neighbour_keys[key]
                      for index in self.bins.get(cell, ()))


def collision_blockers(candidate, bodies, length, a, c, gap, skip=-1,
                       spatial=None, limit=None):
    """Return blocking IDs, optionally stopping once a shortlist cannot improve."""
    position, _, normal = candidate
    indices = spatial.neighbors(position) if spatial is not None else range(len(bodies))
    blockers = []
    for j in indices:
        if j == skip:
            continue
        other = bodies[j]
        for dr in image_displacements(position, other[0], length, 2*a+gap):
            if not separated(dr, normal, other[2], a, c, gap):
                blockers.append(j)
                if limit is not None and len(blockers) >= limit:
                    return blockers
                break
    return blockers


def fits(candidate, bodies, length, a, c, gap, skip=-1, spatial=None):
    return not collision_blockers(candidate, bodies, length, a, c, gap,
                                  skip, spatial, limit=1)


def random_candidate(rng, length, a, c, gap):
    return particle(tuple(rng.uniform(0, length) for _ in range(3)),
                    unit_quaternion(rng))


def periodic_distance2(p, q, length):
    return sum((d-length*math.floor(d/length+0.5))**2
               for d in (q[i]-p[i] for i in range(3)))


def local_relax(rng, bodies, spatial, trigger, blockers, length, a, c, gap,
                max_particles, sweeps):
    """Bounded local moves of blockers and their nearest periodic neighbours.

    Every accepted pose is valid against ALL nearby particles, including fixed
    particles outside the movable set. Accepted moves persist even if the
    trigger still cannot be inserted; there is no deletion or overlap repair.
    """
    if len(blockers) > max_particles:
        return 0, 0
    selected = set(blockers)
    nearest = sorted(range(len(bodies)),
                     key=lambda i: (periodic_distance2(trigger[0], bodies[i][0], length), i))
    for index in nearest:
        if len(selected) >= min(max_particles, len(bodies)):
            break
        selected.add(index)
    accepted = attempted = 0
    for _ in range(sweeps):
        acc, att = randomize(rng, bodies, length, a, c, gap, 1,
                             spatial=spatial, indices=sorted(selected), local_only=True)
        accepted += acc
        attempted += att
        if fits(trigger, bodies, length, a, c, gap, spatial=spatial):
            break
    return accepted, attempted


def random_insertion(rng, count, length, a, c, gap, attempts, batches,
                     stagnation_attempts=128, local_particles=16, local_sweeps=8,
                     stats=None):
    """Continuous random insertion, with local moves only after stagnation."""
    bodies = []
    spatial = PeriodicCellList(length, 2*a+gap)
    accepted = attempted = proposals = 0
    if stats is None:
        stats = {}
    stats.update(local_relaxations=0, local_mc_accepted=0, local_mc_attempted=0,
                 max_local_particles=0)
    growth_steps = max(0, (local_particles-1).bit_length()-2)
    start = last_report = time.monotonic()
    for index in range(count):
        inserted = False
        consecutive = rescues = 0
        best = None
        for trial in range(attempts*batches):
            proposals += 1
            candidate = random_candidate(rng, length, a, c, gap)
            # Only a strictly better shortlist can replace the stored trigger.
            # Truncated blocker lists are never used to choose a movable set.
            limit = len(best[1]) if best is not None else local_particles+1
            blocked = collision_blockers(candidate, bodies, length, a, c, gap,
                                         spatial=spatial, limit=limit)
            if not blocked:
                bodies.append(candidate)
                spatial.add(index, candidate[0])
                inserted = True
                break
            consecutive += 1
            if len(blocked) < limit:
                best = (candidate, blocked)
            if consecutive >= stagnation_attempts:
                if best is not None:
                    trigger, blockers = best
                    size = min(local_particles,
                               max(len(blockers), 4*(2**min(rescues, growth_steps))))
                    acc, att = local_relax(rng, bodies, spatial, trigger, blockers,
                                         length, a, c, gap, size, local_sweeps)
                    accepted += acc
                    attempted += att
                    stats['local_relaxations'] += 1
                    stats['max_local_particles'] = max(stats['max_local_particles'],
                                                       min(size, len(bodies)))
                    rescues += 1
                    if fits(trigger, bodies, length, a, c, gap, spatial=spatial):
                        bodies.append(trigger)
                        spatial.add(index, trigger[0])
                        inserted = True
                        break
                best = None
                consecutive = 0
            now = time.monotonic()
            if now-last_report >= 10:
                print('Random placement: %d/%d particles; candidate %d/%d; '
                      'local relaxations %d; elapsed %.1f s'
                      % (index, count, trial+1, attempts*batches,
                         stats['local_relaxations'], now-start), file=sys.stderr, flush=True)
                last_report = now
        if not inserted:
            raise RuntimeError('Random preparation placed %d/%d particles after %d '
                               'candidates and %d local relaxations for this particle; '
                               'no overlapping or underfilled state was written. '
                               'Increase particles.insertion_batches_per_particle or '
                               'particles.insertion_local_sweeps explicitly.'
                               % (index, count, attempts*batches, rescues))
        if (index + 1) % 12 == 0 or index + 1 >= count-12:
            print('Random placement: %d/%d particles; candidates %d; '
                  'local relaxations %d; elapsed %.1f s'
                  % (index+1, count, proposals, stats['local_relaxations'],
                     time.monotonic()-start),
                  file=sys.stderr, flush=True)
            last_report = time.monotonic()
    stats.update(local_mc_accepted=accepted, local_mc_attempted=attempted)
    return bodies, accepted, attempted, proposals


def randomize(rng, bodies, length, a, c, gap, sweeps, spatial=None,
              indices=None, local_only=False, progress=False):
    if spatial is None:
        spatial = PeriodicCellList(length, 2*a+gap, bodies)
    accepted, attempted = 0, 0
    start = last_report = time.monotonic()
    for sweep in range(sweeps):
        order = list(range(len(bodies))) if indices is None else list(indices)
        rng.shuffle(order)
        for index in order:
            attempted += 1
            old = bodies[index]
            # Occasional whole-box relocations use independent uniform SO(3).
            if not local_only and rng.random() < .1:
                candidate = random_candidate(rng, length, a, c, gap)
            else:
                displacement = .12*a
                position = tuple(old[0][i]+rng.uniform(-displacement, displacement)
                                 for i in range(3))
                position = tuple(x % length for x in position)
                axis = tuple(rng.gauss(0, 1) for _ in range(3))
                norm = math.sqrt(dot(axis, axis))
                angle = rng.uniform(-.30, .30)
                dq = (math.cos(angle/2),) + tuple(math.sin(angle/2)*x/norm for x in axis)
                candidate = particle(position, quat_product(dq, old[1]))
            if fits(candidate, bodies, length, a, c, gap, index, spatial):
                bodies[index] = candidate
                spatial.move(index, candidate[0])
                accepted += 1
        now = time.monotonic()
        if progress and (now-last_report >= 10 or sweep+1 == sweeps):
            print('Final geometric MC: %d/%d sweeps; accepted %d/%d; elapsed %.1f s'
                  % (sweep+1, sweeps, accepted, attempted, now-start),
                  file=sys.stderr, flush=True)
            last_report = now
    return accepted, attempted


def validate(bodies, length, a, c, gap):
    for i, body in enumerate(bodies):
        if not fits(body, bodies, length, a, c, gap, i):
            raise RuntimeError("Final overlap/clearance validation failed at particle %d" % i)
    # This implementation intentionally excludes particle self-image contact.
    # The validated restriction below guarantees all such images are separated.
    if length <= 2*a+gap:
        raise ValueError("Box must exceed particle diameter plus clearance.")


def generate(config, output_path):
    geometry, values = config["geometry"], config["particles"]
    physical_length = float(geometry["box_length_m"])
    physical_a, physical_c = float(values["diameter_m"])/2, float(values["thickness_m"])/2
    count, seed = int(values["count"]), int(values.get("seed", 1729))
    physical_gap = float(values.get("minimum_gap_m", 0.0))
    if not (all(math.isfinite(v) for v in (physical_length, physical_a, physical_c, physical_gap))
            and physical_length > 0 and physical_a >= physical_c > 0 and physical_gap >= 0):
        raise ValueError("Require positive box/oblate semiaxes and nonnegative clearance.")
    if count < 0 or count != values["count"]:
        raise ValueError("particles.count must be a nonnegative integer.")
    if physical_length <= 2*(physical_a+physical_gap):
        raise ValueError("Box too small for particle size and periodic clearance.")
    # Normalize by major semiaxis to avoid micron-scale determinant roundoff.
    a, c, length, gap = 1.0, physical_c/physical_a, physical_length/physical_a, physical_gap/physical_a
    rng = random.Random(seed)
    start = time.monotonic()
    attempts = int(values.get("insertion_attempts_per_batch", 32))
    batches = int(values.get("insertion_batches_per_particle", 2048))
    sweeps = int(values.get("initialization_mc_sweeps", 100))
    stagnation = int(values.get("insertion_stagnation_attempts", 128))
    local_particles = int(values.get("insertion_local_particles", 16))
    local_sweeps = int(values.get("insertion_local_sweeps", 8))
    if attempts <= 0 or batches <= 0 or sweeps < 0:
        raise ValueError("Insertion counts must be positive; MC sweeps nonnegative.")
    if stagnation <= 0 or local_particles <= 0 or local_sweeps <= 0:
        raise ValueError('Stagnation attempts, local particles and local sweeps must be positive.')
    insertion_stats = {}
    bodies, accepted, attempted, proposals = random_insertion(
        rng, count, length, a, c, gap, attempts, batches,
        stagnation, local_particles, local_sweeps, insertion_stats)
    insertion_seconds = time.monotonic()-start
    if sweeps:
        print('Placement complete; starting %d final geometric MC sweeps.' % sweeps,
              file=sys.stderr, flush=True)
    final_start = time.monotonic()
    acc, att = randomize(rng, bodies, length, a, c, gap, sweeps, progress=True)
    final_seconds = time.monotonic()-final_start
    accepted += acc
    attempted += att
    print('Validating all particle pairs and periodic images independently of the cell list.',
          file=sys.stderr, flush=True)
    validation_start = time.monotonic()
    validate(bodies, length, a, c, gap)
    validation_seconds = time.monotonic()-validation_start
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("id", "x_m", "y_m", "z_m", "angle_x_deg", "angle_y_deg", "angle_z_deg"))
        for i, (position, q, _) in enumerate(bodies):
            writer.writerow((i,) + tuple(format(v*physical_a, ".17g") for v in position)
                            + tuple(format(v, ".17g") for v in euler_degrees(q)))
    orientation = [[sum(body[2][i]*body[2][j] for body in bodies)/count if count else 0.0
                    for j in range(3)] for i in range(3)]
    metadata = {
        "particle_count": count, "seed": seed,
        "initialization_method": "cell_list_random_insertion_with_stagnation_local_MC",
        "initialization_algorithm_version": 2,
        "equilibrated": False, "isotropy_guaranteed": False,
        "rotation_convention": "body-to-world Rz(angle_z) Ry(angle_y) Rx(angle_x); degrees; body thin axis z",
        "periodic_axes": ["x", "y", "z"], "initial_shear_phase": 0.0,
        "minimum_gap_m": physical_gap,
        "geometric_solid_volume_fraction": count*(4*math.pi/3)*physical_a**2*physical_c/physical_length**3,
        "orientation_second_moment": orientation,
        "random_insertion_particles": count, "random_proposals": proposals,
        "geometric_mc_accepted": accepted, "geometric_mc_attempted": attempted,
        "insertion_stagnation_attempts": stagnation,
        "insertion_local_particles": local_particles,
        "insertion_local_sweeps": local_sweeps,
        "local_relaxation": insertion_stats,
        "final_mc_sweeps": sweeps, "final_mc_accepted": acc, "final_mc_attempted": att,
        "insertion_seconds": insertion_seconds,
        "final_mc_seconds": final_seconds,
        "validation_seconds": validation_seconds,
        "elapsed_seconds": time.monotonic()-start,
        "validation": "all particle pairs and relevant xyz periodic images; PW contact and separating-plane clearance",
    }
    output.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Graphite simulation JSON config")
    parser.add_argument("--output", required=True, help="Output particle CSV")
    args = parser.parse_args()
    with open(args.config, encoding="utf-8") as stream:
        config = json.load(stream)
    generate(config, args.output)


if __name__ == "__main__":
    main()
