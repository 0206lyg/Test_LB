#!/usr/bin/env python3
"""Generate random nonoverlapping oblate spheroids in a fully periodic box.

Python standard library only. Positions and orientations are proposed randomly.
A single progressive insertion algorithm interleaves candidate batches with
geometric Monte Carlo moves of already inserted particles. These are overlap-
rejecting coordinate proposals, with no physical time, dynamics or aging.
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
    ranges = [range(math.ceil((-cutoff-d)/length),
                    math.floor((cutoff-d)/length)+1) for d in dr]
    for image in itertools.product(*ranges):
        yield tuple(dr[i]+image[i]*length for i in range(3))


def fits(candidate, bodies, length, a, c, gap, skip=-1):
    position, _, normal = candidate
    for j, other in enumerate(bodies):
        if j == skip:
            continue
        for dr in image_displacements(position, other[0], length, 2*a+gap):
            if not separated(dr, normal, other[2], a, c, gap):
                return False
    return True


def random_candidate(rng, length, a, c, gap):
    return particle(tuple(rng.uniform(0, length) for _ in range(3)),
                    unit_quaternion(rng))


def random_insertion(rng, count, length, a, c, gap, attempts, batches):
    """One progressive random protocol; never substitute a crystalline seed."""
    bodies = []
    accepted = attempted = proposals = 0
    for index in range(count):
        inserted = False
        for _ in range(batches):
            for __ in range(attempts):
                proposals += 1
                candidate = random_candidate(rng, length, a, c, gap)
                if fits(candidate, bodies, length, a, c, gap):
                    bodies.append(candidate)
                    inserted = True
                    break
            # Geometric moves are part of every insertion batch, even successful
            # ones. No simulation time or forces are advanced here.
            acc, att = randomize(rng, bodies, length, a, c, gap, 1)
            accepted += acc
            attempted += att
            if inserted:
                break
        if not inserted:
            raise RuntimeError('Random preparation placed %d/%d particles; increase '
                               'particles.insertion_batches_per_particle explicitly.'
                               % (index, count))
        if (index + 1) % 12 == 0 or index + 1 >= count-12:
            print('Random placement: %d/%d particles' % (index+1, count),
                  file=sys.stderr, flush=True)
    return bodies, accepted, attempted, proposals


def randomize(rng, bodies, length, a, c, gap, sweeps):
    accepted, attempted = 0, 0
    for _ in range(sweeps):
        order = list(range(len(bodies)))
        rng.shuffle(order)
        for index in order:
            attempted += 1
            old = bodies[index]
            # Occasional whole-box relocations use independent uniform SO(3).
            if rng.random() < .1:
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
            if fits(candidate, bodies, length, a, c, gap, index):
                bodies[index] = candidate
                accepted += 1
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
    if not (physical_length > 0 and physical_a >= physical_c > 0 and physical_gap >= 0):
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
    if attempts <= 0 or batches <= 0 or sweeps < 0:
        raise ValueError("Insertion counts must be positive; MC sweeps nonnegative.")
    bodies, accepted, attempted, proposals = random_insertion(
        rng, count, length, a, c, gap, attempts, batches)
    acc, att = randomize(rng, bodies, length, a, c, gap, sweeps)
    accepted += acc
    attempted += att
    validate(bodies, length, a, c, gap)
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
        "particle_count": count, "seed": seed, "initialization_method": "progressive_random_insertion_with_geometric_MC",
        "equilibrated": False, "isotropy_guaranteed": False,
        "rotation_convention": "body-to-world Rz(angle_z) Ry(angle_y) Rx(angle_x); degrees; body thin axis z",
        "periodic_axes": ["x", "y", "z"], "initial_shear_phase": 0.0,
        "minimum_gap_m": physical_gap,
        "geometric_solid_volume_fraction": count*(4*math.pi/3)*physical_a**2*physical_c/physical_length**3,
        "orientation_second_moment": orientation,
        "random_insertion_particles": count, "random_proposals": proposals,
        "geometric_mc_accepted": accepted, "geometric_mc_attempted": attempted,
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
