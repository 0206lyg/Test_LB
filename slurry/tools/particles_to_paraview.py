#!/usr/bin/env python3
"""Export saved graphite poses as ellipsoid surfaces for ParaView (stdlib only).

Usage: python particles_to_paraview.py
Reads particles.csv and effective_config.json in the current result directory
(or pass another result directory as the first argument).
By default, export the closest saved frame in each rounded 0.05-strain bin:
0, 0.05, 0.10, ... . Empty bins are skipped; ties keep the earlier frame.
Strain = time_s * flow.shear_rate_s_inv, as in the constant-shear solver.
Use --strain-interval 0 to export every saved frame.
Open paraview_particles/particles.pvd and optionally box.vtp in ParaView.
Coordinates are SI metres; PVD times are the recorded physical time_s values.
Euler angles are radians with R = Rz @ Ry @ Rx. Colours particle_id and
normal_y_squared identify particles and their short-axis orientation to y.

Only recorded particle states are exported: no added/interpolated states, and
no pore/porosity, velocity or pressure field is reconstructed. This works with
old particle CSV output even when fluid VTK output was disabled. Whole ellipsoids
are drawn at their recorded periodic centres; portions outside x/z box faces
are intentionally not clipped and periodic copies are not generated.
An incomplete final frame is skipped with a warning; invalid earlier data fail.
"""

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET


POSE_KEYS = ("x_m", "y_m", "z_m", "angle_x_rad", "angle_y_rad", "angle_z_rad")
REQUIRED = ("step", "time_s", "id") + POSE_KEYS


def rotation(ax, ay, az):
    sx, cx, sy, cy, sz, cz = math.sin(ax), math.cos(ax), math.sin(ay), math.cos(ay), math.sin(az), math.cos(az)
    return ((cz*cy, cz*sy*sx-sz*cx, cz*sy*cx+sz*sx),
            (sz*cy, sz*sy*sx+cz*cx, sz*sy*cx-cz*sx),
            (-sy, cy*sx, cy*cx))


def template(a, c, longitude, latitude):
    vertices = [(0., 0., c)]
    for j in range(1, latitude):
        theta = math.pi*j/latitude
        for i in range(longitude):
            phi = 2*math.pi*i/longitude
            vertices.append((a*math.sin(theta)*math.cos(phi), a*math.sin(theta)*math.sin(phi), c*math.cos(theta)))
    vertices.append((0., 0., -c))
    faces = [(0, 1+i, 1+(i+1) % longitude) for i in range(longitude)]
    for j in range(latitude-2):
        for i in range(longitude):
            u, v = 1+j*longitude+i, 1+j*longitude+(i+1) % longitude
            faces.extend(((u, u+longitude, v), (v, u+longitude, v+longitude)))
    bottom, last_ring = len(vertices)-1, 1+(latitude-2)*longitude
    faces.extend((bottom, last_ring+(i+1) % longitude, last_ring+i) for i in range(longitude))
    return vertices, faces


def frames(path, count):
    """Stream complete frames, validating IDs, times and ordering."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or not set(REQUIRED).issubset(reader.fieldnames):
            raise ValueError("particles.csv lacks required pose columns")
        current, rows, seen, last_time = None, [], set(), -1.
        iterator = iter(reader)
        raw = next(iterator, None)
        while raw is not None:
            following = next(iterator, None)
            # A process stopped while writing can leave a truncated final row.
            if any(raw.get(k) in (None, "") for k in REQUIRED):
                if following is not None:
                    raise ValueError("Incomplete row before the end of particles.csv")
                try:
                    partial_step = int(raw["step"])
                except (TypeError, ValueError, KeyError):
                    partial_step = None
                if rows and len(rows) == count and partial_step is not None and partial_step > current:
                    yield current, rows[0][1], rows
                print("Warning: incomplete trailing frame ignored.", file=sys.stderr)
                return
            try:
                step, ident = int(raw["step"]), int(raw["id"])
                time = float(raw["time_s"])
                pose = tuple(float(raw[k]) for k in POSE_KEYS)
            except (ValueError, TypeError) as exc:
                raise ValueError("Invalid numeric value in particles.csv") from exc
            if step < 0 or time < 0 or not all(math.isfinite(v) for v in (time,) + pose):
                raise ValueError("Nonfinite or negative step/time in particles.csv")
            if current is not None and step != current:
                if step <= current or time <= last_time:
                    raise ValueError("Frame steps/times must increase strictly")
                if len(rows) != count:
                    raise ValueError(f"Incomplete middle frame at step {current}: {len(rows)}/{count} particles")
                yield current, rows[0][1], rows
                rows, seen = [], set()
            if rows and time != rows[0][1]:
                raise ValueError(f"Inconsistent time_s within step {step}")
            if ident not in range(count) or ident in seen:
                raise ValueError(f"Invalid or duplicate particle ID {ident} at step {step}")
            current, last_time = step, time
            seen.add(ident)
            rows.append((ident, time, pose))
            raw = following
        if rows:
            if len(rows) != count:
                print(f"Warning: incomplete trailing frame {current} ignored ({len(rows)}/{count} particles).", file=sys.stderr)
            else:
                yield current, rows[0][1], rows


def strain_frames(saved_frames, shear_rate, interval):
    """Keep the closest recorded frame per rounded strain bin, without interpolation."""
    best, current_bin, best_error = None, None, math.inf
    for frame in saved_frames:
        scaled = frame[1] * shear_rate / interval
        if not math.isfinite(scaled):
            raise ValueError("Nonfinite strain/interval in particle data")
        target_bin = math.floor(scaled + 0.5)
        error = abs(scaled - target_bin)
        if target_bin != current_bin:
            if best is not None:
                yield best
            best, current_bin, best_error = frame, target_bin, error
        elif error < best_error:
            best, best_error = frame, error
    if best is not None:
        yield best


def data_array(parent, name, values, kind="Float64", components=1):
    attrs = {"type": kind, "format": "ascii", "NumberOfComponents": str(components)}
    if name:
        attrs["Name"] = name
    element = ET.SubElement(parent, "DataArray", attrs)
    element.text = " ".join(str(v) if kind.startswith("Int") else format(v, ".17g") for v in values)


def write_vtp(path, points, cells, ids=None, directions=None, time=None, step=None, lines=False):
    root = ET.Element("VTKFile", type="PolyData", version="0.1", byte_order="LittleEndian")
    polydata = ET.SubElement(root, "PolyData")
    if time is not None:
        fields = ET.SubElement(polydata, "FieldData")
        data_array(fields, "TimeValue", [time])
        data_array(fields, "step", [step], "Int64")
        for array in fields:
            array.set("NumberOfTuples", "1")
    piece = ET.SubElement(polydata, "Piece", NumberOfPoints=str(len(points)), NumberOfVerts="0",
                          NumberOfLines=str(len(cells) if lines else 0), NumberOfStrips="0",
                          NumberOfPolys=str(0 if lines else len(cells)))
    for tag, size in (("PointData", len(points)), ("CellData", len(cells))):
        parent = ET.SubElement(piece, tag)
        if ids is not None:
            copies = size // len(ids)
            data_array(parent, "particle_id", (ident for ident in ids for _ in range(copies)), "Int32")
            data_array(parent, "normal_y_squared", (value for value in directions for _ in range(copies)))
    data_array(ET.SubElement(piece, "Points"), None, (x for p in points for x in p), components=3)
    parent = ET.SubElement(piece, "Lines" if lines else "Polys")
    data_array(parent, "connectivity", (i for cell in cells for i in cell), "Int64")
    data_array(parent, "offsets", (len(cells[0])*(i+1) for i in range(len(cells))), "Int64")
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def export(args):
    source = args.result_dir.resolve()
    config = json.loads((source / "effective_config.json").read_text(encoding="utf-8"))
    count = config["particles"]["count"]
    a, c = float(config["particles"]["diameter_m"])/2, float(config["particles"]["thickness_m"])/2
    length = float(config["geometry"]["box_length_m"])
    if not isinstance(count, int) or count <= 0 or not all(math.isfinite(v) and v > 0 for v in (a, c, length)):
        raise ValueError("Config requires positive dimensions and a positive integer particle count")
    saved_frames = frames(source / "particles.csv", count)
    if args.strain_interval > 0:
        shear_rate = float(config["flow"]["shear_rate_s_inv"])
        if not math.isfinite(shear_rate) or shear_rate <= 0:
            raise ValueError("Config requires a finite positive flow.shear_rate_s_inv")
        saved_frames = strain_frames(saved_frames, shear_rate, args.strain_interval)
    output = (args.output if args.output else source / "paraview_particles").resolve()
    if output.exists() and not args.overwrite:
        raise ValueError(f"Output already exists: {output}; use --overwrite to replace it")
    if output.exists() and not (output / "CONVERTER_OUTPUT.txt").is_file():
        raise ValueError("Refusing to replace a directory not created by this converter")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".particle_export_", dir=output.parent))
    try:
        vertices, faces = template(a, c, args.longitude, args.latitude)
        master = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
        collection = ET.SubElement(master, "Collection")
        total = 0
        for step, time, rows in saved_frames:
            points, cells, ids, directions = [], [], [], []
            for ident, _, pose in sorted(rows):
                center, angles = pose[:3], pose[3:]
                matrix = rotation(*angles)
                offset = len(points)
                points.extend(tuple(center[j]+sum(matrix[j][k]*v[k] for k in range(3)) for j in range(3)) for v in vertices)
                cells.extend(tuple(offset+i for i in face) for face in faces)
                ids.append(ident)
                directions.append(matrix[1][2]**2)
            name = f"particles_{step:012d}.vtp"
            write_vtp(staging / name, points, cells, ids, directions, time, step)
            ET.SubElement(collection, "DataSet", timestep=format(time, ".17g"), group="", part="0", file=name)
            total += 1
        if not total:
            raise ValueError("No complete particle frames found")
        corners = [(x*length, y*length, z*length) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
        edges = [(i, i ^ bit) for i in range(8) for bit in (1, 2, 4) if i < (i ^ bit)]
        write_vtp(staging / "box.vtp", corners, edges, lines=True)
        ET.ElementTree(master).write(staging / "particles.pvd", encoding="utf-8", xml_declaration=True)
        (staging / "CONVERTER_OUTPUT.txt").write_text(
            __doc__ + f"\nstrain_interval = {args.strain_interval:g} (0 = all saved frames)\n", encoding="utf-8")
        if output.exists():
            shutil.rmtree(output)
        staging.rename(output)
        print(f"Exported {total} saved frames, {count} particles/frame, SI metres.\nOpen: {output / 'particles.pvd'}\nBox:  {output / 'box.vtp'}")
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("result_dir", type=Path, nargs="?", default=Path("."))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--strain-interval", type=float, default=0.05,
                        help="Strain spacing (default: 0.05); 0 exports all saved frames")
    parser.add_argument("--longitude", type=int, default=24, help="Surface angular resolution (default: 24)")
    parser.add_argument("--latitude", type=int, default=12, help="Polar subdivisions (default: 12)")
    args = parser.parse_args()
    if not math.isfinite(args.strain_interval) or args.strain_interval < 0:
        parser.error("strain-interval must be finite and nonnegative")
    if args.longitude < 4 or args.latitude < 2:
        parser.error("longitude must be >=4 and latitude >=2")
    try:
        export(args)
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
