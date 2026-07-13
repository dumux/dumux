#!/usr/bin/env pvpython
"""Render phi=0 (meniscus) isosurfaces at several times in one 3D scene.

Must be run with ParaView's pvpython (plain python3 has no vtk/paraview
bindings in this environment):

    /Applications/ParaView-6.0.1.app/Contents/bin/pvpython \
        plot_levelset_3d.py <pvd> --out snapshot.png

The 3D stefantube runs simulate only one quarter of the square duct
(y,z >= 0, using mirror symmetry). This script reflects the data across
y=0 and z=0 to reconstruct the full physical duct before rendering, and
overlays a wireframe outline plus translucent wall planes for context.
Each requested time is baked into a static isosurface (via ForceTime) so
all snapshots can be shown together in a single non-animated image,
colored early (blue) to late (red) and increasingly opaque over time.
"""

import argparse
from pathlib import Path


def select_snapshots(times, count):
    if len(times) <= count:
        return list(times)
    indices = sorted({round(i*(len(times) - 1)/(count - 1)) for i in range(count)})
    return [times[i] for i in indices]


def lerp(a, b, t):
    return tuple(ai + (bi - ai)*t for ai, bi in zip(a, b))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pvd")
    parser.add_argument("--out", default="stefantube_levelset_3d.png")
    parser.add_argument("--count", type=int, default=6, help="number of time snapshots to overlay")
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=1350)
    parser.add_argument("--azimuth", type=float, default=-125.0)
    parser.add_argument("--elevation", type=float, default=20.0)
    parser.add_argument("--zoom", type=float, default=1.0)
    parser.add_argument("--early-color", type=float, nargs=3, default=[0.16, 0.44, 0.86])
    parser.add_argument("--late-color", type=float, nargs=3, default=[0.85, 0.22, 0.12])
    args = parser.parse_args()

    from paraview.simple import (
        PVDReader, Reflect, ForceTime, Contour, Outline, Plane,
        Show, GetActiveViewOrCreate, Render, SaveScreenshot,
        ResetCamera, Text,
    )

    reader = PVDReader(FileName=str(Path(args.pvd).resolve()))
    reader.UpdatePipeline()
    times = list(reader.TimestepValues)
    snaps = select_snapshots(times, args.count)

    reflY = Reflect(Input=reader)
    reflY.ReflectionPlane.Normal = [0, 1, 0]
    reflY.ReflectionPlane.Origin = [0, 0, 0]
    reflZ = Reflect(Input=reflY)
    reflZ.ReflectionPlane.Normal = [0, 0, 1]
    reflZ.ReflectionPlane.Origin = [0, 0, 0]
    reflZ.UpdatePipeline()

    xmin, xmax, ymin, ymax, zmin, zmax = reflZ.GetDataInformation().GetBounds()

    view = GetActiveViewOrCreate('RenderView')
    view.ViewSize = [args.width, args.height]
    view.UseColorPaletteForBackground = 0
    view.Background = [1.0, 1.0, 1.0]
    view.Background2 = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 1
    view.DepthPeeling = 1
    view.MaximumNumberOfPeels = 200

    # ---- duct outline (wireframe box of the full reconstructed duct) ----
    outline = Outline(Input=reflZ)
    outlineDisp = Show(outline, view)
    outlineDisp.SetRepresentationType('Wireframe')
    outlineDisp.AmbientColor = [0.05, 0.05, 0.05]
    outlineDisp.DiffuseColor = [0.05, 0.05, 0.05]
    outlineDisp.LineWidth = 2.5

    # ---- translucent duct walls (glass-tube look) ----
    wall_defs = [
        ([xmin, ymax, zmin], [xmax, ymax, zmin], [xmin, ymax, zmax]),  # y = ymax
        ([xmin, ymin, zmin], [xmax, ymin, zmin], [xmin, ymin, zmax]),  # y = ymin
        ([xmin, ymin, zmax], [xmax, ymin, zmax], [xmin, ymax, zmax]),  # z = zmax
        ([xmin, ymin, zmin], [xmax, ymin, zmin], [xmin, ymax, zmin]),  # z = zmin
    ]
    for origin, p1, p2 in wall_defs:
        wall = Plane()
        wall.Origin = origin
        wall.Point1 = p1
        wall.Point2 = p2
        wallDisp = Show(wall, view)
        wallDisp.AmbientColor = [0.72, 0.80, 0.90]
        wallDisp.DiffuseColor = [0.72, 0.80, 0.90]
        wallDisp.Opacity = 0.12

    # ---- meniscus isosurfaces, one per snapshot time ----
    # Contour the RAW quarter-domain data first, then reflect the resulting
    # surface -- reflecting first and contouring the stitched volume instead
    # left a visible ridge along the mirror seam (marching cubes sees a tiny
    # mismatch between the original and mirrored copies right at y=0/z=0);
    # mirroring an already-extracted surface is exact and seam-free.
    n = max(1, len(snaps) - 1)
    for i, t in enumerate(snaps):
        frac = i/n
        forced = ForceTime(Input=reader)
        forced.ForcedTime = t
        contour = Contour(Input=forced)
        contour.ContourBy = ['POINTS', 'phi']
        contour.Isosurfaces = [0.0]
        contour.UpdatePipeline()

        surfReflY = Reflect(Input=contour)
        surfReflY.ReflectionPlane.Normal = [0, 1, 0]
        surfReflY.ReflectionPlane.Origin = [0, 0, 0]
        surfReflZ = Reflect(Input=surfReflY)
        surfReflZ.ReflectionPlane.Normal = [0, 0, 1]
        surfReflZ.ReflectionPlane.Origin = [0, 0, 0]
        surfReflZ.UpdatePipeline()

        disp = Show(surfReflZ, view)
        color = lerp(tuple(args.early_color), tuple(args.late_color), frac)
        disp.AmbientColor = list(color)
        disp.DiffuseColor = list(color)
        disp.Opacity = 0.30 + 0.65*frac
        disp.Specular = 0.15

        label = Text()
        label.Text = f"t = {t:.4g}"
        labelDisp = Show(label, view)
        labelDisp.WindowLocation = 'Any Location'
        labelDisp.Position = [0.015, 0.94 - 0.045*i]
        labelDisp.Color = list(color)
        labelDisp.FontSize = 16
        labelDisp.Bold = 1

    ResetCamera(view)
    cam = view.GetActiveCamera()
    cam.Azimuth(args.azimuth)
    cam.Elevation(args.elevation)
    cam.Zoom(args.zoom)

    Render(view)
    SaveScreenshot(str(Path(args.out).resolve()), view, ImageResolution=[args.width, args.height])
    print(f"Wrote {Path(args.out).resolve()}")


if __name__ == "__main__":
    main()
