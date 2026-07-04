#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Convert a VTK dataset into an embeddable interactive 3D scene.

Reads a VTK dataset, optionally warps it to visualize deformation, and writes
it in a form that can be embedded into the Doxygen documentation. The output
format is chosen from the OUTPUT extension:

  * ``.js``   the (warped) surface geometry embedded in a small JavaScript file
              that registers it under ``window.dumuxScenes`` when loaded with a
              plain ``<script src>``. This is the default used by the
              documentation: the scenes are rendered in the browser by
              ``doc/doxygen/dumux-scene.js`` using vtk.js. Because the geometry
              is loaded via a script tag (not fetched), the scene also works
              from ``file://`` without a web server. Requires only PyVista.

  * ``.vtp``  the same surface as a plain VTK PolyData file (loaded by fetching
              the URL, so it needs the docs to be served over http(s)).
              Requires only PyVista.

  * ``.html`` a single, self-contained interactive HTML file that embeds both
              the geometry and the vtk.js viewer (via PyVista's ``export_html``).
              Needs no server and no network at view time, but is large (~1 MB).
              Additionally requires ``trame``/``trame-vtk``/``trame-vuetify``
              and ``nest_asyncio`` (only for generation, not for viewing).

Usage:
  python3 vtk_to_scene.py INPUT OUTPUT[.js|.vtp|.html] [options]

  INPUT   a .vtu/.vtp/.vtk file or a .pvd time-series collection
  OUTPUT  the .js, .vtp or .html file to write

Example (warp a plate by its vertical deformation ``w`` and color by it):
  python3 vtk_to_scene.py test.pvd scene_plate.js --scalar w --warp-scalar w
"""

import argparse
import base64
import json
import os
import sys


def _read_dataset(path, time):
    """Read INPUT into a single PyVista dataset, selecting a time step for .pvd."""
    import pyvista as pv

    if path.lower().endswith(".pvd"):
        reader = pv.get_reader(path)
        times = list(reader.time_values)
        if times:
            if time == "last":
                reader.set_active_time_value(times[-1])
            elif time == "first":
                reader.set_active_time_value(times[0])
            else:
                reader.set_active_time_value(times[int(time)])
        dataset = reader.read()
    else:
        dataset = pv.read(path)

    # Flatten a MultiBlock (e.g. from a .pvd) into one dataset
    if isinstance(dataset, pv.MultiBlock):
        dataset = dataset.combine()

    return dataset


def _auto_warp_factor(mesh, scalar, relative):
    """Scale so the maximum absolute warp equals ``relative`` * bounding-box diagonal."""
    import numpy as np

    values = np.abs(np.asarray(mesh[scalar]))
    peak = float(values.max()) if values.size else 0.0
    if peak <= 0.0:
        return 1.0

    bounds = mesh.bounds
    diag = float(
        np.sqrt(
            (bounds[1] - bounds[0]) ** 2
            + (bounds[3] - bounds[2]) ** 2
            + (bounds[5] - bounds[4]) ** 2
        )
    )
    if diag <= 0.0:
        diag = 1.0
    return relative * diag / peak


def _surface_with_scalar(mesh, scalar):
    """Extract the surface of ``mesh`` with ``scalar`` set as the active field."""
    surface = mesh.extract_surface()
    if scalar in surface.point_data:
        surface.set_active_scalars(scalar)
    return surface


def _write_vtp(mesh, scalar, output_path):
    """Write the surface of ``mesh`` as PolyData with ``scalar`` as active field."""
    _surface_with_scalar(mesh, scalar).save(output_path)


def _write_js(mesh, scalar, output_path):
    """Embed the surface PolyData in a JS file registering it on ``window.dumuxScenes``.

    The geometry is stored as base64-encoded VTK PolyData (.vtp) bytes. Loading
    the file with a plain ``<script src>`` avoids any fetch, so the scene renders
    from ``file://`` as well as over http(s).
    """
    surface = _surface_with_scalar(mesh, scalar)
    tmp_vtp = output_path + ".tmp.vtp"
    surface.save(tmp_vtp)
    try:
        with open(tmp_vtp, "rb") as vtp_file:
            encoded = base64.b64encode(vtp_file.read()).decode("ascii")
    finally:
        os.remove(tmp_vtp)

    key = os.path.splitext(os.path.basename(output_path))[0]
    with open(output_path, "w", encoding="utf-8") as js_file:
        js_file.write("window.dumuxScenes = window.dumuxScenes || {};\n")
        js_file.write(f"window.dumuxScenes[{json.dumps(key)}] = {json.dumps(encoded)};\n")


def _write_html(mesh, scalar, component, cmap, title, camera, zoom, show_edges,
                window_size, output_path):
    """Write a self-contained interactive HTML scene (requires the trame stack)."""
    import pyvista as pv

    plotter = pv.Plotter(off_screen=True, window_size=list(window_size))
    plotter.add_mesh(
        mesh,
        scalars=scalar,
        component=component,
        cmap=cmap,
        show_edges=show_edges,
        scalar_bar_args={"title": title},
    )
    views = {
        "iso": plotter.view_isometric,
        "xy": plotter.view_xy,
        "xz": plotter.view_xz,
        "yz": plotter.view_yz,
    }
    views.get(camera, plotter.view_isometric)()
    if zoom and zoom != 1.0:
        plotter.camera.zoom(zoom)
    plotter.export_html(output_path)
    plotter.close()


def export_scene(
    input_path,
    output_path,
    scalar=None,
    component=None,
    warp_scalar=None,
    warp_vector=None,
    warp_factor="auto",
    warp_relative=0.15,
    warp_normal=(0.0, 0.0, 1.0),
    cmap="RdBu_r",
    title=None,
    camera="iso",
    zoom=1.0,
    show_edges=False,
    window_size=(800, 600),
    time="last",
):
    """Read ``input_path``, optionally warp, and write ``output_path`` (.vtp or .html)."""
    mesh = _read_dataset(input_path, time)

    # Default the coloring scalar to the first available point field
    if scalar is None:
        point_scalars = list(mesh.point_data.keys())
        if not point_scalars:
            raise SystemExit(f"No point data found in '{input_path}' to color by.")
        scalar = point_scalars[0]
    if scalar not in mesh.point_data and scalar not in mesh.cell_data:
        raise SystemExit(f"Field '{scalar}' not found in '{input_path}'.")

    if title is None:
        title = scalar

    # Optionally warp the geometry to visualize deformation
    if warp_vector is not None:
        factor = 1.0 if warp_factor == "auto" else float(warp_factor)
        mesh = mesh.warp_by_vector(warp_vector, factor=factor)
    elif warp_scalar is not None:
        if warp_factor == "auto":
            factor = _auto_warp_factor(mesh, warp_scalar, warp_relative)
        else:
            factor = float(warp_factor)
        mesh = mesh.warp_by_scalar(warp_scalar, factor=factor, normal=list(warp_normal))

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    lower = output_path.lower()
    if lower.endswith(".html"):
        _write_html(mesh, scalar, component, cmap, title, camera, zoom, show_edges,
                    window_size, output_path)
    elif lower.endswith(".vtp"):
        _write_vtp(mesh, scalar, output_path)
    else:
        _write_js(mesh, scalar, output_path)
    return output_path


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        description=(
            "Convert a VTK dataset (.vtu/.vtp/.vtk/.pvd) into an embeddable "
            "interactive 3D scene (.js or .vtp data for vtk.js, or a "
            "self-contained .html)."
        )
    )
    parser.add_argument("input", help="input .vtu/.vtp/.vtk file or .pvd collection")
    parser.add_argument("output", help="output .js, .vtp or .html file")
    parser.add_argument("--scalar", help="field to color by (default: first point field)")
    parser.add_argument(
        "--component", type=int, default=None,
        help="component index when coloring by a vector field (html output only)",
    )
    parser.add_argument("--warp-scalar", help="warp geometry along a normal by this scalar field")
    parser.add_argument("--warp-vector", help="warp geometry by this vector field")
    parser.add_argument(
        "--warp-factor", default="auto",
        help="warp scale factor, or 'auto' to scale by the bounding box (default: auto)",
    )
    parser.add_argument(
        "--warp-relative", type=float, default=0.15,
        help="target max warp as a fraction of the bounding-box diagonal for --warp-factor auto",
    )
    parser.add_argument(
        "--warp-normal", type=float, nargs=3, default=[0.0, 0.0, 1.0],
        metavar=("X", "Y", "Z"), help="direction for scalar warping (default: 0 0 1)",
    )
    parser.add_argument("--cmap", default="RdBu_r", help="color map name for html output (default: RdBu_r)")
    parser.add_argument("--title", help="scalar bar title (default: the scalar name)")
    parser.add_argument(
        "--camera", default="iso", choices=["iso", "xy", "xz", "yz"],
        help="initial camera orientation for html output (default: iso)",
    )
    parser.add_argument("--zoom", type=float, default=1.0, help="camera zoom factor for html output")
    parser.add_argument("--show-edges", action="store_true", help="draw mesh edges (html output)")
    parser.add_argument(
        "--window-size", type=int, nargs=2, default=[800, 600],
        metavar=("W", "H"), help="viewer size in pixels for html output (default: 800 600)",
    )
    parser.add_argument(
        "--time", default="last",
        help="time step for .pvd input: 'last', 'first', or an index (default: last)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    try:
        import pyvista  # noqa: F401
    except ImportError:
        sys.stderr.write("PyVista is required. Install it with 'pip install pyvista'.\n")
        return 2

    if args.output.lower().endswith(".html"):
        try:
            import trame_vtk  # noqa: F401
            import nest_asyncio  # noqa: F401
        except ImportError:
            sys.stderr.write(
                "Self-contained HTML export requires 'trame', 'trame-vtk', "
                "'trame-vuetify' and 'nest_asyncio'. Install them with:\n"
                "  pip install trame trame-vtk trame-vuetify nest_asyncio\n"
                "Alternatively, write a .vtp file (needs only PyVista).\n"
            )
            return 2

    output = export_scene(
        args.input,
        args.output,
        scalar=args.scalar,
        component=args.component,
        warp_scalar=args.warp_scalar,
        warp_vector=args.warp_vector,
        warp_factor=args.warp_factor,
        warp_relative=args.warp_relative,
        warp_normal=tuple(args.warp_normal),
        cmap=args.cmap,
        title=args.title,
        camera=args.camera,
        zoom=args.zoom,
        show_edges=args.show_edges,
        window_size=tuple(args.window_size),
        time=args.time,
    )
    print(f"Wrote scene: {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
