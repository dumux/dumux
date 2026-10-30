# Part 2: Postprocessing/Visualization of a pore-network model

After compiling and running the code, both the quasi-static two-phase sweep's invasion/trapping
state and each single-phase upscaling solve's pressure field can be visualized with the open
source visualization application ParaView.

Three separate `.pvd`/`.vtp` sequences are written, one per stage/fluid:
```bash
paraview build-cmake/examples/porenetwork_upscaling_twophase/static.pvd &            # stage 1: invasion/trapping
paraview build-cmake/examples/porenetwork_upscaling_twophase/upscaling_liquid.pvd &  # stage 2: wetting-phase solve
paraview build-cmake/examples/porenetwork_upscaling_twophase/upscaling_gas.pvd &     # stage 2: non-wetting-phase solve
```

For general instructions on how to represent pore bodies as spheres/boxes and pore throats as
lines/tubes (with radii scaled by `poreInscribedRadius`/`throatInscribedRadius`), see [Part 4 of
the pore-network upscaling example](../../porenetwork_upscaling/doc/postprocessing.md) -- the
same technique applies here; only the fields of interest differ, as described below.

## Fields specific to this example

**`static.pvd`** (stage 1) writes the following per-timestep fields, which are the ones most
relevant for inspecting drainage/imbibition/trapping:

* **`invaded`** (cell data, per throat) -- `1` if the throat currently carries non-wetting phase, `0` otherwise.
* **`trapped`** (cell data, per throat) -- `1` if the throat is invaded but has lost its connected path back to the inlet (i.e. is part of a non-wetting-phase cluster whose local capillary pressure is frozen). Only ever `1` where `invaded` is also `1`.
* **`pc`** / **`sw`** (point data, per pore) -- the local capillary pressure and wetting-phase saturation of each pore body.
* **`pcEntry`** / **`pcSnapOff`** (cell data, per throat) -- the precomputed entry and snap-off threshold capillary pressures.

A useful way to see trapping directly in ParaView: color the throats (e.g. as tubes, see the
linked instructions above) by `trapped`, and step through the time series with the animation
controls -- invaded-but-trapped clusters (value `1`) appear as isolated red pockets once the
imbibition leg starts retreating.

**`upscaling_liquid.pvd` / `upscaling_gas.pvd`** (stage 2, one sequence per fluid, one frame per
saturation step): the standard single-phase pore-network output fields (pressure at every pore,
flow rate through every throat), computed on that step's frozen, imported conductance field.
Comparing the pressure field's spatial pattern between early (nearly all-conductive) and late
(mostly zero-conductance) steps of a given fluid's sequence illustrates why its apparent
permeability drops as that fluid's connected pathway shrinks.
