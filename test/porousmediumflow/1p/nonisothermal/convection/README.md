# 1D Heat Convection in a Porous Medium

This test verifies the single-phase, non-isothermal model (`OnePModel` with `NIModel`)
for advective heat transport in a fully saturated porous medium.

Water enters a 1D domain of $20\ \mathrm{m}$ at the left boundary with a fixed Darcy velocity
and a higher temperature than the initial one. The numerical temperature is compared to
the exported field `temperatureExact`: a sharp temperature front that moves with the
retarded front velocity, which accounts for the heat stored in the solid.

The test runs with the box and the cell-centered TPFA scheme. A restart test checks that
the box simulation can be continued from an intermediate VTK output file.

```bash
cd dumux/build-cmake/test/porousmediumflow/1p/nonisothermal/convection
make test_1pni_convection_tpfa
./test_1pni_convection_tpfa params.input
```
