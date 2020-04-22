# Biomineralization - Fluids

The main idea of biomineralization revolves around
biologically-induced mineral precipitation.
I our example, it is about precipitation of calcite
due to urea hydrolysis.
The resulting overall reaction equation is:

```math
\mathrm{CO(NH_2)_2} + 2\mathrm{H_2O} + \mathrm{Ca^{2+}} \xrightarrow{urease}
2\mathrm{NH_{4}^{+}} + \mathrm{CaCO_{3}} \downarrow.
```

To describe the processes, the minimum set of components is:
water (w), dissolved inorganic carbon (C<sub>tot</sub>),
sodium (Na), chloride (Cl), calcium (Ca), urea (u), glucose as a substrate (s), oxygen (O<sub>2</sub>), and suspended biomass (b),
as well as the solid components biofilm and calcite.
Thus, there are many components involved and keeping track of the components and their interactions,
necessitates a sophisticated handling.
This is why the material folder in this example is both separated from the regular material folder in dumux/dumux/material
and also documented separately.
For further specialization, the overview over the material subfolder is split into solid and fluid, this description considering the fluids.

## Fluids in the folder `material`

As this example is about biomineralization involving many components with complex inteactions, some specific fluid material files are necessary.
A CO_2-Table file provides tabulated CO_2 properties according to @Span1996 in `material/co2tableslaboratory.hh`
In the component subfolder, `material/components/suspendedbiomass.hh` defines the component suspended biomass, which is the mobile form of biomass being transported suspended in the aqueous fluid phase.
In the fluidsystem subfolder, the biomineralization fluidsystem `material/fluidsystems/biominsimplechemistry.hh` as well as
the complex salinity brine adapter `material/fluidsystems/icpcomplexsalinitybrine.hh` can be found.
The biomineralization fluidsystem `material/fluidsystems/biominsimplechemistry.hh`
contains the fluid related properties and the interactions of the components in the fluids.
The complex salinity brine adapter `material/fluidsystems/icpcomplexsalinitybrine.hh`
adapts the brine fluidsystem (dumux/dumux/material/fluidsystems/brine.hh) expecting a single non-aqueous component defining salinity to be reused for biomineralization with three ions (calcium, sodium, chloride) being assumed to contribute to salinity.


The subsequent documentation is structured as follows:

[[_TOC_]]


[@Span1996]: https://aip.scitation.org/doi/abs/10.1063/1.555991 "A new equation of state for carbon dioxide covering the fluid region from the triple-point temperature to 1100 K at pressures up to 800 MPa"
