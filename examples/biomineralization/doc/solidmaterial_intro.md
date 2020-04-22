# Biomineralization - Solids

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
For further specialization, the overview over the material subfolder is split into solid and fluid, this description considering the solids.

## Solids in the folder `material`

As this example is about biomineralization involving many components with complex inteactions, some specific solid material files are necessary.
First, `material/components/biofilm.hh` defines the solid component biofilm, which plays an essential role in biomineralization by providing the essential catalytic activity for biomineralization.
Second, as biofilm grows or calcite precipitates, the volume fractions of those non-inert solids change, influencing the overall properties of the solids in the system.
The specific solidsystem `material/solidsystems/biominsolids.hh` gives the relations on how to
calculate the resulting average solid properties based on the solids volume fractions.


The subsequent documentation is structured as follows:

[[_TOC_]]


[@Span1996]: https://aip.scitation.org/doi/abs/10.1063/1.555991 "A new equation of state for carbon dioxide covering the fluid region from the triple-point temperature to 1100 K at pressures up to 800 MPa"
