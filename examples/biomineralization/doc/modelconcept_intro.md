# Model Concept for Biomineralization

The conceptual model for biomineralization follows the one presented by [@Ebigbo2012] and [@Hommel2015].
It accounts for two-phase multi-component reactive transport on the continuum scale, including biofilm and calcite as solid phases.
The reactions considered are pH-dependent dissociation reactions, microbial growth, and decay as well as microbially-catalyzed ureolysis and mass-transfer reactions between the different phases.
A mass transfer may occur between both fluid phases by the mutual dissolution of water and CO<sub>2</sub> in the gas or the aqueous phase. It may also occur between
the aqueous phase and the two *solid* phases, biofilm and calcite, denoted by subscripts (f) and (c) respectively, by attachment or
detachment of biomass and precipitation or dissolution of calcite.

The mobile components, denoted by superscripts $`\kappa`$, are water (w), dissolved inorganic carbon (C<sub>tot</sub>),
sodium (Na), chloride (Cl), calcium (Ca), urea (u), substrate (s), oxygen (O<sub>2</sub>), and suspended biomass (b).
A substrate is the carbon and energy source of the bacteria and O<sub>2</sub> is the electron acceptor.

Contrary to [@Ebigbo2012] and [@Hommel2015], for the sake of simplicity,
the dissociation processes are not modeled in this example.
It is assumed that the system has reached a steady state were the calcite precipitation rate is equal to the rate of ureolysis.
Thus, in this example, we consider only the "simplified chemistry case" in which it is assumed that the precipitation rate is equal to the ureolysis rate to simplify the chemical processes considered ($`r_\mathrm{prec} = r_\mathrm{urea}`$),
described in detail in  Chapter 6 of [@Hommel2016].

The documentation provided for the model concept is structured as follows:

[[_TOC_]]


## Primary variables and mass balance equations

The primary variables of this system of equations are chosen as the aqueous-phase pressure, mole fractions of the components in the water phase, and
for the solid phases biofilm and calcite, volume fractions.
However, the CO<sub>2</sub>-phase saturation is used as the primary variable instead of the mole fraction of total inorganic carbon in water
whenever both fluid phases are present within the same control volume ([@Class2002]).
All reactive and mass-transfer processes are incorporated in the mass balance equations for the components by component-specific source and sink terms:

```math
\sum\limits_{\alpha} \left[\frac{\partial}{\partial t}\left(\phi \rho_\mathrm{\alpha,\,mol} x^\kappa_\alpha S_\alpha \right) + \nabla\cdotp \left(\rho_\mathrm{\alpha,\,mol} x^\kappa_\alpha \mathbf{v}_\alpha \right) - \nabla\cdotp \left(\rho_\mathrm{\alpha,\,mol} \mathbf{D}^\kappa_\mathrm{pm,\alpha} \nabla x^\kappa_\alpha \right) \right] = q^\kappa,\:\alpha\in \mathrm{\{n;w\}} .
\tag{2}
```

However, all components except water, CO<sub>2</sub> , and O<sub>2</sub> are assumed to be restricted to the water phase.

The mass balances for the solid phases calcite and biofilm contain
only a storage and source term since they are immobile:

```math
\frac{\partial}{\partial t} \left(\phi_\lambda \rho_\lambda \right) = q^\lambda,\:\lambda\in \mathrm{\{c;f\}}.
\tag{3}
```

Here, $`\phi_\lambda`$ and $`\rho_\lambda`$ are volume fraction and mass density of
the solid phase $`\lambda`$, and $`q^{\lambda}`$ is the source term of phase $`\lambda`$ due to biochemical reactions.
The sources and sinks due to reactions $`q^{\kappa}$ and $q^{\lambda}`$ are specific to the components and are discussed in details in the subsequent section.

## Component-specific reactive source and sink terms

The source and sink terms account for the biogeochemical reactions occurring during MICP and the presence of CO<sub>2</sub>:
ureolysis, calcite precipitation, and dissolution, biomass growth under consumption of oxygen and substrate, biomass decay, as well as attachment and detachment of biomass.

## Water, sodium and chloride
Sodium and chloride do not participate in the reactions and water is the solvent and, thus, abundant, which is why its consumption by the hydrolysis of urea (Eq. [1](#mjx-eqn-eq:q_w_na_cl)) is considered negligible.
Thus, the reactive source terms for water, sodium, and chloride are zero.

## Urea and total nitrogen
The source term for ammonia/ammonium, $`q^\mathrm{N_{tot}}`$, and the sink term for urea $`q^\mathrm{u}`$ result from ureolysis (Eq. [1](#mjx-eqn-eq:q_ntot)).
For each mole of urea hydrolyzed, 2 moles of ammonia/ammonium are generated. The source terms are thus:

```math
q^\mathrm{u}=-r_\mathrm{urea}; \\
q^\mathrm{N_{tot}}=2r_\mathrm{urea},
```

where $`r_\mathrm{urea}`$ is the ureolysis rate calculated according to [@Lauchnor2015],
who investigated the influences of urea,
NH<sup>4+</sup>, and cell concentration, and pH of the medium on the ureolysis of whole cells of _S. pasteurii_:

```math
r_\mathrm{urea} = k_\mathrm{urease}
\:k_\mathrm{ub}\:\rho_\mathrm{f}\:\phi_\mathrm{f}
\:\frac{m^\mathrm{u}}{m^\mathrm{u}+K_\mathrm{u}}.
```

$`r_\mathrm{urea}`$ represents the revised rate of ureolysis according to [@Lauchnor2015],
$`k_\mathrm{urease}`$ the revised maximum activity of urease adapted from  [@Lauchnor2015],
$`\rho_\mathrm{f}`$ and $`\phi_\mathrm{f}`$ the density and volume fraction of biofilm respectively,
$`k_\mathrm{ub}`$ the mass ratio of urease to biofilm,
$`m^\mathrm{u}`$ the molality of urea calculated from the water phase composition,
and $`K_\mathrm{u}`$ is the half saturation constant for urea adapted from [@Lauchnor2015].

## Calcium and calcite
The source terms of calcium $`q^\mathrm{Ca}`$ and calcite $`q^\mathrm{c}`$
are determined by the rates of precipitation and dissolution.
When the aqueous phase is oversaturated with respect to calcite, it precipitates.
In the opposite case, calcite dissolves until the solution is saturated or all calcite is already dissolved:

```math
q^\mathrm{Ca} = r_\mathrm{diss} - r_\mathrm{prec}; \\
q^\mathrm{c} = - r_\mathrm{diss} + r_\mathrm{prec}.
```

In the presented simplified chemistry system, dissolution is neglected
and instead of calculating the complex geochemistry, e.g. dissociation of inorganic carbon into carbonate and bicarbonate,
it is assumed that the system has reached a steady state, were the precipitation rate is equal to the ureolysis rate and, thus, it yields:

```math
r_\mathrm{prec} = r_\mathrm{urea},
```
as long as there is calcium to form calcite.
Note: The "simplified chemistry case" ($`r_\mathrm{prec} = r_\mathrm{urea}`$),
is a simplifying assumption, its effects and validity are discussed in more detail in Chapter 6 of [@Hommel2016].

## Dissolved inorganic carbon
Dissolved inorganic carbon is generated by the hydrolysis of urea
as well as by the dissolution of calcite while it is consumed by the precipitation of calcite.
Thus, the source term of dissolved inorganic carbon $`q^\mathrm{C_{tot}}`$ results in:

```math
q^\mathrm{C_{tot}} =  r_\mathrm{urea} + r_\mathrm{diss} - r_\mathrm{prec},
```


## Suspended and attached biomass

The source and sink terms of suspended and attached biomass (biofilm), $`q^\mathrm{b}`$ and $`q^\mathrm{f}`$,
include four reaction rates each, corresponding to the biomass-related processes the model accounts for.
These processes are growth and decay increasing and decreasing the suspended or attached biomass as well as
attachment and detachment describing the transfer of biomass from the suspended to the attached state and vice versa:

```math
q^\mathrm{b} = \frac{r^\mathrm{b}_\mathrm{g} - r^\mathrm{b}_\mathrm{b} - r_\mathrm{a} + r_\mathrm{d}}{M^\mathrm{b}}; \\
q^\mathrm{f} = \frac{r^\mathrm{f}_\mathrm{g} - r^\mathrm{f}_\mathrm{b} + r_\mathrm{a} - r_\mathrm{d}}{M^\mathrm{f}}
```

where $`r^\mathrm{b}_\mathrm{g}`$ is the growth rate and $`r^\mathrm{b}_\mathrm{b}`$ the decay rate of suspended biomass, $`r_\mathrm{a}`$ the attachment rate, $`r_\mathrm{d}`$ the detachment rate.
Accordingly, $`r^\mathrm{f}_\mathrm{g}`$ and $`r^\mathrm{f}_\mathrm{b}`$ are the growth and decay of biofilm.
All rates influencing both attached and suspended biomass are assumed to be of a first-order type.

The growth rates of suspended and attached biomass are as given below,
with the specific growth rate calculated using double Monod kinetics to reproduce the dependence of
the microbial growth on both substrate and oxygen.

```math
r_\mathrm{g}^\mathrm{b} = \mu_\mathrm{g} C_\mathrm{w}^\mathrm{b}S_\mathrm{w}\phi; \\
  r_\mathrm{g}^\mathrm{f} = \mu_\mathrm{g} \phi_\mathrm{f}\rho_\mathrm{f}; \\
 \mu_\mathrm{g} = k_\mathrm{\mu}
 \frac{C_\mathrm{w}^\mathrm{s}}{K_\mathrm{s} + C_\mathrm{w}^\mathrm{s}}
 \frac{C_\mathrm{w}^\mathrm{O_2}}{K_\mathrm{O_2} + C_\mathrm{w}^\mathrm{O_2}}.
```

Here, $`k_\mathrm{\mu}`$ is the maximum specific growth rate, according to [@Connolly2014],
$`C_\mathrm{w}^\mathrm{s}`$ and $`C_\mathrm{w}^\mathrm{O_2}`$ are the mass concentrations of substrate and oxygen in the water phase and
K<sub>s</sub> and  K<sub>O<sub>2</sub></sub>  are the half-saturation coefficients for substrate and oxygen respectively.

The decay rates are calculated similarly to the growth rates,
except that the specific decay rates of suspended and attached biomass
take different processes into account, increasing inactivation.
For suspended biomass, non-optimal acidic pH is assumed to increase inactivation.
As calcite precipitates mainly in or close to the biofilm, cells may be covered with calcite precipitates
or disrupted by crystals inactivating the affected cells ([@Dupraz2009a] and [@Whiffin2007]).
Consequently, the precipitation rate is assumed to increase the specific decay rate of attached biomass.


```math
r_\mathrm{b}^\mathrm{b} = k_\mathrm{b}^\mathrm{b} C_\mathrm{w}^\mathrm{b}S_\mathrm{w}\phi; \\
 r_\mathrm{b}^\mathrm{f} =  k_\mathrm{b}^\mathrm{f}\phi_\mathrm{f}\rho_\mathrm{f}; \\
k_\mathrm{b}^\mathrm{b}=b_0\left(1+\frac{K_\mathrm{pH}}{m_\mathrm{H^{+}}^2}\right); \\
 k_\mathrm{b}^\mathrm{f} = b_0 + \frac{r_\mathrm{prec} M^\mathrm{c}}{\rho_\mathrm{c}\left(\phi_0 - \phi_\mathrm{c}\right)}
```

The attachment rate quantifies the biomass transfer from
the suspended to the attached state.
As attachment is modeled assuming a first-order kinetic rather than a
sorption-type behavior, with preferential attachment to previously attached biomass:

```math
r_\mathrm{a}= k_\mathrm{a} C^\mathrm{b}_\mathrm{w} \phi S_\mathrm{w};\\
k_\mathrm{a}=c_\mathrm{a,1} \phi_\mathrm{f} + c_\mathrm{a,2},
```

Detachment of biomass from biofilm is assumed to be proportional to the shear stress.
Additionally, the growth contributes to the detachment rate, as vigorously growing
biofilm is typically weaker and as such more susceptible to detachment.
As the model of [@Ebigbo2012] is defined on the Darcy scale but the shear stress is a micro-scale property,
it is approximated using the absolute value of the water-phase potential gradient:

```math
r_\mathrm{d}=k_\mathrm{d} \phi_\mathrm{f} \rho_\mathrm{f};\\
k_\mathrm{d}=c_\mathrm{d}
\left( \phi S_\mathrm{w} \left|\nabla p_\mathrm{w} -  \rho_\mathrm{w} \mathbf{g} \right| \right)^{0.58}
+ \frac{\phi_\mathrm{f}}{\phi_0-\phi_\mathrm{c}} \mu_\mathrm{g}.
```

Here, $`c_\mathrm{d}`$ is a coefficient for the shear-stress-dependent detachment,
$`\left|\nabla p_\mathrm{w} -  \rho_\mathrm{w} \mathbf{g} \right|`$ the absolute value of
the water-phase potential gradient, and $`\phi_0`$ the initial porosity.


## Substrate and oxygen
The consumption of substrate and oxygen is linked to the growth of both suspended and attached biomass
by the yield coefficient Y of substrate. In the case of oxygen, the coefficient F, which is the ration of oxygen consumed per substrate consumed,
is used to express the biomass yield per oxygen consumed:

```math
q^\mathrm{s} = - \frac{r^\mathrm{b}_\mathrm{g} + r^\mathrm{f}_\mathrm{g}}{M^\mathrm{s}Y}; \\
 q^\mathrm{O_2} = F q^\mathrm{s} = - F \frac{ r^\mathrm{b}_\mathrm{g} + r^\mathrm{f}_\mathrm{g}} {M^\mathrm{O_2}Y}.
```

## Supplementary equations
\subsubsection{Permeability and porosity}\label{sec:perm_poro_pw}
The permeability decreases due to biofilm growth and calcite precipitation.
In the model, the reduction of permeability is calculated based on the
reduction of porosity using a Kozenzy-Carman-type relationship:

```math
\frac{K}{K_\mathrm{0}}=\left(\frac{\phi}{\phi_\mathrm{0}}\right)^3 \left(\frac{1-\phi_0}{1-\phi}\right)^2.
```

Here, $`K_\mathrm{0}`$ is the initial permeability,
and $`\phi_0`$ is the initial porosity. The porosity decreases as
the volume fractions of biofilm and calcite increase:

```math
\phi=\phi_\mathrm{0}-\phi_\mathrm{c}-\phi_\mathrm{f}.
```


## Fluid properties
The density of the CO<sub>2</sub> phase is calculated using
an ideal gas assumption as the pressures are low enough.
The viscosity of the CO<sub>2</sub> phase is calculated using
the relation given by [@Fenghour1998].
In these calculations, the effects of the small amounts of water and oxygen are neglected.
The density and the viscosity of the aqueous phase are calculated
according to [@Batzle1992] as a function of salinity.
Sodium, chloride and calcium are considered to contribute to the salinity.


## Phase partitioning of components
The dissolution of CO<sub>2</sub> in the aqueous phase is calculated according
to [@Duan2003] as a function of temperature, pressure and salinity.
In the two-phase case, the concentration of carbonic acid is assumed to be
equal to the solubility while, in the case of the exclusive presence of the aqueous phase,
the concentration of inorganic carbon
is exclusively dependent on the precipitation, dissolution and ureolysis reactions.
In this case, the solubility represents the maximum possible concentration.
For the solubility of oxygen in the aqueous phase, Henry's law is used.





[@Batzle1992]: https://library.seg.org/doi/10.1190/1.1443207  "Seismic Properties of Pore Fluids"
[@Class2002]: https://www.sciencedirect.com/science/article/pii/S0309170802000155 "Numerical simulation of non-isothermal multiphase multicomponent processes in porous media. 2. Applications for the injection of steam and air"
[@Connolly2014]: http://www.ncbi.nlm.nih.gov/pubmed/23835134 "Construction of two ureolytic model organisms for the study of microbially induced calcium carbonate precipitation"
[@Duan2003]: https://www.sciencedirect.com/science/article/pii/S0009254102002632 "An improved model calculating CO<sub>2</sub> solubility in pure water and aqueous NaCl solutions from 273 to 533 K and from 0 to 2000 bar"
[@Chou1989]: https://www.sciencedirect.com/science/article/pii/0009254189900636 "Comparative study of the kinetics and mechanisms of dissolution of carbonate minerals"
[@Compton1989]: https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2427.1989.tb01101.x "The dissolution of calcite in acid waters: Mass transport versus surface control"
[@Dupraz2009a]: http://linkinghub.elsevier.com/retrieve/pii/S0009254109002265 "Experimental and numerical modeling of bacterially induced pH increase and calcite precipitation in saline aquifers"
[@Whiffin2007]: https://www.tandfonline.com/doi/full/10.1080/01490450701436505 "Microbial Carbonate Precipitation as a Soil Improvement Technique"
[@Ebigbo2012]: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011WR011714 "Darcy-scale modeling of microbially induced carbonate mineral precipitation in sand columns"
[@Fenghour1998]: https://aip.scitation.org/doi/10.1063/1.556013 "The Viscosity of Carbon Dioxide"
[@Hommel2015]: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014WR016503 "A revised model for microbially induced calcite precipitation: Improvements	and new insights based on recent experiments"
[@Hommel2016]: https://elib.uni-stuttgart.de/handle/11682/8787 "Modelling biogeochemical and mass transport processes in the subsurface: investigation of microbially induced calcite precipitation"
[@Lauchnor2015]: http://dx.doi.org/10.1111/jam.12804 "Whole cell kinetics of ureolysis by Sporosarcina pasteurii"
[@Mitchell2008]: https://www.sciencedirect.com/science/article/pii/S0896844608002040 "Resilience of planktonic and biofilm cultures to supercritical CO<sub>2</sub>"
