# Flow and transport in porous media

[TOC]

This section gives a brief introduction to the general physics concepts
used for single and multi-phase flow in porous media in many DuMux models.
The following description mostly concerns models on the Darcy (homogenized) continuum scale
and assumes the existence of a representative elementary volume (REV) for
which average quantities are defined. For a more comprehensive treatment of
the mathematical modeling and description of the physical processes, we recommend
the references @cite CLASS2007, @cite helmig1997multiphase. (Before we start, we remark that DuMux can be used to solve general
conservation equations such as the Navier-Stokes equations and also support other modeling concepts
such as pore-networks models to simulate pore-scale processes. Here, we focus on models in the scope of Darcy's law.)

![](multiphase_processes.svg){html: width=70%}

We start with a few basic definitions and introduce some notation.

__Phase.__ A _phase_ is defined as a continuum having distinct properties (e.g. density and viscosity).
If phases are miscible, they contain dissolved portions of the substance of the other phase.
Fluid and solid phases are distinguished. The fluid phases have different affinities to the solid phases. The phase, which has a higher affinity to the solid phases is referred to as the (more) wetting phase. In the case of two phases, the less wetting one is called the nonwetting phase.

For compositional multi-phase models, fluid phases may be composed of several components, while the solid phases are assumed to consist exclusively of a single component.

__Component.__ The term _component_ stands for constituents of the phases which
can be associated with a unique chemical species or, more generally, with
a group of species exploiting similar physical behavior. For example,
the figure at the beginning of the section shows
shows a water-gas-NAPL system composed of the phases water (subscript
$\text{w}$), gas ($\text{g}$), and NAPL ($\text{n}$). These phases are
composed of the components water (superscript $\text{w}$), the pseudo-component
air ($\text{a}$), and an organic contaminant ($\text{c}$).

The composition of the components in a phase can influence the phase properties. Furthermore, for mass transfer, the phase behavior is quite different from the component behavior.

__Thermodynamic equilibrium.__ For the non-isothermal, multi-phase, multi-component processes in porous media
we state the assumption of __local thermodynamic equilibrium__.
_Chemical equilibrium_ means that the mass/mole fractions of a component in
different phases are in equilibrium.
_Thermal equilibrium_ assumes the same temperature for all considered phases.
_Mechanical equilibrium_ means that the forces at the fluid-fluid and solid-fluid phase
interfaces are in balance and the interfaces are not moving. Assuming that one of these
conditions hold locally (within an REV) can be used to simplified the governing equations.
For instance, it can often be assumed due to slow transport processes in porous media that
all phases have the same temperature. In other words, reaching thermal equilibrium in an REV
is a process much faster than time scales of interest.

__Notation (convention).__ The subscript index $\alpha$, e.g. $\text{w}$, $\text{n}$, and $\text{g}$ refers
to the phase, while the superscript $\kappa$, e.g. $\text{w}$, $\text{a}$, and $\text{c}$
refers to the component. The following __symbols__ are often used in the DuMux documentation of
the mathemetical models for flow in porous media.

| Symbol    | Description | Symbol | Description |
| -------- | ------- | ------- | ------- |
$p_\alpha$ | phase pressure | $\phi$ | porosity |
$T$ | temperature | $K$ | intrinsic permeability tensor |
$S_\alpha$ | phase saturation | $\tau$ | tortuosity |
$x_\alpha^\kappa$ | mole fraction of component $\kappa$ in phase $\alpha$ | $\boldsymbol{g}$ | gravitational acceleration |
$X_\alpha^\kappa$ | mass fraction of component $\kappa$ in phase $\alpha$ | $q^\kappa_\alpha$ | volume source term of $\kappa$ in $\alpha$ |
$\varrho_{\text{mol},\alpha}$ | molar density of phase $\alpha$ | $u_\alpha$ | specific internal energy |
$\varrho_{\alpha}$ | mass density of phase $\alpha$ | $h_\alpha$ | specific enthalpy |
$M$ | molar mass of a phase or component | $c_\text{s}$ | specific heat enthalpy |
$k_{\text{r}\alpha}$ | relative permeability | $\lambda_\text{pm}$ | heat conductivity |
$\mu_\alpha$ | phase viscosity | $q^h$ | heat source term |
$D_\alpha^\kappa$ | diffusivity of component $\kappa$ in phase $\alpha$ | $\boldsymbol{v}_{\alpha}$  | velocity |

## Porous medium properties

The __porosity__ $\phi$ is defined as the fraction of the volume occupied by fluids in an REV $V_\mathrm{fluid}$
divided by the total volume of the REV $V_\mathrm{total}$.

\begin{equation}
\phi=\frac{V_\mathrm{fluid}}{V_\mathrm{total}}=1-\frac{V_\mathrm{solid}}{V_\mathrm{total}}.
\end{equation}

The __intrinsic permeability__ is a measure on the REV scale of the ease of fluid flow through porous media.
It relates the potential gradient and the resulting flow velocity in the Darcy equation.
As the porous medium may have a structure leading to preferential flow in certain directions,
intrinsic permeability is in general a tensorial quantity $\mathbf{K}$.
For isotropic porous media, it can be reduced to a scalar quantity $K$.

## Mixture properties

The __composition__ of a phase is described by __mass__ or __mole fractions__ of the components.
The mole fraction $x^\kappa_\alpha$ of component $\kappa$ in phase $\alpha$ is defined as:

\begin{equation}
x^\kappa_\alpha = \frac{n^\kappa_\alpha}{\sum_i n^i_\alpha},
\end{equation}
where $n^\kappa_\alpha$ is the number of moles of component $\kappa$ in phase $\alpha$.
The mass fraction $X^\kappa_\alpha$ is defined similarly
using the mass of component $\kappa$ in phase $\alpha$ instead of $n^\kappa_\alpha$,
$X^\kappa_\alpha = \frac{\mathrm{mass}^\kappa_\alpha}{\mathrm{mass}^{total}_\alpha}$.

The molar mass $M^\kappa$ of the component $\kappa$ relates the mass fraction
to the mole fraction and vice versa.

## Fluid properties

The most important fluid properties
to describe fluid flow on the REV scale are density and viscosity.

The __density__ $\rho_\alpha$ of a fluid phase $\alpha$ is defined as the ratio of its mass to its volume
$(\rho_\alpha = \frac{\mathrm{mass_\alpha}}{\mathrm{volume_\alpha}})$ while
the molar density $\rho_{\mathrm{mol},\alpha}$ is defined as the ratio of the number of moles per volume
$(\rho_{\mathrm{mol},\alpha} = \frac{\mathrm{moles_\alpha}}{\mathrm{volume_\alpha}})$.

The __dynamic viscosity__ $\mu_\alpha$ characterizes the resistance of a fluid to flow.
As density, it is a fluid phase property.
For Newtonian fluids, the dynamic viscosity relates the shear stress $\tau_\mathrm{s}$ to the
velocity gradient:

\begin{equation}
\tau_\mathrm{s} = \mu_\alpha \frac{d v_{\alpha,\,x}}{d y}.
\end{equation}

The __kinematic viscosity__ $\nu_\alpha := \frac{\mu_\alpha}{\rho_\alpha}$ is more often used
in the description of incompressible systems.

Both density and viscosity generally depend on pressure, temperature, and phase composition.

## Fluid-fluid and Fluid-solid interactions

If more than a single fluid is present in the porous medium,
the fluids interact with each other and the solids,
which leads to additional properties for multi-phase systems.

The __saturation__ $S_\alpha$ of a phase $\alpha$ is defined as the ratio of the volume occupied
by that phase to the total pore volume within an REV.
As all pores are filled with some fluid, the sum of the saturations of all present phases is equal to one.

__Capillary pressure.__ Immiscible fluids form a sharp interface as a result of differences in their intermolecular forces
translating into different adhesive and cohesive forces at the fluid-fluid and fluid-fluid-solid interfaces
creating interfacial tension on the microscale.
From the mechanical equilibrium which has also to be satisfied at the interface,
a difference between the pressures of the fluid phases results defined as
the capillary pressure (phase difference pressure) $p_\mathrm{c}$:

\begin{equation}
p_\mathrm{c} = p_\mathrm{n} - p_\mathrm{w}.
\end{equation}

On the microscale, $p_\mathrm{c}$ can be calculated from the surface tension
according to the Laplace equation.

On the REV scale, however, capillary pressure needs to be defined by quantities of that scale.
Several empirical relations provide expressions to link $p_\mathrm{c}$ to the wetting-phase saturation $S_\mathrm{w}$.
An example is the relation given by Brooks and Corey @cite brooks1964hydrau
to determine $p_\mathrm{c}$ based on
$S_\mathrm{e}$, which is the effective wetting-phase saturation,
the entry pressure $p_\mathrm{d}$, and the parameter $\lambda$ describing the pore-size distribution:

\begin{equation}
p_\mathrm{c} = p_\mathrm{d} S_\mathrm{e}^{-\frac{1}{\lambda}}, \quad S_\mathrm{e} = \frac{S_\mathrm{w}-S_\mathrm{w,r}}{1-S_\mathrm{w,r}},
\end{equation}

where $S_\mathrm{w,r}$ is the residual wetting phase saturation which cannot be displaced
by another fluid phase and remains in the porous medium.

__Relative permeability.__ The presence of two fluid phases in the porous medium reduces the space available for flow
for each of the fluid phases. This increases the resistance to flow of the phases, which is accounted for by the means of
the relative permeability $k_\mathrm{r,\alpha} \leq 0$, which scales the intrinsic permeability.
The relative permeability strongly and nonlinearly depends on the saturation.
The relations describing the relative permeabilities of the wetting and nonwetting phase are different
as the wetting phase predominantly occupies small pores and the edges of larger pores while the
nonwetting phases occupies large pores.
The relative permeabilities for the wetting phase $k_\mathrm{r,w}$ and the nonwetting phase
in a two-fluid-phase system can be modeled, for instance following Brooks and Corey @cite brooks1964hydrau, by

\begin{equation}
k_\mathrm{r,w} = S_\mathrm{e}^{\frac{2+3\lambda}{\lambda}}, \quad k_\mathrm{r,n} = \left( 1- S_\mathrm{e}\right)^2 \left( 1- S_\mathrm{e}^{\frac{2+\lambda}{\lambda}}\right).
\end{equation}

Also see Dumux::FluidMatrix::BrooksCorey for where these constitutive relations are implemented in DuMux.

## Transport processes in porous media

On the macro-scale, mass transport can be characterized by the driving force of the
transport process. Pressure gradients result in the advective transport of a fluid phase
and all the components constituting the phase,
while concentration gradients result in the diffusion of a component within a phase.

__Advective transport__ is determined by the flow field.
On the macro-scale, the Darcy or filter velocity $\mathbf{v}$ is calculated using the Darcy equation
depending on the potential gradient $(\nabla p_\alpha - \rho_\alpha \mathbf{g})$,
accounting for both pressure difference and gravitation,
the intrinsic permeability of the porous medium,
and the viscosity $\mu$ of the fluid phase:

\begin{equation}
\mathbf{v}=-\frac{\mathbf{K}}{\mu}(\nabla p - \rho \mathbf{g}).
\end{equation}

$\mathbf{v}$ is proportional to $(\nabla p - \rho \mathbf{g})$ with the proportional factor $\frac{\mathbf{K}}{\mu}$.
This equation can be extended to multi-phase flow by considering phase velocities $\mathbf{v}_{\alpha}$ of phase $\alpha$
and modeling phase interaction through the relative permeability $k_\mathrm{r,\alpha}$,

\begin{equation}
\mathbf{v}_{\alpha}=-\frac{k_\mathrm{r,\alpha}\mathbf{K}}{\mu_{\alpha}}(\nabla p_{\alpha} - \rho_{\alpha} \mathbf{g})
\end{equation}

__Molecular diffusion__ is a process determined by the concentration gradient.
It is commonly modeled as Fickian diffusion following Fick's first law:

\begin{equation}
\mathbf{j_d}=-\rho_{\alpha} D^\kappa_\alpha \nabla X^\kappa_\alpha,
\end{equation}

where $D^\kappa_\alpha$ is the molecular diffusion coefficient of component $\kappa$ in phase $\alpha$.
In a porous medium, the actual path lines are tortuous due to the impact of the solid matrix.
This tortuosity and the impact of the presence of multiple fluid phases
is accounted for by using an effective diffusion coefficient $D^\kappa_\mathrm{pm, \alpha}$:

\begin{equation}
D^\kappa_\mathrm{pm, \alpha}= \phi \tau_\alpha S_\alpha D^\kappa_\alpha,
\end{equation}

where $\tau_\alpha$ is the tortuosity of phase $\alpha$ (see Dumux::DiffusivityConstantTortuosity).

## Gas mixing laws

Prediction of the $p$-$\varrho$-$T$ behavior of gas mixtures is typically based on one of two concepts: Dalton's law or Amagat's law.
Both laws make the same predictions for ideal gases but differ for non-ideal gas mixtures.
In the following the two concepts will be explained in more detail.

__Dalton's law__ assumes that the gases in the mixture are non-interacting (with each other) and each gas independently applies its own pressure (partial pressure), the sum of which is the total pressure:

\begin{equation}
p = \sum_{i}^{}p_i.
\end{equation}
Here $p_i$ refers to the partial pressure of component i.
As an example, if two equal volumes of gas A and gas B are mixed, the volume of the mixture stays the same but the pressures add up (see figure)

![](dalton1.svg)

The density of the mixture, $\varrho$, can be calculated as follows:
\begin{equation}
\varrho = \frac{m}{V} = \frac{m_\mathrm{A} + m_\mathrm{B}}{V} = \frac{\varrho_\mathrm{A} V + \varrho_\mathrm{B} V}{V} = \varrho_\mathrm{A} + \varrho_\mathrm{B},
\end{equation}

or for an arbitrary number of gases:
\begin{equation}
\varrho = \sum_{i}^{} \varrho_i, \quad \varrho_m = \sum_{i}^{} \varrho_{m,i}.
\end{equation}

__Amagat's law__ assumes that the volumes of the component gases are additive. The interactions of the different gases are the same as the average interactions of the components. This is known as Amagat's law:

\begin{equation}
V = \sum_{i}^{}V_i.
\end{equation}

As an example, if two volumes of gas A and B at equal pressure are mixed, the pressure of the mixture stays the same, but the volumes add up.

![](dalton2.svg)

The density of the mixture, $\varrho$, can be calculated as follows:
\begin{equation}
\varrho = \frac{m}{V} = \frac{m}{V_\mathrm{A} + V_\mathrm{B}} = \frac{m}{\frac{m_\mathrm{A}}{\varrho_\mathrm{A}} + \frac{m_\mathrm{B}}{\varrho_\mathrm{B}}} =
\frac{m}{\frac{X_\mathrm{A} m}{\varrho_\mathrm{A}} + \frac{X_\mathrm{B} m}{\varrho_\mathrm{B}}} = \frac{1}{\frac{X_\mathrm{A}}{\varrho_\mathrm{A}} + \frac{X_\mathrm{B}}{\varrho_\mathrm{B}}},
\end{equation}

or for an arbitrary number of gases:

\begin{equation}
\varrho = \frac{1}{\sum_{i}^{}\frac{X_i}{\varrho_i}}, \quad  \varrho_m = \frac{1}{\sum_{i}^{}\frac{x_i}{\varrho_{m,i}}}.
\end{equation}

__Ideal gases__. An ideal gas is defined as a gas whose molecules are spaced so far apart that the behavior of a molecule is not influenced by the presence of other molecules.
This assumption is usually valid at low pressures and high temperatures. The ideal gas law states that, for one gas:

\begin{equation}
p = \varrho \frac{RT}{M} = \varrho_m RT.
\end{equation}

Using the assumption of ideal gases and Dalton's law (or equivalently Amagat's law)
leads to the following expression for the (mass) density and molar density of the mixture:

\begin{equation}
\varrho = \frac{p}{RT} \sum_{i}^{}M_i x_i, \quad \varrho_m = \frac{p}{RT}.
\end{equation}

Also see Dumux::IdealGas for where this is implemented in DuMux.

## Porous medium flow models

A list of porous medium flow models implemented in DuMux can be found in @ref PorousmediumflowModels.
The module description found under the link features a description of the governing equations
describing each mathematical model.
