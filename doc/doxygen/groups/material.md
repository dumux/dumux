# Doxygen subgroups of the material group

<!-- Binarycoefficients -->

@defgroup Binarycoefficients Relations for binary mixtures
@brief Binary coefficients such as binary diffusion coefficients, Henry coefficients.
@details Binary coefficients describe the relations of a mixture of two components. Typical binary coefficients are Henry coefficients or binary molecular diffusion coefficients.
@ingroup Material

<!-- Chemistry -->

@defgroup Chemistry Chemical constitutive models
@brief Chemical reaction models
@details Chemical reactions can be relevant for all thermodynamic relations for the liquid and gas phase of multiple chemical species. The main purpose is to provide a convenient way to access these relationships via source or sink terms.
@ingroup Material

<!-- Components -->

@defgroup Components Thermodynamical properties of chemical species
@brief Thermodynamical properties of single chemical species or fixed mixtures of species ($ \text{CO}_2, \text{H}_2\text{O}, \text{Air}, ... $)
@details Components provide the thermodynamic relations for the liquid, gaseous and/or solid state of a single
chemical species or a _fixed_ mixture of species. Fluid systems use components to compute thermodynamic quantities of phases. An example would be the dynamic viscosity at different temperatures and pressures.
@ingroup Material

<!-- Components subgroups -->

@defgroup IAPWS IAPWS
@brief Tabulated values according to the International Association for the Properties of Water and Steam (IAPWS)
@ingroup Components

<!-- ConstraintSolvers -->

@defgroup ConstraintSolvers Constraint solvers for thermodynamic constraints of mixtures
@brief Constraint solvers converting primary to secondary variables
@details Constraint solvers are auxiliary tools to make sure that a fluid state is consistent with some thermodynamic constraints. All constraint solvers specify a well defined set of input variables and make sure that the resulting fluid state is consistent with a given set of thermodynamic equations. Constraint solvers connect the thermodynamic relations expressed by fluid systems with the thermodynamic quantities stored by fluid states. Using them is not mandatory for models, but given the fact that some thermodynamic constraints can be quite complex to solve, sharing this code between models makes sense.
@ingroup Material

<!-- EOS -->

@defgroup EOS Equations of State
@brief Thermodynamic equations relating state variables (e.g. temperature, pressure, density)
@details Equations of state (EOS) are provided in the form of auxiliary classes which provide interfaces to query relations between a fluid phase's temperature, pressure, composition and density.
@ingroup Material

<!-- Fluidmatrixinteractions -->

@defgroup Fluidmatrixinteractions Fluid-matrix interactions
@brief Constitutive models for interaction of fluids and solids
@details This module includes constitutive relations such as pc-Sw relations, kr-Sw relations, effective diffusion coefficients, friction laws. The relations depend on both the fluid state as well as material parameters of the solid matrix. For example, in porous media theory, the effective heat conductivity depends on the solid heat conductivity, the fluid heat conductivity, as well as the porosity of the solid and the fluid saturation.
@ingroup Material

<!-- Fluidmatrixinteractions subgroups  -->

@defgroup EffectiveDiffusivity Effective diffusivity in porous media
@brief Laws for calculating effective diffusion coefficients.
@details When averaging over a given volume of a porous medium, diffusion appears effectively restricted since not all volume is accessible to particles and diffusion is hindered by the solid matrix acting as obstacles. Effective diffusivity laws provide constitutive relations for the effective diffusion coefficients based on the solid matrix material parameters and the fluid configuration in the pore space.

The effective diffusion coefficient of component $ \kappa $
in phase $ \alpha $ can be modeled as
\begin{equation}
D^\kappa_{\text{eff},\alpha} = \phi S_\alpha \tau D^\kappa_\alpha,
\end{equation}
where
$ \phi $ is the porosity (volume fraction of the pore space),
$ S_\alpha $ is the saturation of phase $ \alpha $ (the volume
fraction of phase $ \alpha $ being $ n_\alpha = \phi S_\alpha $),
$ D^\kappa_\alpha $ denotes the binary diffusion coefficient of
component $ \kappa $ in phase $ \alpha $, and $ \tau $
is the tortuosity coefficient.

Bear \cite bear1972 reports values of $\tau$ in the range of 0.4 to 0.8.
Note that in some literature the tortuosity $ \lambda $ is used instead
of the tortuosity coefficient $\tau$. The two
quantities are related by $ \lambda = \sqrt{1/\tau} $.

The following laws are implemented:

@ingroup Fluidmatrixinteractions

@defgroup EffectiveHeatConductivity Effective heat conductivity in porous media
@brief Laws for calculating effective heat conductivity coefficients.
@details In porous media, the effective heat conductivity depends on the solid-fluid conductivity ratio, the volume fractions of the constituent phases, and the geometry of the solid-fluid interface \cite aichlmayr2006effective. The following laws are implemented:
@ingroup Fluidmatrixinteractions

@defgroup FrictionLaws Friction Laws
@brief Friction Laws for calculating bottom shear stress
@details Friction laws calculate the stress between the flowing fluid and the bottom,
which is called bottom shear stress. The bottom shear stress is
needed to calculate on the one hand the loss of momentum due to
bottom friction and on the other hand the bedload transport rate.
@ingroup Fluidmatrixinteractions

@defgroup DispersionTensors Dispersion Tensors
@brief Dispersion Tensors for calculating the dispersion in all spatial dimensions
@details Dispersion is caused by particles traveling with different velocity through the porous medium. A particle in the middle of a pore has a higher velociy than a particle close to the grain. Additionaly, due to the tortuosity of the porous medium particles have to take different paths through it. This leads to an effect, which is like an enhanced diffusion, equalizing concentration gradients. 
@ingroup Fluidmatrixinteractions

@defgroup PoreNetwork Pore Network
@brief Constitutive Relations for pore networks models.
@details In the pore-network model, a porous medium is represented as large void spaces called pore body (pore) connected to each other by narrow void spaces called pore throat (throat). The primary variables like pressure and saturation are located at the pore bodies and pore throats determine flow conductivity. Constitutive relations for pore bodies describe capillary pressure-saturation relationships for different pore shapes to be used in two-phase flow simulations. For pore throats, threshold capillary pressures that indicate when a throat is invaded by the non-wetting phase and when it is filled with the wetting phase are provided. Furthermore, relations for computing single-phase and two-phase transmissibility at the pore throats are given.
@ingroup Fluidmatrixinteractions

<!-- FluidStates -->

@defgroup FluidStates Fluid States
@brief Fluid states represent the thermodynamic configuration of a system
@details A fluid state always provides access methods to __all__ thermodynamic quantities, but the concept of a fluid state does not mandate what assumptions are made to store these thermodynamic quantities. Not that fluid states do __not__ make sure that the thermodynamic state which they represent is physically possible, they simply store the variables specified.
@ingroup Material

<!-- FluidSystems -->

@defgroup FluidSystems Fluid Systems
@brief Fluid systems express the thermodynamic relations (functions).
@ingroup Material

<!-- SolidStates -->

@defgroup SolidStates Solid States
@brief Solid states are responsible for representing all relevant
thermodynamic quantities of solid systems.
@details A solid state provides access methods to __all__ thermodynamic quantities, but the concept of a solid state does not mandate what assumptions are made to store these thermodynamic quantities. What solid states also do __not__ do is to make sure that the thermodynamic state which they represent is physically possible.
@ingroup Material

<!-- SolidSystems -->

@defgroup SolidSystems Solid Systems
@brief Solid systems express the thermodynamic relations (functions).
@ingroup Material
