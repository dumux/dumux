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
@brief Thermodynamical properties of single chemical species or fixed mixtures of species (\f$ \text{CO}_2, \text{H}_2\text{O}, \text{Air}, ... \f$)
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
@brief Constitutive relations such as pc-Sw relations, kr-Sw relations, effective diffusion coefficients, friction laws
@details Constitutive models for interaction of fluids and solids. The relations depend on the fluid state as well as material parameters of the matrix. For example, in porous media theory, capillary pressure is often expressed as a function of the phase saturation and some shape parameter \f$\lambda\f$ which is dependent on the material (Brooks-Corey model).
@ingroup Material

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
