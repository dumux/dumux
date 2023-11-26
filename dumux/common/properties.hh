// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Properties
 *
 * \brief _Defines_ all properties used in Dumux.
 * \note Include this to forward declare properties in your headers.
 * \note This is guaranteed to also include the property system
 */

#ifndef DUMUX_PROPERTIES_HH
#define DUMUX_PROPERTIES_HH

// explicitly guard the include so that the property system
// header doesn't need to be opened and checked all the time
// this include is guaranteed to be here for users of this header
#ifndef DUMUX_PROPERTY_SYSTEM_HH
#include <dumux/common/properties/propertysystem.hh>
#endif // DUMUX_PROPERTY_SYSTEM_HH

namespace Dumux::Properties {

///////////////////////////////////////
// Basic properties of numeric models:
///////////////////////////////////////
DUMUX_DEFINE_PROPERTY(Scalar)                 //!< Property to specify the type of scalar values.
DUMUX_DEFINE_PROPERTY(ModelDefaultParameters) //!< Property which defines the group that is queried for parameters by default
DUMUX_DEFINE_PROPERTY(Grid)                   //!< The DUNE grid type
DUMUX_DEFINE_PROPERTY(PrimaryVariables)       //!< A vector of primary variables
DUMUX_DEFINE_PROPERTY(ModelTraits)            //!< Traits class encapsulating model specifications
DUMUX_DEFINE_PROPERTY(BaseModelTraits)        //!< Model traits to be used as a base for nonisothermal, mineralization ... models
DUMUX_DEFINE_PROPERTY(Problem)                //!< Property to specify the type of a problem which has to be solved
DUMUX_DEFINE_PROPERTY(PointSource)            //!< Property defining the type of point source used
DUMUX_DEFINE_PROPERTY(PointSourceHelper)      //!< Property defining the class that computes which sub control volume point sources belong to
DUMUX_DEFINE_PROPERTY(IOFields)               //!< A class helping models to define input and output fields
DUMUX_DEFINE_PROPERTY(BaseLocalResidual)      //!< The type of the base class of the local residual (specific to a discretization scheme)
DUMUX_DEFINE_PROPERTY(JacobianMatrix)         //!< Type of the global jacobian matrix
DUMUX_DEFINE_PROPERTY(SolutionVector)         //!< Vector containing all primary variable vector of the grid

//! The type of the local residual function, i.e. the equation to be solved. Must inherit
//! from the BaseLocalResidual property and fulfill its interfaces.
DUMUX_DEFINE_PROPERTY(LocalResidual)

//! TODO: Remove this property as soon as the decoupled models are integrated
DUMUX_DEFINE_PROPERTY(LinearSolver)

////////////////////////////////////////////////
// Basic properties regarding balance equations
/////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(UseMoles)                      //!< Property whether to use moles or kg as amount unit for balance equations
DUMUX_DEFINE_PROPERTY(ReplaceCompEqIdx)              //!< The component balance index that should be replaced by the total mass/mole balance
DUMUX_DEFINE_PROPERTY(BalanceEqOpts)                 //!< A class that collects options for the evaluation of the balance equations
DUMUX_DEFINE_PROPERTY(EnableCompositionalDispersion) //!< Property whether to include compositional dispersion
DUMUX_DEFINE_PROPERTY(EnableThermalDispersion)       //!< Property whether to include thermal dispersion

/////////////////////////////////////////////
// Properties used by finite volume schemes:
/////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(ElementBoundaryTypes)           //!< Stores the boundary types on an element
DUMUX_DEFINE_PROPERTY(GridGeometry)                   //!< Grid wrapper creating finite volume geometry and connectivity
DUMUX_DEFINE_PROPERTY(GridVariables)                  //!< The grid variables object managing variable data on the grid (volvars/fluxvars cache)
DUMUX_DEFINE_PROPERTY(VolumeVariables)                //!< The secondary variables within a sub-control volume
DUMUX_DEFINE_PROPERTY(GridVolumeVariables)            //!< The type for a global container for the volume variables
DUMUX_DEFINE_PROPERTY(FluxVariables)                  //!< Container storing the different types of flux variables
DUMUX_DEFINE_PROPERTY(FluxVariablesCache)             //!< Stores data associated with flux vars
DUMUX_DEFINE_PROPERTY(FluxVariablesCacheFiller)       //!< The engine behind the global flux cache (how to fill caches for the stencil)
DUMUX_DEFINE_PROPERTY(GridFluxVariablesCache)         //!< The global vector of flux variable containers

DUMUX_DEFINE_PROPERTY(EnableGridGeometryCache)        //!< Whether to store finite volume geometry for the entire grid
DUMUX_DEFINE_PROPERTY(EnableGridVolumeVariablesCache) //!< If disabled, the volume variables are not stored (reduces memory, but is slower)
DUMUX_DEFINE_PROPERTY(EnableGridFluxVariablesCache)   //!< Specifies if data on flux vars should be stores (faster, but more memory consuming)

/////////////////////////////////////////////////////////////////
// Additional properties used by the cell-centered mpfa schemes:
/////////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(PrimaryInteractionVolume)       //!< The primary interaction volume type
DUMUX_DEFINE_PROPERTY(SecondaryInteractionVolume)     //!< The secondary interaction volume type used e.g. on the boundaries
DUMUX_DEFINE_PROPERTY(DualGridNodalIndexSet)          //!< The type used for the nodal index sets of the dual grid

/////////////////////////////////////////////////////////////
// Properties used by models involving flow in porous media:
/////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(EnergyLocalResidual)                 //!< The local residual of the energy equation
DUMUX_DEFINE_PROPERTY(AdvectionType)                       //!< The type for the calculation the advective fluxes
DUMUX_DEFINE_PROPERTY(SolutionDependentAdvection)          //!< specifies if the parameters for the advective fluxes depend on the solution
DUMUX_DEFINE_PROPERTY(MolecularDiffusionType)              //!< The type for the calculation of the molecular diffusion fluxes
DUMUX_DEFINE_PROPERTY(DispersionFluxType)                  //!< The type for the calculation of the dispersive fluxes
DUMUX_DEFINE_PROPERTY(SolutionDependentMolecularDiffusion) //!< specifies if the parameters for the diffusive fluxes depend on the solution
DUMUX_DEFINE_PROPERTY(HeatConductionType)                  //!< The type for the calculation of the heat conduction fluxes
DUMUX_DEFINE_PROPERTY(CompositionalDispersionModel)        //!< The type for the calculation of the compositional dispersion tensor
DUMUX_DEFINE_PROPERTY(ThermalDispersionModel)              //!< The type for the calculation of the thermal dispersion tensor
DUMUX_DEFINE_PROPERTY(SolutionDependentHeatConduction)     //!< specifies if the parameters for the heat conduction fluxes depend on the solution

DUMUX_DEFINE_PROPERTY(SpatialParams)                       //!< The type of the spatial parameters object
DUMUX_DEFINE_PROPERTY(FluidSystem)                         //!< The type of the fluid system to use
DUMUX_DEFINE_PROPERTY(FluidState)                          //!< The type of the fluid state to use
DUMUX_DEFINE_PROPERTY(SolidSystem)                         //!< The type of the solid system to use
DUMUX_DEFINE_PROPERTY(SolidState)                          //!< The type of the solid state to use
DUMUX_DEFINE_PROPERTY(EffectiveDiffusivityModel)           //!< The employed model for the computation of the effective diffusivity
DUMUX_DEFINE_PROPERTY(ThermalConductivityModel)            //!< Model to be used for the calculation of the effective conductivity
DUMUX_DEFINE_PROPERTY(VelocityOutput)                      //!< specifies the velocity calculation module to be used
DUMUX_DEFINE_PROPERTY(Formulation)                         //!< The formulation of the model
// TODO: is this useful? -> everything is a constraint solver just a different type
DUMUX_DEFINE_PROPERTY(UseConstraintSolver)                 //!< Whether to use a constraint solver for computing the secondary variables

// When using the box method in a multi-phase context, an interface solver might be necessary
DUMUX_DEFINE_PROPERTY(EnableBoxInterfaceSolver)

//////////////////////////////////////////////////////////////
// Additional properties used by the 2pnc and 2pncmin models:
//////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(Chemistry)                      //!< The chemistry class with which equilibrium reactions are solved
DUMUX_DEFINE_PROPERTY(SetMoleFractionsForFirstPhase)  //!< Set the mole fraction in the wetting or nonwetting phase
//////////////////////////////////////////////////////////////
// Additional properties used by the 3pwateroil model:
//////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(OnlyGasPhaseCanDisappear)       //!< reduces the phasestates to threePhases and wnPhaseOnly

/////////////////////////////////////////////////////////////
// Properties used by geomechanical models:
/////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(StressType)                     //!< The type used for the evaluation of stress tensors and forces

/////////////////////////////////////////////////////////////
// Properties used by the staggered-grid discretization method
/////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(NumEqCellCenter)                  //!< The number of equations for cell-centered dofs
DUMUX_DEFINE_PROPERTY(NumEqFace)                        //!< The number of equations for face dofs
DUMUX_DEFINE_PROPERTY(CellCenterSolutionVector)         //!< The solution vector type for cell-centered dofs
DUMUX_DEFINE_PROPERTY(FaceSolutionVector)               //!< The solution vector type for face dofs
DUMUX_DEFINE_PROPERTY(GridFaceVariables)                //!< Global vector containing face-related data
DUMUX_DEFINE_PROPERTY(CellCenterPrimaryVariables)       //!< The primary variables container type for cell-centered dofs
DUMUX_DEFINE_PROPERTY(FacePrimaryVariables)             //!< The primary variables container type for face dofs
DUMUX_DEFINE_PROPERTY(IntersectionMapper)               //!< Specifies the intersection mapper
DUMUX_DEFINE_PROPERTY(StaggeredPrimaryVariables)        //!< The hybrid primary variables container type
DUMUX_DEFINE_PROPERTY(BaseEpsilon)                      //!< A base epsilon for numerical differentiation, can contain multiple values
DUMUX_DEFINE_PROPERTY(FaceVariables)                    //!< Class containing local face-related data
DUMUX_DEFINE_PROPERTY(BoundaryValues)                   //!< Class containing local boundary data
DUMUX_DEFINE_PROPERTY(StaggeredFaceSolution)            //!< A vector containing the solution for a face (similar to ElementSolution)
DUMUX_DEFINE_PROPERTY(EnableGridFaceVariablesCache)     //!< Switch on/off caching of face variables
DUMUX_DEFINE_PROPERTY(UpwindSchemeOrder)                //!< Specifies the order of the upwinding scheme (1 == first order, 2 == second order(tvd methods))

/////////////////////////////////////////////////////////////
// Properties used by the mpnc model
/////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(PressureFormulation)              //! the formulation of the pressure e.g most wetting first

/////////////////////////////////////////////////////////////
// Properties used by the nonequilibrium model
/////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(EquilibriumModelTraits)
DUMUX_DEFINE_PROPERTY(EquilibriumLocalResidual)
DUMUX_DEFINE_PROPERTY(EquilibriumIndices)
DUMUX_DEFINE_PROPERTY(EquilibriumIOFields)
DUMUX_DEFINE_PROPERTY(NumEqBalance)
DUMUX_DEFINE_PROPERTY(EnableThermalNonEquilibrium)
DUMUX_DEFINE_PROPERTY(EnableChemicalNonEquilibrium)
DUMUX_DEFINE_PROPERTY(NumEnergyEqFluid)
DUMUX_DEFINE_PROPERTY(NumEnergyEqSolid)

DUMUX_DEFINE_PROPERTY(NusseltFormulation)
DUMUX_DEFINE_PROPERTY(SherwoodFormulation)

/////////////////////////////////////////////////////////////
// Properties used by free flow models
/////////////////////////////////////////////////////////////

DUMUX_DEFINE_PROPERTY(NormalizePressure) //!<  Returns whether to normalize the pressure term in the momentum balance or not
DUMUX_DEFINE_PROPERTY(ViscousFluxType)   //!< The type for the calculation of the (turbulent) viscous (momentum) fluxes

/////////////////////////////////////////////////////////////
// Properties used by multidomain simulations
/////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(CouplingManager)

/////////////////////////////////////////////////////////////
// Basic properties of by old/deprecated sequential models:
// Do not use this unless you are dealing with such old code
////////////////////////////////////////////////////////////
DUMUX_DEFINE_PROPERTY(TimeManager)

} // end namespace Dumux::Properties

#endif
