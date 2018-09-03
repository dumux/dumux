// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Common
 * \file
 *
 * \brief _Declares_ all properties used in Dumux.
 * \note Include this to forward declare properties in your headers.
 */

#ifndef DUMUX_PROPERTIES_HH
#define DUMUX_PROPERTIES_HH

#ifndef DUMUX_PROPERTY_SYSTEM_HH
#include <dumux/common/properties/propertysystem.hh>
#include <dumux/common/properties/propertysystemmacros.hh> // remove this once all macros are gone
#endif

namespace Dumux {
namespace Properties {

///////////////////////////////////////
// Basic properties of numeric models:
///////////////////////////////////////
NEW_PROP_TAG(Scalar);                 //!< Property to specify the type of scalar values.
NEW_PROP_TAG(ModelDefaultParameters); //!< Property which defines the group that is queried for parameters by default
NEW_PROP_TAG(Grid);                   //!< The DUNE grid type
NEW_PROP_TAG(PrimaryVariables);       //!< A vector of primary variables
NEW_PROP_TAG(NumEqVector);            //!< A vector of size number equations that can be used for Neumann fluxes, sources, residuals, ...
NEW_PROP_TAG(GridView);               //!< The type of the grid view according to the grid type
NEW_PROP_TAG(ModelTraits);            //!< Traits class encapsulating model specifications
NEW_PROP_TAG(Problem);                //!< Property to specify the type of a problem which has to be solved
NEW_PROP_TAG(PointSource);            //!< Property defining the type of point source used
NEW_PROP_TAG(PointSourceHelper);      //!< Property defining the class that computes which sub control volume point sources belong to
// TODO: Remove deprecated property VtkOutputFields
NEW_PROP_TAG(VtkOutputFields);        //!< A class helping models to define default vtk output parameters
NEW_PROP_TAG(IOFields);               //!< A class helping models to define input and output fields
NEW_PROP_TAG(BaseLocalResidual);      //!< The type of the base class of the local residual (specific to a discretization scheme)
NEW_PROP_TAG(JacobianMatrix);         //!< Type of the global jacobian matrix
NEW_PROP_TAG(SolutionVector);         //!< Vector containing all primary variable vector of the grid
NEW_PROP_TAG(BoundaryTypes);          //!< Stores the boundary types of a single degree of freedom

//! The type of the local residual function, i.e. the equation to be solved. Must inherit
//! from the BaseLocalResidual property and fulfill its interfaces.
NEW_PROP_TAG(LocalResidual);

//! TODO: Remove this property as soon as the decoupled models are integrated
NEW_PROP_TAG(LinearSolver);

////////////////////////////////////////////////
// Basic properties regarding balance equations
/////////////////////////////////////////////////
// TODO: Integrate UseMoles into BalanceEqOpts
NEW_PROP_TAG(UseMoles);               //!< Property whether to use moles or kg as amount unit for balance equations
NEW_PROP_TAG(ReplaceCompEqIdx);       //!< The component balance index that should be replaced by the total mass/mole balance
NEW_PROP_TAG(BalanceEqOpts);          //!< A class that collects options for the evaluation of the balance equations

/////////////////////////////////////////////
// Properties used by finite volume schemes:
/////////////////////////////////////////////
NEW_PROP_TAG(ElementBoundaryTypes);                //!< Stores the boundary types on an element

NEW_PROP_TAG(FVGridGeometry);                      //!< The type of the global finite volume geometry
NEW_PROP_TAG(EnableFVGridGeometryCache);           //!< specifies if geometric data is saved (faster, but more memory consuming)

NEW_PROP_TAG(VolumeVariables);                     //!< The secondary variables within a sub-control volume
NEW_PROP_TAG(GridVolumeVariables);                 //!< The type for a global container for the volume variables
NEW_PROP_TAG(EnableGridVolumeVariablesCache);      //!< If disabled, the volume variables are not stored (reduces memory, but is slower)
NEW_PROP_TAG(FluxVariables);                       //!< Container storing the different types of flux variables
NEW_PROP_TAG(FluxVariablesCache);                  //!< Stores data associated with flux vars
NEW_PROP_TAG(GridFluxVariablesCache);              //!< The global vector of flux variable containers
NEW_PROP_TAG(EnableGridFluxVariablesCache);        //!< specifies if data on flux vars should be saved (faster, but more memory consuming)
NEW_PROP_TAG(GridVariables);                       //!< The grid variables object managing variable data on the grid (volvars/fluxvars cache)

/////////////////////////////////////////////////////////////////
// Additional properties used by the cell-centered mpfa schemes:
/////////////////////////////////////////////////////////////////
NEW_PROP_TAG(PrimaryInteractionVolume);            //!< The primary interaction volume type
NEW_PROP_TAG(SecondaryInteractionVolume);          //!< The secondary interaction volume type used e.g. on the boundaries
NEW_PROP_TAG(DualGridNodalIndexSet);               //!< The type used for the nodal index sets of the dual grid

/////////////////////////////////////////////////////////////
// Properties used by models involving flow in porous media:
/////////////////////////////////////////////////////////////
NEW_PROP_TAG(EnergyLocalResidual);                 //!< The local residual of the energy equation
NEW_PROP_TAG(AdvectionType);                       //!< The type for the calculation the advective fluxes
NEW_PROP_TAG(SolutionDependentAdvection);          //!< specifies if the parameters for the advective fluxes depend on the solution
NEW_PROP_TAG(MolecularDiffusionType);              //!< The type for the calculation of the molecular diffusion fluxes
NEW_PROP_TAG(SolutionDependentMolecularDiffusion); //!< specifies if the parameters for the diffusive fluxes depend on the solution
NEW_PROP_TAG(HeatConductionType);                  //!< The type for the calculation of the heat conduction fluxes
NEW_PROP_TAG(SolutionDependentHeatConduction);     //!< specifies if the parameters for the heat conduction fluxes depend on the solution

NEW_PROP_TAG(SpatialParams);                       //!< The type of the spatial parameters object
NEW_PROP_TAG(FluidSystem);                         //!< The type of the fluid system to use
NEW_PROP_TAG(FluidState);                          //!< The type of the fluid state to use
NEW_PROP_TAG(SolidSystem);                         //!< The type of the solid system to use
NEW_PROP_TAG(SolidState);                           //!< The type of the solid state to use
NEW_PROP_TAG(PrimaryVariableSwitch);               //!< The primary variable switch needed for compositional models
NEW_PROP_TAG(EffectiveDiffusivityModel);           //!< The employed model for the computation of the effective diffusivity
NEW_PROP_TAG(ThermalConductivityModel);            //!< Model to be used for the calculation of the effective conductivity
NEW_PROP_TAG(VelocityOutput);                      //!< specifies the velocity calculation module to be used
NEW_PROP_TAG(Formulation);                         //!< The formulation of the model
// TODO: is this useful? -> everything is a constraint solver just a different type
NEW_PROP_TAG(UseConstraintSolver);                 //!< Whether to use a contraint solver for computing the secondary variables

// When using the box method in a multi-phase context, an interface solver might be necessary
NEW_PROP_TAG(EnableBoxInterfaceSolver);

//////////////////////////////////////////////////////////////
// Additional properties used by the 2pnc and 2pncmin models:
//////////////////////////////////////////////////////////////
NEW_PROP_TAG(Chemistry);                           //!< The chemistry class with which solves equlibrium reactions
NEW_PROP_TAG(SetMoleFractionsForFirstPhase);       //!< Set the mole fraction in the wetting or non-wetting phase

//////////////////////////////////////////////////////////////
// Additional properties used by the richards model
//////////////////////////////////////////////////////////////
NEW_PROP_TAG(EnableWaterDiffusionInAir); //!< Property for turning Richards into extended Richards

//////////////////////////////////////////////////////////////
// Additional properties used by the 3pwateroil model:
//////////////////////////////////////////////////////////////
NEW_PROP_TAG(OnlyGasPhaseCanDisappear); //!< reduces the phasestates to threePhases and wnPhaseOnly

/////////////////////////////////////////////////////////////
// Properties used by geomechanical models:
/////////////////////////////////////////////////////////////
NEW_PROP_TAG(StressType);       //!< The type used for the evaluation of stress tensors and forces

/////////////////////////////////////////////////////////////
// Properties used by the staggered-grid discretization method
/////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumEqCellCenter);                     //!< The number of equations for cell-centered dofs
NEW_PROP_TAG(NumEqFace);                           //!< The number of equations for face dofs
NEW_PROP_TAG(CellCenterSolutionVector);            //!< The solution vector type for cell-centered dofs
NEW_PROP_TAG(FaceSolutionVector);                  //!< The solution vector type for face dofs
NEW_PROP_TAG(GridFaceVariables);                   //!< Global vector containing face-related data
NEW_PROP_TAG(CellCenterPrimaryVariables);          //!< The primary variables container type for cell-centered dofs
NEW_PROP_TAG(FacePrimaryVariables);                //!< The primary variables container type for face dofs
NEW_PROP_TAG(IntersectionMapper);                  //!< Specifies the intersection mapper
NEW_PROP_TAG(StaggeredPrimaryVariables);           //!< The hybrid primary variables container type
NEW_PROP_TAG(BaseEpsilon);                         //!< A base epsilon for numerical differentiation, can contain multiple values
NEW_PROP_TAG(FaceVariables);                       //!< Class containing local face-related data
NEW_PROP_TAG(BoundaryValues);                      //!< Class containing local boundary data
NEW_PROP_TAG(StaggeredFaceSolution);               //!< A vector containing the solution for a face (similar to ElementSolution)
NEW_PROP_TAG(EnableGridFaceVariablesCache);      //!< Switch on/off caching of face variables

/////////////////////////////////////////////////////////////
// Properties used by the mpnc model
/////////////////////////////////////////////////////////////

NEW_PROP_TAG(PressureFormulation); //! the formulation of the pressure e.g most wetting first

/////////////////////////////////////////////////////////////
// Properties used by the nonequilibrium model
/////////////////////////////////////////////////////////////
NEW_PROP_TAG(EquilibriumModelTraits);
NEW_PROP_TAG(EquilibriumLocalResidual);
NEW_PROP_TAG(EquilibriumIndices);
NEW_PROP_TAG(EquilibriumIOFields);
NEW_PROP_TAG(NumEqBalance);
NEW_PROP_TAG(EnableThermalNonEquilibrium);
NEW_PROP_TAG(EnableChemicalNonEquilibrium);
NEW_PROP_TAG(NumEnergyEqFluid);
NEW_PROP_TAG(NumEnergyEqSolid);

NEW_PROP_TAG(AwnSurface);
NEW_PROP_TAG(AwsSurface);
NEW_PROP_TAG(AnsSurface);
NEW_PROP_TAG(NusseltFormulation);
NEW_PROP_TAG(SherwoodFormulation);

/////////////////////////////////////////////////////////////
// Properties used by free flow models
/////////////////////////////////////////////////////////////

NEW_PROP_TAG(NormalizePressure); //!<  Returns whether to normalize the pressure term in the momentum balance or not

/////////////////////////////////////////////////////////////
// Properties used by multidomain simulations
/////////////////////////////////////////////////////////////
NEW_PROP_TAG(CouplingManager);

///////////////////////////////////////
// Basic properties of sequential models:
///////////////////////////////////////
NEW_PROP_TAG(TimeManager);

} // end namespace Properties
} // end namespace Dumux

#endif // DUMUX_PROPERTIES_HH
