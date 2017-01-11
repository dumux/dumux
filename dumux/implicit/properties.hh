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
#ifndef DUMUX_IMPLICIT_PROPERTIES_HH
#define DUMUX_IMPLICIT_PROPERTIES_HH

#include <dumux/common/propertysystem.hh>

#include <dumux/common/basicproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/implicit/adaptive/gridadaptproperties.hh>

/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the ImplicitModel.
 */
namespace Dumux
{

namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for fully-implicit models
NEW_TYPE_TAG(ImplicitBase, INHERITS_FROM(NewtonMethod, LinearSolverTypeTag, NumericModel, GridAdapt));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(BaseModel); //!< The type of the base class of the model
NEW_PROP_TAG(NumEq); //!< Number of equations in the system of PDEs
NEW_PROP_TAG(BaseLocalResidual); //!< The type of the base class of the local residual
NEW_PROP_TAG(LocalResidual); //!< The type of the local residual function
NEW_PROP_TAG(EnergyLocalResidual); //! The local residual of the energy equation
NEW_PROP_TAG(LocalJacobian); //!< The type of the local jacobian operator

NEW_PROP_TAG(SubControlVolume);//!< The type of the sub control volume
NEW_PROP_TAG(SubControlVolumeFace); //!< The type of the sub control volume face
NEW_PROP_TAG(FVElementGeometry); //!< The type of the local finite volume geometry (iterators over scvs, scvfs)
NEW_PROP_TAG(GlobalFVGeometry); //!< The type of the global finite volume geometry
NEW_PROP_TAG(EnableGlobalFVGeometryCache); //! specifies if geometric data is be saved (faster, but more memory consuming)

NEW_PROP_TAG(JacobianAssembler); //!< Assembles the global jacobian matrix
NEW_PROP_TAG(JacobianMatrix); //!< Type of the global jacobian matrix
NEW_PROP_TAG(BoundaryTypes); //!< Stores the boundary types of a single degree of freedom
NEW_PROP_TAG(ElementBoundaryTypes); //!< Stores the boundary types on an element

NEW_PROP_TAG(PrimaryVariables); //!< A vector of primary variables within a sub-control volume
NEW_PROP_TAG(SolutionVector); //!< Vector containing all primary variable vector of the grid
NEW_PROP_TAG(ElementSolutionVector); //!< A vector of primary variables within a sub-control volume

NEW_PROP_TAG(EnableGlobalVolumeVariablesCache); //!< If disabled, the volume variables are not stored (reduces memory, but is slower)
NEW_PROP_TAG(VolumeVariables);  //!< The secondary variables within a sub-control volume
NEW_PROP_TAG(ElementVolumeVariables);  //!< The type for a local (element/stencil) container for the volume variables
NEW_PROP_TAG(GlobalVolumeVariables);  //!< The type for a global container for the volume variables
NEW_PROP_TAG(FluxVariables); //!< Container storing the different types of flux variables
NEW_PROP_TAG(FluxVariablesCache); //!< Stores data associated with flux vars (if enabled)
NEW_PROP_TAG(GlobalFluxVariablesCache); //!< The global vector of flux variable containers
NEW_PROP_TAG(ElementFluxVariablesCache); //!< A local vector of flux variable caches per element
NEW_PROP_TAG(EnableGlobalFluxVariablesCache); //! specifies if data on flux vars should be saved (faster, but more memory consuming)
NEW_PROP_TAG(BoundaryVariables); //!< Data required to calculate fluxes over boundary faces in cc models(outflow)
NEW_PROP_TAG(ConstantBoundaryConditions); //!< specifies whether or not the boundary conditions are constant over time
NEW_PROP_TAG(DiscretizationMethod); //!< Property for the used discretization method
NEW_PROP_TAG(PrimaryVariableSwitch); //!< The primary variable switch needed for compositional models

// Specify the forms of fluxes that should be considered in the model
// also, specify their corresponding flux variables
NEW_PROP_TAG(EnableAdvection); //! specifies if advection is considered in the model
NEW_PROP_TAG(AdvectionType); //! The type for the calculation the advective fluxes
NEW_PROP_TAG(SolutionDependentAdvection); //!< specifies if the parameters for the advective fluxes depend on the solution
NEW_PROP_TAG(EnableMolecularDiffusion); //! specifies if molecular diffusive fluxes are considered in the model
NEW_PROP_TAG(MolecularDiffusionType); //! The type for the calculation of the molecular diffusion fluxes
NEW_PROP_TAG(SolutionDependentMolecularDiffusion); //!< specifies if the parameters for the diffusive fluxes depend on the solution
NEW_PROP_TAG(EnableEnergyBalance); //! Specifies if the model solves an energy equation
NEW_PROP_TAG(HeatConductionType); //! The type for the calculation of the heat conduction fluxes
NEW_PROP_TAG(SolutionDependentHeatConduction); //!< specifies if the parameters for the heat conduction fluxes depend on the solution

// vtk output
NEW_PROP_TAG(VtkAddVelocity); //!< specifies if an element velocity it reconstructed for the output
NEW_PROP_TAG(VtkAddProcessRank); //!< specifies if the process rank should be added the output

// stencils
NEW_PROP_TAG(StencilsVector); //! The type of the global vector of stencils per element

// high level simulation control
NEW_PROP_TAG(TimeManager);  //!< Manages the simulation time
NEW_PROP_TAG(NewtonMethod);     //!< The type of the newton method
NEW_PROP_TAG(NewtonController); //!< The type of the newton controller

//! Specify whether the jacobian matrix of the last iteration of a
//! time step should be re-used as the jacobian of the first iteration
//! of the next time step.
NEW_PROP_TAG(ImplicitEnableJacobianRecycling);

//! Specify whether the jacobian matrix should be only reassembled for
//! elements where at least one vertex is above the specified
//! tolerance
NEW_PROP_TAG(ImplicitEnablePartialReassemble);

/*!
 * \brief Specify which kind of method should be used to numerically
 * calculate the partial derivatives of the residual.
 *
 * -1 means backward differences, 0 means central differences, 1 means
 * forward differences. By default we use central differences.
 */
NEW_PROP_TAG(ImplicitNumericDifferenceMethod);

// mappers from local to global indices

//! mapper for vertices
NEW_PROP_TAG(VertexMapper);
//! mapper for elements
NEW_PROP_TAG(ElementMapper);
//! mapper for degrees of freedom
NEW_PROP_TAG(DofMapper);

//! the maximum allowed number of timestep divisions for the
//! Newton solver
NEW_PROP_TAG(ImplicitMaxTimeStepDivisions);

//! indicate whether discretization is box or not
NEW_PROP_TAG(ImplicitIsBox);

//! the upwind weight for the advective term
NEW_PROP_TAG(ImplicitUpwindWeight);

//! the local fem space used for the AMG backend
NEW_PROP_TAG(ImplicitLocalFemMap);

//! the PDELab backend used for the AMG backend
NEW_PROP_TAG(ImplicitPDELabBackend);

}
}

// \}

#endif
