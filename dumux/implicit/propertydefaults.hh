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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \file
 * \brief Default properties for implicit models
 */
#ifndef DUMUX_IMPLICIT_PROPERTY_DEFAULTS_HH
#define DUMUX_IMPLICIT_PROPERTY_DEFAULTS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>
#include <dumux/nonlinear/newtonconvergencewriter.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/timemanager.hh>
#include <dumux/linear/amgbackend.hh>

#include <dumux/porousmediumflow/implicit/fluxvariables.hh>
#include <dumux/porousmediumflow/implicit/fluxvariablescache.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/localresidual.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

#include <dumux/discretization/volumevariables.hh>
#include <dumux/discretization/darcyslaw.hh>
#include <dumux/discretization/fickslaw.hh>
#include <dumux/discretization/fourierslaw.hh>
#include <dumux/discretization/fluxstencil.hh>
#include <dumux/discretization/upwindscheme.hh>

#include <dumux/io/vtkoutputmodulebase.hh>
#include <dumux/porousmediumflow/implicit/velocityoutput.hh>

#include "properties.hh"
#include "model.hh"
#include "assembler.hh"
#include "localjacobian.hh"

#include <dune/common/version.hh>

namespace Dumux {

// forward declarations
template <class TypeTag> class NewtonController;
template <class TypeTag> class AMGBackend;

namespace Properties {
//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for the time manager
SET_TYPE_PROP(ImplicitBase, TimeManager, TimeManager<TypeTag>);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(ImplicitBase,
              GridView,
              typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! use the plain newton method by default
SET_TYPE_PROP(ImplicitBase, NewtonMethod, NewtonMethod<TypeTag>);

//! use the plain newton controller by default
SET_TYPE_PROP(ImplicitBase, NewtonController, NewtonController<TypeTag>);

//! use the plain newton convergence writer by default
SET_TYPE_PROP(ImplicitBase, NewtonConvergenceWriter, NewtonConvergenceWriter<TypeTag>);

//! Mapper for the grid view's vertices.
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
SET_TYPE_PROP(ImplicitBase,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView)>);
#else
SET_TYPE_PROP(ImplicitBase,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);
#endif


//! Mapper for the grid view's elements.
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
SET_TYPE_PROP(ImplicitBase,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView)>);
#else
SET_TYPE_PROP(ImplicitBase,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);
#endif

//! Set the BaseModel to ImplicitModel
SET_TYPE_PROP(ImplicitBase, BaseModel, ImplicitModel<TypeTag>);

//! The volume variable class, to be overloaded by the model
SET_TYPE_PROP(ImplicitBase, VolumeVariables, ImplicitVolumeVariables<TypeTag>);

//! The class that contains the different flux variables (i.e. darcy, diffusion, energy)
//! by default, we set the flux variables to ones for porous media
SET_TYPE_PROP(ImplicitBase, FluxVariables, PorousMediumFluxVariables<TypeTag>);

//! The flux variables cache class, by default the one for porous media
SET_TYPE_PROP(ImplicitBase, FluxVariablesCache, PorousMediumFluxVariablesCache<TypeTag>);

//! The class computing the flux stencil
SET_TYPE_PROP(ImplicitBase, FluxStencil, FluxStencil<TypeTag>);

//! The class applying the upwind scheme
SET_TYPE_PROP(ImplicitBase, UpwindScheme, UpwindScheme<TypeTag>);

//! We use darcy's law as the default for the advective fluxes
SET_TYPE_PROP(ImplicitBase, AdvectionType, DarcysLaw<TypeTag>);

//! We use fick's law as the default for the diffusive fluxes
SET_TYPE_PROP(ImplicitBase, MolecularDiffusionType, FicksLaw<TypeTag>);

//! We use fourier's law as the default for heat conduction fluxes
SET_TYPE_PROP(ImplicitBase, HeatConductionType, FouriersLaw<TypeTag>);

//! The type of a solution for the whole grid at a fixed time
SET_TYPE_PROP(ImplicitBase,
              SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! The type of a solution for a whole element
SET_TYPE_PROP(ImplicitBase,
              ElementSolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! A vector of primary variables
SET_TYPE_PROP(ImplicitBase,
              PrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

//! Boundary types at a single degree of freedom
SET_TYPE_PROP(ImplicitBase,
              BoundaryTypes,
              BoundaryTypes<GET_PROP_VALUE(TypeTag, NumEq)>);

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(ImplicitBase, ImplicitNumericDifferenceMethod, +1);

//! disable jacobian matrix recycling by default
SET_BOOL_PROP(ImplicitBase, ImplicitEnableJacobianRecycling, false);

//! disable partial reassembling by default
SET_BOOL_PROP(ImplicitBase, ImplicitEnablePartialReassemble, false);

//! We do not store the FVGeometry by default
SET_BOOL_PROP(ImplicitBase, EnableGlobalFVGeometryCache, false);

//! We do not store the volume variables by default
SET_BOOL_PROP(ImplicitBase, EnableGlobalVolumeVariablesCache, false);

//! disable flux variables data caching by default
SET_BOOL_PROP(ImplicitBase, EnableGlobalFluxVariablesCache, false);

//! by default, parameters are solution-dependent
SET_BOOL_PROP(ImplicitBase, SolutionDependentAdvection, true);
SET_BOOL_PROP(ImplicitBase, SolutionDependentMolecularDiffusion, true);
SET_BOOL_PROP(ImplicitBase, SolutionDependentHeatConduction, true);

//! specify if we evaluate the permeability in the volume (for discontinuous fields, default)
//! or at the scvf center for analytical permeability fields (e.g. convergence studies)
SET_BOOL_PROP(ImplicitBase, EvaluatePermeabilityAtScvfIP, false);

//! by default, boundary conditions are not constant over time
SET_BOOL_PROP(ImplicitBase, ConstantBoundaryConditions, false);

SET_TYPE_PROP(ImplicitBase, PrimaryVariableSwitch, NoPrimaryVariableSwitch<TypeTag> );

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(ImplicitBase, JacobianMatrix)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

//! use the AMGBackend solver by default
SET_TYPE_PROP(ImplicitBase, LinearSolver, Dumux::AMGBackend<TypeTag> );

// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
SET_SCALAR_PROP(ImplicitBase, LinearSolverResidualReduction, 1e-6);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(ImplicitBase, LinearSolverMaxIterations, 250);

//! set number of equations of the mathematical model as default
SET_INT_PROP(ImplicitBase, LinearSolverBlockSize, GET_PROP_VALUE(TypeTag, NumEq));

//! set number of maximum timestep divisions to 10
SET_INT_PROP(ImplicitBase, ImplicitMaxTimeStepDivisions, 10);

//! Per default we have assume isothermal problems. Set this to true to solve an energy equation
SET_BOOL_PROP(ImplicitBase, EnableEnergyBalance, false);

SET_TYPE_PROP(ImplicitBase, EnergyLocalResidual, EnergyLocalResidual<TypeTag> );

//! Set the upwind weight for the advective term
SET_SCALAR_PROP(ImplicitBase, ImplicitUpwindWeight, 1.0);

//! vtk output
SET_BOOL_PROP(ImplicitBase, VtkAddVelocity, false); //!< Don't reconstruct velocity per default
SET_BOOL_PROP(ImplicitBase, VtkAddProcessRank, true); //!< Add process rank to output per default
SET_TYPE_PROP(ImplicitBase, VtkOutputModule, VtkOutputModuleBase<TypeTag>);
SET_TYPE_PROP(ImplicitBase, VelocityOutput, ImplicitVelocityOutput<TypeTag>);

} // namespace Properties

} // namespace Dumux

#endif
