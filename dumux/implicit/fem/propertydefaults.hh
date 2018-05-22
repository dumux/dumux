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
 * \ingroup BoxModel
 * \file
 * \brief Default properties for fem models
 */
#ifndef DUMUX_FEM_PROPERTY_DEFAULTS_HH
#define DUMUX_FEM_PROPERTY_DEFAULTS_HH

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/timemanager.hh>
#include <dumux/linear/amgbackend.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fem/ipdata.hh>

#include "localresidual.hh"
#include "localjacobian.hh"
#include "assembler.hh"
#include "properties.hh"
#include "model.hh"

//#include <dumux/freeflow/stokes_fem/secondaryvariables.hh>

namespace Dumux
{
// forward declarations
template<class TypeTag> class FemLocalResidual;
template<class TypeTag> class AMGBackend;

namespace Properties
{
//! define the VolumeVariables
//SET_TYPE_PROP(BoxStokes, SecondaryVariables, StokesSecondaryVariables<TypeTag>);


//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for the time manager
SET_TYPE_PROP(FemModel, TimeManager, TimeManager<TypeTag>);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(FemModel,
              GridView,
              typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! use the plain newton method by default
SET_TYPE_PROP(FemModel, NewtonMethod, NewtonMethod<TypeTag>);

//! use the plain newton controller by default
SET_TYPE_PROP(FemModel, NewtonController, NewtonController<TypeTag>);

//! Mapper for the grid view's vertices.
SET_TYPE_PROP(FemModel,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);

//! Mapper for the grid view's elements.
SET_TYPE_PROP(FemModel,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);

//! Set the BaseModel to FemImplicitModel
SET_TYPE_PROP(FemModel, BaseModel, FemImplicitModel<TypeTag>);

//! The type of a solution for the whole grid at a fixed time
SET_TYPE_PROP(FemModel,
              SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! The type of a solution for a whole element
SET_TYPE_PROP(FemModel,
              ElementSolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! A vector of primary variables
SET_TYPE_PROP(FemModel,
              PrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

//! Boundary types at a single degree of freedom
SET_TYPE_PROP(FemModel,
              BoundaryTypes,
              BoundaryTypes<GET_PROP_VALUE(TypeTag, NumEq)>);

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(FemModel, ImplicitNumericDifferenceMethod, +1);

//! disable jacobian matrix recycling by default
SET_BOOL_PROP(FemModel, ImplicitEnableJacobianRecycling, false);

//! disable partial reassembling by default
SET_BOOL_PROP(FemModel, ImplicitEnablePartialReassemble, false);

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(FemModel, JacobianMatrix)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

//! use the AMGBackend solver by default
SET_TYPE_PROP(FemModel, LinearSolver, Dumux::AMGBackend<TypeTag>);

// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
SET_SCALAR_PROP(FemModel, LinearSolverResidualReduction, 1e-6);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(FemModel, LinearSolverMaxIterations, 250);

//! set number of equations of the mathematical model as default
SET_INT_PROP(FemModel, LinearSolverBlockSize, GET_PROP_VALUE(TypeTag, NumEq));

//! set number of maximum timestep divisions to 10
SET_INT_PROP(FemModel, ImplicitMaxTimeStepDivisions, 10);

//! Set the corresponding discretization method property
SET_PROP(FemModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::Fem;
};

//! The ip data
SET_TYPE_PROP(FemModel, FemIntegrationPointData, FemIntegrationPointData<TypeTag>);

//! Set the BaseLocalResidual to the FemLocalResidual
SET_TYPE_PROP(FemModel, BaseLocalResidual, FemLocalResidual<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(FemModel, JacobianAssembler, FemAssembler<TypeTag>);

//! The local jacobian operator
SET_TYPE_PROP(FemModel, LocalJacobian, FemLocalJacobian<TypeTag>);

//! By default we use an approximation order of 1
SET_INT_PROP(FemModel, FemBasisOrder, 1);

//! By default we use a quadratic integration rule
SET_INT_PROP(FemModel, FemQuadratureOrder, 2);

//! The global finite element basis
SET_TYPE_PROP(FemModel, FeBasis, Dune::Functions::PQkNodalBasis<typename GET_PROP_TYPE(TypeTag, GridView),
                                                                GET_PROP_VALUE(TypeTag, FemBasisOrder)>);

} // namespace Properties
} // namespace Dumux

#endif
