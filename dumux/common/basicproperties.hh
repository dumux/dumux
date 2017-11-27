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
 * \file
 *
 * \brief Defines a type tags and some fundamental properties for
 *        fully coupled and decoupled models
 */
#ifndef DUMUX_BASIC_PROPERTIES_HH
#define DUMUX_BASIC_PROPERTIES_HH

#include <dune/common/parametertree.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/common/balanceequationopts.hh>
#include <dumux/common/propertysystem.hh>
#include <dumux/common/properties.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/pointsource.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/io/defaultvtkoutputfields.hh>
#include <dumux/io/gridcreator.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for numeric models.
NEW_TYPE_TAG(NumericModel);

//! Set the default type of scalar values to double
SET_TYPE_PROP(NumericModel, Scalar, double);

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(NumericModel, GridView, typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! Set the default vector with size number of equations to a field vector
SET_TYPE_PROP(NumericModel, NumEqVector, Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar), GET_PROP_VALUE(TypeTag, NumEq)>);

//! Set the default primary variable vector to a vector of size of number of equations
SET_TYPE_PROP(NumericModel, PrimaryVariables, typename GET_PROP_TYPE(TypeTag, NumEqVector));

//! use the global group as default for the model's parameter group
SET_STRING_PROP(NumericModel, ModelParameterGroup, "");

//! do not specific any model-specific default parameters here
SET_PROP(NumericModel, ModelDefaultParameters)
{
    static void defaultParams(Dune::ParameterTree& tree, const std::string& group = "") { }
};

//! Use the DgfGridCreator by default
SET_TYPE_PROP(NumericModel, GridCreator, GridCreator<TypeTag>);

//! Use the minimal point source implementation as default
SET_TYPE_PROP(NumericModel, PointSource, PointSource<TypeTag>);

//! Use the point source helper using the bounding box tree as a default
SET_TYPE_PROP(NumericModel, PointSourceHelper, BoundingBoxTreePointSourceHelper<TypeTag>);

//! Set default output level to 0 -> only primary variables are added to output
SET_INT_PROP(NumericModel, VtkOutputLevel, 0);

//! Set the default to a function throwing a NotImplemented error
SET_TYPE_PROP(NumericModel, VtkOutputFields, DefaultVtkOutputFields);

//! Mapper for the grid view's vertices.
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
SET_TYPE_PROP(NumericModel,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView)>);
#else
SET_TYPE_PROP(NumericModel,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);
#endif

//! Mapper for the grid view's elements.
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
SET_TYPE_PROP(NumericModel,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView)>);
#else
SET_TYPE_PROP(NumericModel,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);
#endif

//! The type of a solution for the whole grid at a fixed time
SET_TYPE_PROP(NumericModel, SolutionVector, Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(NumericModel, JacobianMatrix)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

//! Boundary types at a single degree of freedom
SET_TYPE_PROP(NumericModel, BoundaryTypes, BoundaryTypes<GET_PROP_VALUE(TypeTag, NumEq)>);

//! Set the default class for the balance equation options
SET_TYPE_PROP(NumericModel, BalanceEqOpts, BalanceEquationOptions<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
