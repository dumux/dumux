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

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/pointsource.hh>
#include <dumux/io/defaultvtkoutputfields.hh>
#include <dumux/io/gridcreator.hh>

namespace Dumux
{
namespace Properties
{
///////////////////////////////////
// Type tag definitions:
///////////////////////////////////

//! Type tag for all models.
NEW_TYPE_TAG(NumericModel);

/////////////////////////////////////////////
// Property names which are always available:
/////////////////////////////////////////////

//! Property to specify the type of scalar values.
NEW_PROP_TAG(Scalar);

//! Property which provides a Dune::ParameterTree.
NEW_PROP_TAG(ParameterTree);

//! Property which defines the group that is queried for parameters by default
NEW_PROP_TAG(ModelParameterGroup);

//! Property which defines the group that is queried for parameters by default
NEW_PROP_TAG(ModelDefaultParameters);

//! Property which defines the group that is queried for grid (creator) parameters by default
NEW_PROP_TAG(GridParameterGroup);

//! Property which provides a GridCreator (manages grids)
NEW_PROP_TAG(GridCreator);

//! The DUNE grid type
NEW_PROP_TAG(Grid);

//! The number of equations to solve (equal to number of primary variables)
NEW_PROP_TAG(NumEq);

//! A vector of primary variables
NEW_PROP_TAG(PrimaryVariables);

//! A vector of size number equations that can be used for
//! Neumann fluxes, sources, residuals, ...
NEW_PROP_TAG(NumEqVector);

//! The type of the grid view according to the grid type
NEW_PROP_TAG(GridView);

//! Property to specify the type of a problem which has to be solved
NEW_PROP_TAG(Problem);

//! Property to specify the name of the problem
NEW_PROP_TAG(ProblemName);

//! Property defining the type of the model which is used to solve the problem
NEW_PROP_TAG(Model);

//! Property defining the type of point source used
NEW_PROP_TAG(PointSource);

//! Property defining the class that computes which sub control volume point sources belong to
NEW_PROP_TAG(PointSourceHelper);

//! The grid variables object managing variable data on the grid
NEW_PROP_TAG(GridVariables);

/*!
 * \brief Specify the maximum size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
NEW_PROP_TAG(TimeLoopMaxTimeStepSize);
NEW_PROP_TAG(TimeManagerMaxTimeStepSize);

//! the maximum allowed number of timestep divisions for the
//! Newton solver
NEW_PROP_TAG(TimeLoopMaxTimeStepDivisions);

//! Property to define the output level
NEW_PROP_TAG(VtkOutputLevel);

//! A class helping models to define default vtk output parameters
NEW_PROP_TAG(VtkOutputFields);

///////////////////////////////////
// Default values for properties:
///////////////////////////////////

//! Set the default problem name to dumuxsim
SET_STRING_PROP(NumericModel, ProblemName, "dumuxsim");

//! Set the default type of scalar values to double
SET_TYPE_PROP(NumericModel, Scalar, double);

//! Set the default vector with size number of equations to a field vector
SET_TYPE_PROP(NumericModel, NumEqVector, Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                           GET_PROP_VALUE(TypeTag, NumEq) >);

//! Set the default primary variable vector to a vector of size of number of equations
SET_TYPE_PROP(NumericModel, PrimaryVariables, typename GET_PROP_TYPE(TypeTag, NumEqVector));

//! use an unlimited time step size by default
SET_SCALAR_PROP(NumericModel, TimeLoopMaxTimeStepSize, std::numeric_limits<typename GET_PROP_TYPE(TypeTag,Scalar)>::max());
SET_SCALAR_PROP(NumericModel, TimeManagerMaxTimeStepSize, std::numeric_limits<typename GET_PROP_TYPE(TypeTag,Scalar)>::max());

//! set number of maximum timestep divisions to 10
SET_INT_PROP(NumericModel, TimeLoopMaxTimeStepDivisions, 10);

//! Set the ParameterTree property
SET_PROP(NumericModel, ParameterTree)
{
    typedef Dune::ParameterTree type;

    static Dune::ParameterTree &tree()
    {
        static Dune::ParameterTree obj_;
        return obj_;
    }

    static Dune::ParameterTree &compileTimeParams()
    {
        static Dune::ParameterTree obj_;
        return obj_;
    }

    static Dune::ParameterTree &runTimeParams()
    {
        static Dune::ParameterTree obj_;
        return obj_;
    }

    static Dune::ParameterTree &deprecatedRunTimeParams()
    {
        static Dune::ParameterTree obj_;
        return obj_;
    }

    static Dune::ParameterTree &unusedNewRunTimeParams()
    {
        static Dune::ParameterTree obj_;
        return obj_;
    }
};

//! use the global group as default for the model's parameter group
SET_STRING_PROP(NumericModel, ModelParameterGroup, "");

//! do not specific any model-specific default parameters here
SET_PROP(NumericModel, ModelDefaultParameters)
{
    static void defaultParams(Dune::ParameterTree& tree, const std::string& group = "") { }
};

SET_STRING_PROP(NumericModel, GridParameterGroup, "Grid");

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

} // namespace Properties
} // namespace Dumux

#endif
