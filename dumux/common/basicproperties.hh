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

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/pointsource.hh>
#include <dumux/io/gridcreator.hh>
#include <dumux/io/vtkmultiwriter.hh>

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

//! Property which defines the group that is queried for grid (creator) parameters by default
NEW_PROP_TAG(GridParameterGroup);

//! Property which provides a GridCreator (manages grids)
NEW_PROP_TAG(GridCreator);

//! The DUNE grid type
NEW_PROP_TAG(Grid);

//! The type of the grid view according to the grid type
NEW_PROP_TAG(GridView);

//! Property to specify the type of a problem which has to be solved
NEW_PROP_TAG(Problem);

//! Property defining the type of the model which is used to solve the problem
NEW_PROP_TAG(Model);

//! Property defining the type of point source used
NEW_PROP_TAG(PointSource);

/*!
 * \brief Specify the maximum size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
NEW_PROP_TAG(TimeManagerMaxTimeStepSize);

//! Property to define the output level
NEW_PROP_TAG(VtkOutputLevel);

//! the type of VTK Writer to be used, i.e. ascii or binary (Dune::VTK::appendraw) format
NEW_PROP_TAG(VtkMultiWriter);

///////////////////////////////////
// Default values for properties:
///////////////////////////////////

//! Set the default type of scalar values to double
SET_TYPE_PROP(NumericModel, Scalar, double);

//! use an unlimited time step size by default
SET_SCALAR_PROP(NumericModel, TimeManagerMaxTimeStepSize, std::numeric_limits<typename GET_PROP_TYPE(TypeTag,Scalar)>::max());

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

//! use the Grid group as default for the grid parameter group
SET_STRING_PROP(NumericModel, GridParameterGroup, "Grid");

//! Use the DgfGridCreator by default
SET_TYPE_PROP(NumericModel, GridCreator, Dumux::GridCreator<TypeTag>);

//! Use the minimal point source implementation as default
SET_TYPE_PROP(NumericModel, PointSource, Dumux::PointSource<TypeTag>);

//! Set default output level to 0 -> only primary variables are added to output
SET_INT_PROP(NumericModel, VtkOutputLevel, 0);

//! set the VtkMultiWriter such that it uses the ascii format by default
SET_PROP(NumericModel, VtkMultiWriter)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef typename Dumux::VtkMultiWriter<GridView> type;
};

} // namespace Properties
} // namespace Dumux

#endif
