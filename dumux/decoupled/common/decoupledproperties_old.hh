// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_DECOUPLED_PROPERTIES_HH
#define DUMUX_DECOUPLED_PROPERTIES_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/basicproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>

/*!
 * \ingroup Sequential
 * \ingroup Properties
 */
/*!
 * \file
 * \brief Base file for properties related to sequential (decoupled) models
 */
namespace Dumux
{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! Create a type tag for all decoupled models
NEW_TYPE_TAG(DecoupledModel, INHERITS_FROM(ExplicitModel, LinearSolverTypeTag));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//! Property tag for types associated with the solution of the PDE.
//! This means vectors of primary variables, solution functions on the
//! grid, and elements, and shape functions.
NEW_PROP_TAG( SolutionTypes);
NEW_PROP_TAG( TransportSolutionType);
NEW_PROP_TAG( PrimaryVariables);

NEW_PROP_TAG( Grid); //!< The type of the DUNE grid
NEW_PROP_TAG( GridView); //!< The type of the grid view
NEW_PROP_TAG( AdaptiveGrid); //!< Defines if the grid is h-adaptive

NEW_PROP_TAG( Problem); //!< The type of the problem
NEW_PROP_TAG( Model); //!< The type of the discretizations
NEW_PROP_TAG( NumEq ); //!< Number of equations in the system of PDEs
NEW_PROP_TAG( NumPhases); //!< Number of phases in the system
NEW_PROP_TAG( NumComponents); //!< Number of components in the system
NEW_PROP_TAG( Variables); //!< The type of the container of global variables
NEW_PROP_TAG(TimeManager);  //!< Manages the simulation time
NEW_PROP_TAG(BoundaryTypes); //!< Stores the boundary types of a single degree of freedom
}
}

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/timemanager.hh>
#include <dumux/common/boundarytypes.hh>
#include<dumux/common/boundaryconditions.hh>

template<class TypeTag>
class VariableClass;

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_PROP(DecoupledModel, GridView)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

public:
    typedef typename Grid::LeafGridView type;
};

//! Disables Grid Adaptivity as standard
SET_BOOL_PROP(DecoupledModel, AdaptiveGrid, false);

//! Use the parent VariableClass
SET_TYPE_PROP(DecoupledModel, Variables, VariableClass<TypeTag>);

NEW_PROP_TAG(MaxIntersections);   //!< maximum number of intersections per element
SET_PROP(DecoupledModel, MaxIntersections)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    static const int dim = Grid::dimension;
public:
    typedef int type;
    static const int value = 2*dim;
};

/*!
 * \brief Specifies the types which are assoicated with a solution.
 *
 * This means shape functions, solution vectors, etc.
 */
SET_PROP(DecoupledModel, SolutionTypes)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, Variables) Variables;

    enum
    {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        maxIntersections = GET_PROP_VALUE(TypeTag, MaxIntersections)
    };

    template<int dim>
    struct VertexLayout
    {
        bool contains (Dune::GeometryType gt) const
        {   return gt.dim() == 0;}
    };

    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt) const
        {   return gt.dim() == dim;}
    };

public:
    /*!
     * \brief Mapper for the grid view's vertices.
     */
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, VertexLayout> VertexMapper;

    /*!
     * \brief Mapper for the grid view's elements.
     */
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, ElementLayout> ElementMapper;

    /*!
     * \brief The type of a solution at a fixed time.
     *
     * This defines the primary and secondary variable vectors at each degree of freedom.
     */
    typedef Dune::FieldVector<Scalar, numEq> PrimaryVariables;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarSolution;//!<type for vector of scalars
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numComponents> ComponentProperty;//!<type for vector of phase properties
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numPhases> PhaseProperty;//!<type for vector of phase properties
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numPhases> FluidProperty;//!<type for vector of fluid properties: Vector[element][phase]
    typedef Dune::BlockVector<Dune::FieldVector<Dune::FieldVector<Scalar, numPhases>, maxIntersections > > PhasePropertyElemFace;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef Dune::BlockVector<Dune::FieldVector<Dune::FieldVector<Scalar, dim>, maxIntersections > > DimVecElemFace;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars
};

SET_TYPE_PROP(DecoupledModel,  PrimaryVariables, typename GET_PROP(TypeTag, SolutionTypes)::PrimaryVariables);

/*!
 * \brief Default implementation for the Vector of the transportet quantity
 *
 * This type defines the data type of the transportet quantity. In case of a
 * immiscible 2p system, this would represent a vector holding the saturation
 * of one phase.
 */
SET_PROP(DecoupledModel, TransportSolutionType)
{
    private:
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionType;

    public:
    typedef typename SolutionType::ScalarSolution type;//!<type for vector of scalar properties
};

//! Set the default type for the time manager
SET_TYPE_PROP(DecoupledModel, TimeManager, Dumux::TimeManager<TypeTag>);

/*!
 * \brief Boundary types at a single degree of freedom.
 */
SET_PROP(DecoupledModel, BoundaryTypes)
{ private:
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
public:
    typedef Dumux::BoundaryTypes<numEq>  type;
};

//! set the default for the reduction of the initial residual
SET_PROP(DecoupledModel, LinearSolverResidualReduction)
{public:
    static constexpr double value = 1e-13;
};

//! set the default number of maximum iterations for the linear solver
SET_PROP(DecoupledModel, LinearSolverMaxIterations)
{public:
    static constexpr int value = 500;
};

//! set the default number of maximum iterations for the linear solver
SET_PROP(DecoupledModel, LinearSolverBlockSize)
{public:
    static constexpr int value = 1;
};
// \}

}
}

#endif
