// $Id$
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

#include <dumux/decoupled/2p/diffusion/mimetic/mimeticoperator.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticgroundwater.hh>

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
NEW_TYPE_TAG(DecoupledModel, INHERITS_FROM(ExplicitModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//! Property tag for types associated with the solution of the PDE.
//! This means vectors of primary variables, solution functions on the
//! grid, and elements, and shape functions.
NEW_PROP_TAG( SolutionTypes);
NEW_PROP_TAG( TransportSolutionType);

NEW_PROP_TAG( Grid); //!< The type of the DUNE grid
NEW_PROP_TAG( GridView); //!< The type of the grid view

NEW_PROP_TAG( Problem); //!< The type of the problem
NEW_PROP_TAG( Model); //!< The type of the discretizations
NEW_PROP_TAG( NumPhases); //!< Number of phases in the system
NEW_PROP_TAG( NumComponents); //!< Number of components in the system
NEW_PROP_TAG( Variables); //!< The type of the container of global variables
NEW_PROP_TAG( LocalStiffness); //!< The type of communication needed for the mimetic operator
}
}

#include <dumux/common/boundarytypes.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bvector.hh>

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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

public:
    typedef typename Grid::LeafGridView type;
};


//! Use the parent VariableClass
SET_TYPE_PROP(DecoupledModel, Variables, VariableClass<TypeTag>);

/*!
 * \brief Specifies the types which are assoicated with a solution.
 *
 * This means shape functions, solution vectors, etc.
 */
SET_PROP(DecoupledModel, SolutionTypes)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    enum
    {
        dim = GridView::dimension,
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents))
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
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarSolution;//!<type for vector of scalars
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numComponents> ComponentProperty;//!<type for vector of phase properties
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numPhases> PhaseProperty;//!<type for vector of phase properties
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numPhases> FluidProperty;//!<type for vector of fluid properties: Vector[element][phase]
    typedef Dune::BlockVector<Dune::FieldVector<Dune::FieldVector<Scalar, numPhases>, 2*dim > > PhasePropertyElemFace;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef Dune::BlockVector<Dune::FieldVector<Dune::FieldVector<Scalar, dim>, 2*dim > > DimVecElemFace;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars
};

/*!
 * \brief Default implementation for the Vector of the transportet quantity
 *
 * This type defines the data type of the transportet quantity. In case of a
 * immiscible 2p system, this would represent a vector holding the saturation
 * of one phase.
 */
SET_PROP_DEFAULT(TransportSolutionType)
{
    private:
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionType;

    public:
    typedef typename SolutionType::ScalarSolution type;//!<type for vector of scalar properties
};

// \}

}
}

#endif
