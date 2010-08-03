// $Id: decoupledproperties.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_DECOUPLED_PROPERTIES_HH
#define DUMUX_DECOUPLED_PROPERTIES_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticoperator.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticgroundwater.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxScheme.
 */
namespace Dumux
{
namespace Properties
{
/*!
 * \addtogroup diffusion
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG( DecoupledModel);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//!< Property tag for scalar values
NEW_PROP_TAG( Scalar);

//! Property tag for types associated with the solution of the PDE.
//! This means vectors of primary variables, solution functions on the
//! grid, and elements, and shape functions.
NEW_PROP_TAG( SolutionTypes);

NEW_PROP_TAG( Grid); //!< The type of the DUNE grid
NEW_PROP_TAG( GridView); //!< The type of the grid view

NEW_PROP_TAG( ReferenceElements); //!< DUNE reference elements to be used

NEW_PROP_TAG( Problem); //!< The type of the problem
NEW_PROP_TAG( Model); //!< The type of the discretizations
NEW_PROP_TAG( NumPhases); //!< Number of phases in the system
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
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for scalar values to double
SET_PROP_DEFAULT(Scalar)
{   typedef double type;};

//! Use the leaf grid view if not defined otherwise
SET_PROP_DEFAULT(GridView)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

public:
    typedef typename Grid::LeafGridView type;
};

/*!
 * \brief Specify the reference elements which we ought to use.
 *
 * We use Dune::ReferenceElements by default (-> old entity
 * numbering).
 *
 * TODO: Some specialization if the grid only supports one kind of
 *       cells would be nice. this would be better fixed inside DUNE,
 *       though. something like:
 *       Dune::GenericReferenceElements<Dune::GeometryType<cube, dim> >
 */
SET_PROP_DEFAULT(ReferenceElements)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

    typedef typename Grid::ctype CoordScalar;
    static const int dim = Grid::dimension;

public:
    typedef Dune::GenericReferenceElements<CoordScalar, dim> Container;
    typedef Dune::GenericReferenceElements<CoordScalar, dim-1> ContainerFaces;
    typedef Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;
};

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_TYPE_PROP(DecoupledModel, Variables, VariableClass<TypeTag>);

SET_PROP_DEFAULT(LocalStiffness)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) Variables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

public:
    typedef MimeticGroundwaterEquationLocalStiffness<GridView,Scalar,Variables, Problem> type;
};

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
        dim = GridView::dimension, numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
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
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, numPhases> > PhaseProperty;//!<type for vector of phase properties
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, numPhases> > FluidProperty;//!<type for vector of fluid properties
    typedef Dune::BlockVector<Dune::FieldVector<Dune::FieldVector<Scalar, numPhases>, 2*dim > > PhasePropertyElemFace;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef Dune::BlockVector<Dune::FieldVector<Dune::FieldVector<Scalar, dim>, 2*dim > > DimVecElemFace;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars
};

// \}

}
}

#endif
