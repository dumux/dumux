// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
#ifndef DUMUX_BOX_PROPERTIES_HH
#define DUMUX_BOX_PROPERTIES_HH

#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include <dumux/fvgeometry/fvelementgeometry.hh>

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>

#include <dumux/auxiliary/properties.hh>
#include <dumux/boxmodels/tags.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxScheme.
 */
namespace Dune
{
namespace Properties
{
/*!
 * \addtogroup BoxScheme
 */
// \{
 
//! Set the default for the FVElementGeometry
SET_PROP(BoxScheme, FVElementGeometry)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

public:
    typedef Dune::FVElementGeometry<Grid>  type;
};

//! use the plain newton method for the box scheme by default
SET_PROP(BoxScheme, NewtonMethod)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))  Model;

public:
    typedef Dune::NewtonMethod<Model> type;
};

//! use the plain newton controller for the box scheme by default
SET_PROP(BoxScheme, NewtonController)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;

public:
    typedef Dune::NewtonController<NewtonMethod> type;
};

/*!
 * \brief Specifies the types which are assoicated with a solution.
 *
 * This means shape functions, solution vectors, etc.
 */
SET_PROP(BoxScheme, SolutionTypes)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))    Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))  GridView;
    typedef typename GridView::Grid                          Grid;
    typedef typename Grid::ctype                             CoordScalar;

    enum {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))
    };
    
public:
    //! A solution function. This is a function with the same domain as the grid. 
    typedef Dune::P1Function<GridView, Scalar, GridView, numEq>       SolutionFunction;

    /*!
     * \brief The type which maps an entity with an attached degree of
     *        freedom to an index in the solution.
     */
    typedef typename SolutionFunction::VM                              DofEntityMapper;

    /*!
     * \brief The type of a solution at a fixed time. 
     *
     * This is the representation of the solution function and defines
     * a primary variable vector at each degree of freedom.
     */
    typedef typename SolutionFunction::RepresentationType             Solution;
    /*!
     * \brief A vector of primary variables.
     */
    typedef typename SolutionFunction::BlockType                      PrimaryVarVector;

    /*!
     * \brief The solution for a single finite element.
     */
    typedef Dune::BlockVector<PrimaryVarVector>                       SolutionOnElement;
    
    /*!
     * \brief Vector of boundary types at a degree of freedom.
     */
    typedef Dune::FieldVector<Dune::BoundaryConditions::Flags, numEq> BoundaryTypeVector;
    
    /*!
     * \brief The shape functions used by the SolutionFunction
     */
    typedef Dune::LagrangeShapeFunctions<CoordScalar, Scalar, dim>    ShapeFunctions;
    /*!
     * \brief The container for shape functions within an finite element.
     */
    typedef Dune::LagrangeShapeFunctionSet<CoordScalar, Scalar, dim>  ShapeFunctionSet;
    /*!
     * \brief Assembler for the global jacobian matrix.
     */
    typedef Dune::P1OperatorAssembler<Grid, Scalar, GridView, GridView, numEq>  JacobianAssembler;
};

// \}

}
}

#endif
