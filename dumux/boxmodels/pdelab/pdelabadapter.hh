/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
/*!
 * \file
 *
 * \brief This file provides wrappers which allow the Dumux box models to
 *        be used with dune pdelab
 */
#ifndef DUMUX_PDELAB_ADAPTER_HH
#define DUMUX_PDELAB_ADAPTER_HH

#if ! HAVE_DUNE_PDELAB
#error "DUNE-PDELab must be available in order to include this file!"
#endif

#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include "pdelabboxistlvectorbackend.hh"
#include "pdelabboxlocaloperator.hh"

namespace Dumux
{
namespace Properties
{
/*!
 * \addtogroup ModelCoupling
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for box problems which ought to be used in conjunction
//! with PDELab
NEW_TYPE_TAG(BoxPDELab);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//! Specifies the host grid
NEW_PROP_TAG(Grid);

//! Specifies the scalar grid function space used for sub-problems
NEW_PROP_TAG(ScalarGridFunctionSpace);

//! Specifies the grid function space used for sub-problems
NEW_PROP_TAG(GridFunctionSpace);

//! Specifies the grid operator used for sub-problems
NEW_PROP_TAG(GridOperator);

//! Specifies the grid operator space used for sub-problems
NEW_PROP_TAG(GridOperatorSpace);

//! Specifies the type of the constraints
NEW_PROP_TAG(Constraints);

//! Specifies the type of the constraints transformation
NEW_PROP_TAG(ConstraintsTrafo);

//! Specifies the local finite element space
NEW_PROP_TAG(LocalFEMSpace);

//! Specifies the local operator
NEW_PROP_TAG(LocalOperator);

}
}

namespace Dumux
{
namespace Properties
{
// set the grid functions space for the sub-models
SET_PROP(BoxPDELab, ScalarGridFunctionSpace)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalFEMSpace)) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Constraints)) Constraints;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
public:
    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, Constraints, Dumux::PDELab::BoxISTLVectorBackend<TypeTag> > type;
};

// set the grid functions space for the sub-models
SET_PROP(BoxPDELab, GridFunctionSpace)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ScalarGridFunctionSpace)) ScalarGridFunctionSpace;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
public:
    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, numEq, Dune::PDELab::GridFunctionSpaceBlockwiseMapper> type;
};

// set the grid function space for the sub-models
SET_TYPE_PROP(BoxPDELab, Constraints, Dune::PDELab::NoConstraints);

// set the grid functions space for the sub-models
SET_PROP(BoxPDELab, ConstraintsTrafo)
{private:
  typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
  typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
public:
    typedef typename GridFunctionSpace::template ConstraintsContainer<Scalar>::Type type;
};

// set the local operator used for submodels
SET_PROP(BoxPDELab, LocalOperator)
{ typedef Dumux::PDELab::BoxLocalOperator<TypeTag> type; };

// set the grid operator space used for submodels
// DEPRECATED: use GridOperator instead
SET_PROP(BoxPDELab, GridOperatorSpace)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConstraintsTrafo)) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;

#warning HACK: first line does not work. but why???
    //typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalOperator)) LocalOperator;
    typedef Dumux::PDELab::BoxLocalOperator<TypeTag> LocalOperator;

    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};

public:

    typedef Dune::PDELab::GridOperatorSpace<GridFunctionSpace,
        GridFunctionSpace,
        LocalOperator,
        ConstraintsTrafo,
        ConstraintsTrafo,
        Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
        true
        > type;
};

//! use the local FEM space associated with cubes by default
SET_PROP(BoxPDELab, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim>  type;
};

}
}

#endif
