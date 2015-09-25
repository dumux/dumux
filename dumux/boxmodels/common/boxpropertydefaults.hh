// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
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
 * \brief Default properties for box models
 */
#ifndef DUMUX_BOX_PROPERTY_DEFAULTS_HH
#define DUMUX_BOX_PROPERTY_DEFAULTS_HH

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>

#include "pdelabboxistlvectorbackend.hh"
#include "boxassembler.hh"
#include "boxfvelementgeometry.hh"
#include "boxelementboundarytypes.hh"
#include "boxlocaljacobian.hh"
#include "boxelementvolumevariables.hh"
#include "boxvolumevariables.hh"

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#endif // HAVE_DUNE_PDELAB

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/timemanager.hh>

#include "boxproperties.hh"

namespace Dumux {

namespace Properties {
//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for the time manager
SET_TYPE_PROP(BoxModel, TimeManager, Dumux::TimeManager<TypeTag>);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_PROP(BoxModel, GridView)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

public:
    typedef typename Grid::LeafGridView type;
};

//! Set the default for the FVElementGeometry
SET_PROP(BoxModel, FVElementGeometry)
{
    typedef Dumux::BoxFVElementGeometry<TypeTag>  type;
};

//! Set the default for the ElementBoundaryTypes
SET_PROP(BoxModel, ElementBoundaryTypes)
{ typedef Dumux::BoxElementBoundaryTypes<TypeTag> type; };

//! use the plain newton method for the box scheme by default
SET_PROP(BoxModel, NewtonMethod)
{public:
    typedef Dumux::NewtonMethod<TypeTag> type;
};

//! use the plain newton controller for the box scheme by default
SET_PROP(BoxModel, NewtonController)
{public:
    typedef Dumux::NewtonController<TypeTag> type;
};

//! Mapper for the grid view's vertices.
SET_PROP(BoxModel, VertexMapper)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    template<int dim>
    struct VertexLayout {
        bool contains (Dune::GeometryType gt) const
        { return gt.dim() == 0; }
    };
public:
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, VertexLayout> type;
};

//! Mapper for the grid view's elements.
SET_PROP(BoxModel, ElementMapper)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    template<int dim>
    struct ElementLayout {
        bool contains (Dune::GeometryType gt) const
        { return gt.dim() == dim; }
    };
public:
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, ElementLayout> type;
};

//! Mapper for the degrees of freedoms.
SET_PROP(BoxModel, DofMapper)
{ typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) type; };

//! The local jacobian operator for the box scheme
SET_TYPE_PROP(BoxModel, LocalJacobian, Dumux::BoxLocalJacobian<TypeTag>);

/*!
 * \brief The type of a solution for the whole grid at a fixed time.
 */
#if HAVE_DUNE_PDELAB
SET_PROP(BoxModel, SolutionVector)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
public:
    typedef typename GridFunctionSpace::template VectorContainer<Scalar>::Type type;
};
#else
SET_PROP(BoxModel, SolutionVector)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)) };
public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, numEq> > type;
};
#endif

/*!
 * \brief The type of a solution for a whole element.
 */
SET_PROP(BoxModel, ElementSolutionVector)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
public:
    typedef Dune::BlockVector<PrimaryVariables> type;
};

/*!
 * \brief A vector of primary variables.
 */
SET_PROP(BoxModel, PrimaryVariables)
{ typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector))::block_type type; };

/*!
 * \brief The volume variable class.
 *
 * This should almost certainly be overloaded by the model...
 */
SET_TYPE_PROP(BoxModel, VolumeVariables, Dumux::BoxVolumeVariables<TypeTag>);

/*!
 * \brief An array of secondary variable containers.
 */
SET_TYPE_PROP(BoxModel, ElementVolumeVariables, Dumux::BoxElementVolumeVariables<TypeTag>);

/*!
 * \brief Boundary types at a single degree of freedom.
 */
SET_PROP(BoxModel, BoundaryTypes)
{ private:
    enum { numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)) };
public:
    typedef Dumux::BoundaryTypes<numEq>  type;
};

/*!
 * \brief Assembler for the global jacobian matrix.
 */
SET_TYPE_PROP(BoxModel, JacobianAssembler, Dumux::BoxAssembler<TypeTag>);


//! use central differences to calculate the jacobian by default
SET_INT_PROP(BoxModel, NumericDifferenceMethod, 0);

// disable jacobian matrix recycling by default
SET_BOOL_PROP(BoxModel, EnableJacobianRecycling, false);
// disable partial reassembling by default
SET_BOOL_PROP(BoxModel, EnablePartialReassemble, false);
// disable time-step ramp up by default
SET_BOOL_PROP(BoxModel, EnableTimeStepRampUp, false);

#if ! HAVE_DUNE_PDELAB

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(BoxModel, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)) };
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
public:
    typedef typename Dune::BCRSMatrix<MatrixBlock> type;
};
#endif


#if HAVE_DUNE_PDELAB
//! Extract the type of a global jacobian matrix from the solution types
SET_PROP(BoxModel, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridOperatorSpace)) GridOperatorSpace;
public:
    typedef typename GridOperatorSpace::template MatrixContainer<Scalar>::Type type;
};

//! use the local FEM space associated with cubes by default
SET_PROP(BoxModel, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim>  type;
};

SET_PROP(BoxModel, GridFunctionSpace)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalFEMSpace)) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
public:
    //typedef MyBoxConstraints Constraints;
    typedef Dune::PDELab::NoConstraints Constraints;

    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, Constraints, Dumux::PDELab::BoxISTLVectorBackend<TypeTag> >
        ScalarGridFunctionSpace;

    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, numEq, Dune::PDELab::GridFunctionSpaceBlockwiseMapper>
        type;

    typedef typename type::template ConstraintsContainer<Scalar>::Type
        ConstraintsTrafo;
};

SET_PROP(BoxModel, ConstraintsTrafo)
{ typedef typename GET_PROP(TypeTag, PTAG(GridFunctionSpace))::ConstraintsTrafo type; };
SET_PROP(BoxModel, Constraints)
{ typedef typename GET_PROP(TypeTag, PTAG(GridFunctionSpace))::Constraints type; };
SET_PROP(BoxModel, ScalarGridFunctionSpace)
{ typedef typename GET_PROP(TypeTag, PTAG(GridFunctionSpace))::ScalarGridFunctionSpace type; };

SET_PROP(BoxModel, GridOperatorSpace)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConstraintsTrafo)) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};

public:
    typedef Dumux::PDELab::BoxLocalOperator<TypeTag> LocalOperator;

    typedef Dune::PDELab::GridOperatorSpace<GridFunctionSpace,
        GridFunctionSpace,
        LocalOperator,
        ConstraintsTrafo,
        ConstraintsTrafo,
        Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
        true
        > type;
};


SET_PROP(BoxModel, LocalOperator)
{ typedef typename GET_PROP(TypeTag, PTAG(GridOperatorSpace))::LocalOperator type; };

#endif // HAVE_DUNE_PDELAB

} // namespace Properties
} // namespace Dumux

#endif
