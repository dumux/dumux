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

#include <dumux/common/propertysystem.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxModel.
 */
namespace Dumux
{

namespace Properties
{
/*!
 * \addtogroup BoxModel
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the box-scheme
NEW_TYPE_TAG(BoxModel);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//!< Property tag for scalar vaslues
NEW_PROP_TAG(Scalar);

NEW_PROP_TAG(Grid);     //!< The type of the DUNE grid
NEW_PROP_TAG(GridView); //!< The type of the grid view

NEW_PROP_TAG(ReferenceElements); //!< DUNE reference elements to be used
NEW_PROP_TAG(FVElementGeometry); //! The type of the finite-volume geometry in the box scheme

NEW_PROP_TAG(Problem);       //!< The type of the problem
NEW_PROP_TAG(Model);         //!< The type of the discretization
NEW_PROP_TAG(NumEq);         //!< Number of equations in the system of PDEs
NEW_PROP_TAG(LocalResidual); //!< The type of the local residual function
NEW_PROP_TAG(LocalJacobian); //!< The type of the local jacobian operator

NEW_PROP_TAG(JacobianAssembler); //!< Assembles the global jacobian matrix
NEW_PROP_TAG(JacobianMatrix); //!< Type of the global jacobian matrix
NEW_PROP_TAG(BoundaryTypes); //!< Stores the boundary types of a single degree of freedom
NEW_PROP_TAG(ElementBoundaryTypes); //!< Stores the boundary types on an element

NEW_PROP_TAG(PrimaryVariables); //!< A vector of primary variables within a sub-control volume
NEW_PROP_TAG(SolutionVector); //!< Vector containing all primary variable vector of the grid
NEW_PROP_TAG(ElementSolutionVector); //!< A vector of primary variables within a sub-control volume

NEW_PROP_TAG(VolumeVariables);  //!< The secondary variables within a sub-control volume
NEW_PROP_TAG(ElementVolumeVariables); //!< The secondary variables of all sub-control volumes in an element
NEW_PROP_TAG(FluxVariables); //!< Data required to calculate a flux over a face

// high level simulation control
NEW_PROP_TAG(TimeManager);  //!< Manages the simulation time
NEW_PROP_TAG(NewtonMethod);     //!< The type of the newton method
NEW_PROP_TAG(NewtonController); //!< The type of the newton controller

// properties for the PDELab wrapper
NEW_PROP_TAG(LocalFEMSpace); //!< The local finite element space used for the finite element interpolation
NEW_PROP_TAG(ScalarGridFunctionSpace); //!< The used grid function space for a single finite element function
NEW_PROP_TAG(GridFunctionSpace); //!< The used grid function space
NEW_PROP_TAG(Constraints); //!< The constraints on the grid function space
NEW_PROP_TAG(ConstraintsTrafo); //!< The type of PDELab's constraints transformation
NEW_PROP_TAG(LocalOperator); //!< The type of the local operator used by PDELab
NEW_PROP_TAG(GridOperatorSpace); //!< The used grid operator space

//! Specify whether the jacobian matrix of the last iteration of a
//! time step should be re-used as the jacobian of the first iteration
//! of the next time step.
NEW_PROP_TAG(EnableJacobianRecycling);

//! Specify whether the jacobian matrix should be only reassembled for
//! elements where at least one vertex is above the specified
//! tolerance
NEW_PROP_TAG(EnablePartialReassemble);

// mappers from local to global indices
NEW_PROP_TAG(VertexMapper);
NEW_PROP_TAG(ElementMapper);
NEW_PROP_TAG(DofMapper);
}
}


#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>

#include "boxfvelementgeometry.hh"
#include "pdelabboxassembler.hh"
#include "pdelabboxistlvectorbackend.hh"
//#include <dumux/boxmodels/pdelab/boxdirichletconstraints.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/timemanager.hh>

namespace Dumux {

template<typename TypeTag>
class BoxFVElementGeometry;

template<typename TypeTag>
class BoxElementBoundaryTypes;

template<typename TypeTag>
class BoxElementVolumeVariables;

template<typename TypeTag>
class BoxLocalJacobian;

namespace PDELab {
template<typename TypeTag>
class BoxLocalOperator;

template<typename TypeTag>
class BoxAssembler;

template<typename TypeTag>
class BoxISTLVectorBackend;
};

namespace Properties {
//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for scalar values to double
SET_PROP_DEFAULT(Scalar)
{ typedef double type; };

//! Set the default type for the time manager
SET_PROP_DEFAULT(TimeManager)
{ typedef Dumux::TimeManager<TypeTag> type; };

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
    typedef Dune::GenericReferenceElement<CoordScalar, dim>  ReferenceElement;
};

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
SET_PROP(BoxModel, SolutionVector)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
public:
    typedef typename GridFunctionSpace::template VectorContainer<Scalar>::Type type;
};

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
SET_TYPE_PROP(BoxModel, JacobianAssembler, Dumux::PDELab::BoxAssembler<TypeTag>);

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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
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

// disable jacobian matrix recycling by default
SET_BOOL_PROP(BoxModel, EnableJacobianRecycling, false);
// disable partial reassembling by default
SET_BOOL_PROP(BoxModel, EnablePartialReassemble, false);

// \}

}
}

#endif
