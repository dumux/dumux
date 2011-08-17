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

#include "boxassembler.hh"
#include "boxfvelementgeometry.hh"
#include "boxelementboundarytypes.hh"
#include "boxlocaljacobian.hh"
#include "boxelementvolumevariables.hh"
#include "boxvolumevariables.hh"

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/timemanager.hh>

#include "boxproperties.hh"

#include <limits>

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
SET_PROP(BoxModel, SolutionVector)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)) };
public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, numEq> > type;
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

//! use an unlimited time step size by default
#if 0
// requires GCC 4.6 and above to call the constexpr function of
// numeric_limits
SET_SCALAR_PROP(BoxModel, MaxTimeStepSize, std::numeric_limits<Scalar>::infinity());
#else
SET_SCALAR_PROP(BoxModel, MaxTimeStepSize, 1e100);
#endif

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(BoxModel, NumericDifferenceMethod, +1);

//! do not use hints by default
SET_BOOL_PROP(BoxModel, EnableHints, false);

// disable jacobian matrix recycling by default
SET_BOOL_PROP(BoxModel, EnableJacobianRecycling, false);
// disable partial reassembling by default
SET_BOOL_PROP(BoxModel, EnablePartialReassemble, false);

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

// use the stabilized BiCG solver preconditioned by the ILU-0 by default
SET_TYPE_PROP(BoxModel, LinearSolver, Dumux::BoxBiCGStabILU0Solver<TypeTag> );

// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
SET_PROP(BoxModel, LinearSolverResidualReduction)
{public:
    static constexpr double value = 1e-6;
};

//! set the default number of maximum iterations for the linear solver
SET_PROP(BoxModel, LinearSolverMaxIterations)
{public:
    static constexpr int value = 250;
};

//! set number of equations of the mathematical model as default
SET_PROP_DEFAULT(LinearSolverBlockSize)
{public:
    static constexpr int value = GET_PROP_VALUE(TypeTag, PTAG(NumEq));
};

} // namespace Properties
} // namespace Dumux

#endif
