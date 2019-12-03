// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
#ifndef DUMUX_PRESSURE_PROPERTIES_HH
#define DUMUX_PRESSURE_PROPERTIES_HH

//Dune-includes
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include "properties.hh"
#include <dumux/linear/seqsolverbackend.hh>


/*!
 * \ingroup Sequential
 * \ingroup IMPETProperties
 */
/*!
 * \file
 * \brief Base file for properties related to sequential IMPET algorithms
 */
namespace Dumux
{
namespace Properties
{
/*!
 *
 * \brief General properties for sequential IMPET algorithms
 *
 * This class holds properties necessary for the sequential IMPET solution.
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
// Create new type tags
namespace TTag {
struct Pressure { using InheritsFrom = std::tuple<SequentialModel>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//Properties for linear solvers
template<class TypeTag, class MyTypeTag>
struct PressureRHSVector { using type = UndefinedProperty; };//!< Type of the right hand side vector given to the linear solver
template<class TypeTag, class MyTypeTag>
struct PressureSolutionVector { using type = UndefinedProperty; };//!Type of solution vector or pressure system
template<class TypeTag, class MyTypeTag>
struct VisitFacesOnlyOnce { using type = UndefinedProperty; }; //!< Indicates if faces are only regarded from one side
}
}

#include <dumux/porousmediumflow/sequential/cellcentered/velocitydefault.hh>

namespace Dumux
{
namespace Properties
{
//! Faces are only regarded from one side and not from both cells
template<class TypeTag>
struct VisitFacesOnlyOnce<TypeTag, TTag::Pressure> { static constexpr bool value = false; };

//Set defaults
template<class TypeTag>
struct PressureCoefficientMatrix<TypeTag, TTag::Pressure>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using MB = Dune::FieldMatrix<Scalar, 1, 1>;

public:
    using type = Dune::BCRSMatrix<MB>;
};
template<class TypeTag>
struct PressureRHSVector<TypeTag, TTag::Pressure>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Dune::BlockVector<Dune::FieldVector<Scalar, 1> >;
};

template<class TypeTag>
struct PressureSolutionVector<TypeTag, TTag::Pressure> { using type = typename GetProp<TypeTag, SolutionTypes>::ScalarSolution; };

// use the stabilized BiCG solver preconditioned by the ILU-0 by default
template<class TypeTag>
struct LinearSolver<TypeTag, TTag::Pressure> { using type = ILU0BiCGSTABBackend ; };

template<class TypeTag>
struct Velocity<TypeTag, TTag:: Pressure> { using type = FVVelocityDefault<TypeTag>; };

}
}

#endif
