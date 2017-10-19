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
 * \file
 * \ingroup Properties
 * \ingroup Linear
 *
 * \brief Defines some fundamental properties for the AMG backend.
 */
#ifndef DUMUXAMGPROPERTIES_HH
#define DUMUXAMGPROPERTIES_HH

#include <dune/istl/schwarz.hh>
#include <dune/istl/novlpschwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/common/version.hh>

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/properties.hh>
#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include "linearsolverproperties.hh"

namespace Dumux
{

// Forward declaration for the property definitions
template <class TypeTag> class AMGBackend;

namespace Properties
{
//! The type traits required for using the AMG backend
NEW_PROP_TAG(AmgTraits);

template <class MType, class VType, bool isParallel>
class NonoverlappingSolverTraits
{
public:
    using Comm = Dune::Amg::SequentialInformation;
    using LinearOperator = Dune::MatrixAdapter<MType,VType,VType>;
    using ScalarProduct = Dune::SeqScalarProduct<VType>;
    using Smoother = Dune::SeqSSOR<MType,VType, VType>;
};

#if HAVE_MPI
template <class MType, class VType>
class NonoverlappingSolverTraits<MType, VType, true>
{
public:
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>,int>;
    using LinearOperator = Dune::NonoverlappingSchwarzOperator<MType,VType, VType,Comm>;
    using ScalarProduct = Dune::NonoverlappingSchwarzScalarProduct<VType,Comm>;
    using Smoother = Dune::NonoverlappingBlockPreconditioner<Comm,Dune::SeqSSOR<MType,VType, VType> >;
};
#endif

//! Box: use the non-overlapping AMG
SET_PROP(BoxModel, AmgTraits)
{
public:
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = Grid::dimension,
        isNonOverlapping = true,
        isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
    };
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MType = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> >;
    using VType = Dune::BlockVector<Dune::FieldVector<Scalar,numEq> >;
    using SolverTraits = NonoverlappingSolverTraits<MType, VType, isParallel>;
    using Comm = typename SolverTraits::Comm;
    using LinearOperator = typename SolverTraits::LinearOperator;
    using ScalarProduct = typename SolverTraits::ScalarProduct;
    using Smoother = typename SolverTraits::Smoother;

    using DofMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
};

template <class MType, class VType, bool isParallel>
class OverlappingSolverTraits
{
public:
    using Comm = Dune::Amg::SequentialInformation;
    using LinearOperator = Dune::MatrixAdapter<MType,VType,VType>;
    using ScalarProduct = Dune::SeqScalarProduct<VType>;
    using Smoother = Dune::SeqSSOR<MType,VType, VType>;
};

#if HAVE_MPI
template <class MType, class VType>
class OverlappingSolverTraits<MType, VType, true>
{
public:
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>,int>;
    using LinearOperator = Dune::OverlappingSchwarzOperator<MType,VType, VType,Comm>;
    using ScalarProduct = Dune::OverlappingSchwarzScalarProduct<VType,Comm>;
    using Smoother = Dune::BlockPreconditioner<VType,VType,Comm,Dune::SeqSSOR<MType,VType, VType> >;
};
#endif

//! Cell-centered: use the overlapping AMG
SET_PROP(CCModel, AmgTraits)
{
public:
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = 0,
        isNonOverlapping = false,
        isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
    };
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MType = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> >;
    using VType = Dune::BlockVector<Dune::FieldVector<Scalar,numEq> >;
    using SolverTraits = OverlappingSolverTraits<MType, VType, isParallel>;
    using Comm = typename SolverTraits::Comm;
    using LinearOperator = typename SolverTraits::LinearOperator;
    using ScalarProduct = typename SolverTraits::ScalarProduct;
    using Smoother = typename SolverTraits::Smoother;

    using DofMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);
};

//! Sequential model: use the overlapping AMG
SET_PROP(SequentialModel, AmgTraits)
{
public:
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    enum {
        numEq = JacobianMatrix::block_type::rows,
        dofCodim = 0,
        isNonOverlapping = false,
        isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
    };
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MType = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> >;
    using VType = Dune::BlockVector<Dune::FieldVector<Scalar,numEq> >;
    using SolverTraits = OverlappingSolverTraits<MType, VType, isParallel>;
    using Comm = typename SolverTraits::Comm;
    using LinearOperator = typename SolverTraits::LinearOperator;
    using ScalarProduct = typename SolverTraits::ScalarProduct;
    using Smoother = typename SolverTraits::Smoother;

    using DofMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);
};

} // namespace Properties
} // namespace Dumux
#endif
