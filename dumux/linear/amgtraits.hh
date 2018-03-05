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
 * \ingroup Linear
 *
 * \brief Define traits for the AMG backend.
 */
#ifndef DUMUX_AMG_TRAITS_HH
#define DUMUX_AMG_TRAITS_HH

#include <dune/istl/schwarz.hh>
#include <dune/istl/novlpschwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/grid/common/capabilities.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

//! The implementation is specialized for the different discretizations
template<class TypeTag, DiscretizationMethod discMethod> struct AmgTraitsImpl;

//! The type traits required for using the AMG backend
template<class TypeTag>
using AmgTraits = AmgTraitsImpl<TypeTag, GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod>;

//! NonoverlappingSolverTraits used by discretization with non-overlapping parallel model
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
template<class TypeTag>
struct AmgTraitsImpl<TypeTag, DiscretizationMethod::box>
{
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

    using DofMapper = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::VertexMapper;
};

//! OverlappingSolverTraits used by discretization with overlapping parallel model
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

//! Cell-centered tpfa: use the overlapping AMG
template<class TypeTag>
struct AmgTraitsImpl<TypeTag, DiscretizationMethod::cctpfa>
{
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

    using DofMapper = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::ElementMapper;
};

template<class TypeTag>
struct AmgTraitsImpl<TypeTag, DiscretizationMethod::ccmpfa>
: public AmgTraitsImpl<TypeTag, DiscretizationMethod::cctpfa> {};

} // end namespace Dumux

#endif // DUMUX_AMG_TRAITS_HH
