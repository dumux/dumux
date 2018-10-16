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

#include <dumux/discretization/methods.hh>

// TODO: The following is a temporary solution to make the parallel AMG work
// for UGGrid. Once it is resolved upstream
// (https://gitlab.dune-project.org/core/dune-grid/issues/78),
// it should be guarded by a DUNE_VERSION macro and removed later.

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif // HAVE_UG

namespace Dumux {
namespace Temp {
namespace Capabilities {

template<class Grid, int codim>
struct canCommunicate
{
  static const bool v = false;
};

#if HAVE_UG
template<int dim, int codim>
struct canCommunicate<Dune::UGGrid<dim>, codim>
{
  static const bool v = true;
};
#endif // HAVE_UG

} // namespace Capabilities
} // namespace Temp
} // namespace Dumux
// end workaround

namespace Dumux {

//! The implementation is specialized for the different discretizations
template<class MType, class VType, class FVGridGeometry, DiscretizationMethod discMethod> struct AmgTraitsImpl;

//! The type traits required for using the AMG backend
template<class MType, class VType, class FVGridGeometry>
using AmgTraits = AmgTraitsImpl<MType, VType, FVGridGeometry, FVGridGeometry::discMethod>;

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
template<class Matrix, class Vector, class FVGridGeometry>
struct AmgTraitsImpl<Matrix, Vector, FVGridGeometry, DiscretizationMethod::box>
{
    using Grid = typename FVGridGeometry::GridView::Traits::Grid;
    enum {
        dofCodim = Grid::dimension,
        isNonOverlapping = true,
        // TODO: see above for description of this workaround, remove second line if fixed upstream
        isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
                     || Dumux::Temp::Capabilities::canCommunicate<Grid, dofCodim>::v
    };
    using MType = Matrix;
    using VType = Dune::BlockVector<Dune::FieldVector<typename Vector::block_type::value_type, Vector::block_type::dimension>>;
    using SolverTraits = NonoverlappingSolverTraits<MType, VType, isParallel>;
    using Comm = typename SolverTraits::Comm;
    using LinearOperator = typename SolverTraits::LinearOperator;
    using ScalarProduct = typename SolverTraits::ScalarProduct;
    using Smoother = typename SolverTraits::Smoother;

    using DofMapper = typename FVGridGeometry::VertexMapper;
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
template<class Matrix, class Vector, class FVGridGeometry>
struct AmgTraitsImpl<Matrix, Vector, FVGridGeometry, DiscretizationMethod::cctpfa>
{
    using Grid = typename FVGridGeometry::GridView::Traits::Grid;
    enum {
        dofCodim = 0,
        isNonOverlapping = false,
        // TODO: see above for description of this workaround, remove second line if fixed upstream
        isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
                     || Dumux::Temp::Capabilities::canCommunicate<Grid, dofCodim>::v
    };
    using MType = Matrix;
    using VType = Dune::BlockVector<Dune::FieldVector<typename Vector::block_type::value_type, Vector::block_type::dimension>>;
    using SolverTraits = OverlappingSolverTraits<MType, VType, isParallel>;
    using Comm = typename SolverTraits::Comm;
    using LinearOperator = typename SolverTraits::LinearOperator;
    using ScalarProduct = typename SolverTraits::ScalarProduct;
    using Smoother = typename SolverTraits::Smoother;

    using DofMapper = typename FVGridGeometry::ElementMapper;
};

template<class Matrix, class Vector, class FVGridGeometry>
struct AmgTraitsImpl<Matrix, Vector, FVGridGeometry, DiscretizationMethod::ccmpfa>
: public AmgTraitsImpl<Matrix, Vector, FVGridGeometry, DiscretizationMethod::cctpfa> {};

template<class Matrix, class Vector, class FVGridGeometry>
struct AmgTraitsImpl<Matrix, Vector, FVGridGeometry, DiscretizationMethod::godunov>
: public AmgTraitsImpl<Matrix, Vector, FVGridGeometry, DiscretizationMethod::cctpfa> {};

}// end namespace Dumux

#endif // DUMUX_AMG_TRAITS_HH
