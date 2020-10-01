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
/*!
 * \file
 * \ingroup Linear
 * \brief Define traits for linear solvers.
 */
#ifndef DUMUX_LINEAR_SOLVER_TRAITS_HH
#define DUMUX_LINEAR_SOLVER_TRAITS_HH

#include <dune/istl/schwarz.hh>
#include <dune/istl/novlpschwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/grid/common/capabilities.hh>

#include <dumux/discretization/method.hh>

// TODO: The following is a temporary solution to make the parallel AMG work
// for UGGrid. Once it is resolved upstream
// (https://gitlab.dune-project.org/core/dune-grid/issues/78),
// it should be guarded by a DUNE_VERSION macro and removed later.

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif // HAVE_UG

namespace Dumux::Temp::Capabilities {

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

} // namespace Dumux::Temp::Capabilities
// end workaround

namespace Dumux {

//! The implementation is specialized for the different discretizations
template<class GridGeometry, DiscretizationMethod discMethod>
struct LinearSolverTraitsImpl;

//! The type traits required for using the IstlFactoryBackend
template<class GridGeometry>
using LinearSolverTraits = LinearSolverTraitsImpl<GridGeometry, GridGeometry::discMethod>;

//! sequential solver traits
template<class MType, class VType>
struct SequentialSolverTraits
{
    using Matrix = MType;
    using Vector = VType;
    using LinearOperator = Dune::MatrixAdapter<MType, VType, VType>;
    using ScalarProduct = Dune::SeqScalarProduct<VType>;

    template<class SeqPreconditioner>
    using Preconditioner = SeqPreconditioner;
};

#if HAVE_MPI
template <class MType, class VType>
struct NonoverlappingSolverTraits
{
public:
    using Matrix = MType;
    using Vector = VType;
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    using LinearOperator = Dune::NonoverlappingSchwarzOperator<MType, VType, VType, Comm>;
    using ScalarProduct = Dune::NonoverlappingSchwarzScalarProduct<VType, Comm>;
    static constexpr bool isNonOverlapping = true;

    template<class SeqPreconditioner>
    using Preconditioner = Dune::NonoverlappingBlockPreconditioner<Comm, SeqPreconditioner>;
};

template <class MType, class VType>
struct OverlappingSolverTraits
{
public:
    using Matrix = MType;
    using Vector = VType;
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    using LinearOperator = Dune::OverlappingSchwarzOperator<MType, VType, VType, Comm>;
    using ScalarProduct = Dune::OverlappingSchwarzScalarProduct<VType, Comm>;
    static constexpr bool isNonOverlapping = false;

    template<class SeqPreconditioner>
    using Preconditioner = Dune::BlockPreconditioner<VType, VType, Comm, SeqPreconditioner>;
};
#endif

template<class GridGeometry>
struct LinearSolverTraitsBase
{
    using GridView = typename GridGeometry::GridView;
    using Grid = typename GridGeometry::GridView::Traits::Grid;

    template<class Matrix, class Vector>
    using Sequential = SequentialSolverTraits<Matrix, Vector>;

#if HAVE_MPI
    template<class Matrix, class Vector>
    using ParallelOverlapping = OverlappingSolverTraits<Matrix, Vector>;

    template<class Matrix, class Vector>
    using ParallelNonoverlapping = NonoverlappingSolverTraits<Matrix, Vector>;
#endif
};

//! Box: use overlapping or non-overlapping model depending on the grid
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethod::box>
: public LinearSolverTraitsBase<GridGeometry>
{
    using DofMapper = typename GridGeometry::VertexMapper;
    using Grid = typename GridGeometry::GridView::Traits::Grid;
    static constexpr int dofCodim = Grid::dimension;

    // TODO: see above for description of this workaround, remove second line if fixed upstream
    static constexpr bool canCommunicate =
             Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
             || Dumux::Temp::Capabilities::canCommunicate<Grid, dofCodim>::v;

    template<class GridView>
    static bool isNonOverlapping(const GridView& gridView)
    { return gridView.overlapSize(0) == 0; }
};

//! Cell-centered tpfa: use overlapping model
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethod::cctpfa>
: public LinearSolverTraitsBase<GridGeometry>
{
    using DofMapper = typename GridGeometry::ElementMapper;
    using Grid = typename GridGeometry::GridView::Traits::Grid;
    static constexpr int dofCodim = 0;

    // TODO: see above for description of this workaround, remove second line if fixed upstream
    static constexpr bool canCommunicate =
             Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
             || Dumux::Temp::Capabilities::canCommunicate<Grid, dofCodim>::v;

    template<class GridView>
    static bool isNonOverlapping(const GridView& gridView)
    { return false; }
};

//! Cell-centered mpfa: use overlapping model
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethod::ccmpfa>
: public LinearSolverTraitsImpl<GridGeometry, DiscretizationMethod::cctpfa> {};

//! staggered: use overlapping model TODO provide staggered-specific traits, combining overlapping/non-overlapping
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethod::staggered>
: public LinearSolverTraitsImpl<GridGeometry, DiscretizationMethod::cctpfa> {};

} // end namespace Dumux

#endif
