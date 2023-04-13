// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Define traits for linear solvers.
 */
#ifndef DUMUX_LINEAR_SOLVER_TRAITS_HH
#define DUMUX_LINEAR_SOLVER_TRAITS_HH

#include <bitset>

#include <dune/istl/schwarz.hh>
#include <dune/istl/novlpschwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/grid/common/capabilities.hh>

#include <dumux/common/gridcapabilities.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

//! The implementation is specialized for the different discretizations
template<class GridGeometry, class DiscretizationMethod>
struct LinearSolverTraitsImpl;

/*!
 * \brief The type traits required for using the IstlFactoryBackend
 */
template<class GridGeometry>
using LinearSolverTraits = LinearSolverTraitsImpl<GridGeometry, typename GridGeometry::DiscretizationMethod>;

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

struct SeqLinearSolverTraits
{
    template<class Matrix, class Vector>
    using Sequential = SequentialSolverTraits<Matrix, Vector>;

    static constexpr bool canCommunicate = false;
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
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::Box>
: public LinearSolverTraitsBase<GridGeometry>
{
    using DofMapper = typename GridGeometry::VertexMapper;
    using Grid = typename GridGeometry::GridView::Traits::Grid;
    static constexpr int dofCodim = Grid::dimension;
    static constexpr bool canCommunicate = Dumux::Detail::canCommunicate<Grid, dofCodim>;

    template<class GridView>
    static bool isNonOverlapping(const GridView& gridView)
    { return gridView.overlapSize(0) == 0; }
};

template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::PQ1Bubble>
: public LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::Box>
{
    using Grid = typename GridGeometry::GridView::Traits::Grid;
    using DofMapper = typename GridGeometry::DofMapper;

    static constexpr int dofCodim = Grid::dimension;
    static constexpr std::bitset<Grid::dimension+1> dofCodims{ (1UL << Grid::dimension) + 1UL };

    static const DofMapper& dofMapper(const GridGeometry& gg)
    { return { gg.dofMapper() }; }
};

//! Cell-centered tpfa: use overlapping model
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::CCTpfa>
: public LinearSolverTraitsBase<GridGeometry>
{
    using DofMapper = typename GridGeometry::ElementMapper;
    using Grid = typename GridGeometry::GridView::Traits::Grid;
    static constexpr int dofCodim = 0;
    static constexpr bool canCommunicate = Dumux::Detail::canCommunicate<Grid, dofCodim>;

    template<class GridView>
    static bool isNonOverlapping(const GridView& gridView)
    { return false; }
};

//! Face-centered staggered: use overlapping model
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::FCStaggered>
: public LinearSolverTraitsBase<GridGeometry>
{
    class DofMapper
    {
    public:
        DofMapper(const typename GridGeometry::GridView& gridView)
        : gridView_(gridView) {}

        template<class Entity>
        auto index(const Entity& e) const
        { return gridView_.indexSet().index(e); }

        auto size() const
        { return gridView_.size(1); }

    private:
        typename GridGeometry::GridView gridView_;
    };

    static DofMapper dofMapper(const GridGeometry& gg)
    { return { gg.gridView() }; }

    using Grid = typename GridGeometry::GridView::Traits::Grid;
    static constexpr int dofCodim = 1;

    // TODO: see above for description of this workaround, remove second line if fixed upstream
    static constexpr bool canCommunicate =
        Dune::Capabilities::canCommunicate<Grid, dofCodim>::v
        || Dumux::Temp::Capabilities::canCommunicate<Grid, dofCodim>::v;

    template<class GridView>
    static bool isNonOverlapping(const GridView& gridView)
    {
        assert(gridView.overlapSize(0) > 0);
        return false;
    }
};

//! Face-centered diamond scheme: use overlapping or non-overlapping model depending on the grid
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::FCDiamond>
: public LinearSolverTraitsBase<GridGeometry>
{
    using DofMapper = typename GridGeometry::DofMapper;
    using Grid = typename GridGeometry::GridView::Traits::Grid;
    static constexpr int dofCodim = 1;
    static constexpr bool canCommunicate = Dumux::Detail::canCommunicate<Grid, dofCodim>;

    static const DofMapper& dofMapper(const GridGeometry& gg)
    { return { gg.dofMapper() }; }

    template<class GridView>
    static bool isNonOverlapping(const GridView& gridView)
    { return gridView.overlapSize(0) == 0; }
};

//! Cell-centered mpfa: use overlapping model
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::CCMpfa>
: public LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::CCTpfa> {};

//! staggered: use overlapping model
template<class GridGeometry>
struct LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::Staggered>
: public LinearSolverTraitsImpl<GridGeometry, DiscretizationMethods::CCTpfa> {};

} // end namespace Dumux

#endif
