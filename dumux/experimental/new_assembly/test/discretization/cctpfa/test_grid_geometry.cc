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
 * \ingroup CCTpfaDiscretization
 * \brief Test for the `CCTpfaGridGeometry` class.
 */
#include <iostream>
#include <cmath>

#include <dune/grid/yaspgrid.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/timer.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometry.hh>

template<typename GridGeometry>
auto tolerance(const GridGeometry& gg)
{
    typename GridGeometry::GridView::ctype tol = 1e6;
    using std::min;
    for (int dim = 0; dim < GridGeometry::GridView::dimension; ++dim)
        tol = min(tol, gg.bBoxMax()[dim] - gg.bBoxMin()[dim]);
    return tol*1e-10;
}

template<std::floating_point ctype>
bool equal(const ctype a, const ctype b, const ctype tolerance)
{ return Dune::FloatCmp::eq(a, b, tolerance); }

template<typename ctype, int dim>
bool equal(const Dune::FieldVector<ctype, dim>& p0,
           const Dune::FieldVector<ctype, dim>& p1,
           const ctype tolerance)
{
    return std::ranges::equal(p0, p1, [&] (const ctype a, const ctype b) {
        return equal(a, b, tolerance);
    });
}

template<typename GridGeometry, typename ctype, int dim>
bool isOnDomainBoundary(const GridGeometry& gridGeometry,
                        const Dune::FieldVector<ctype, dim>& pos)
{
    const auto tol = tolerance(gridGeometry);
    const auto& min = gridGeometry.bBoxMin();
    const auto& max = gridGeometry.bBoxMax();

    bool onBoundary = false;
    for (int i = 0; i < dim; ++i)
    {
        if (pos[i] < min[i] - tol) return false;
        if (pos[i] > max[i] + tol) return false;
        if (equal(pos[i], min[i], tol) || equal(pos[i], max[i], tol))
            onBoundary = true;
    }
    return onBoundary;
}

template<typename GridView,
         bool caching = false,
         typename Traits = Dumux::DefaultCCTpfaGridGeometryTraits<GridView>>
void testCCTpfaGridGeometry(const GridView& gridView, const std::size_t expectedNumScvf)
{
    Dune::Timer timer;
    Dumux::CCTpfaGridGeometry<GridView, caching, Traits> gridGeometry{gridView};
    timer.stop();
    std::cout << "Construction took " << timer.elapsed() << " seconds\n";

    std::size_t numScv = 0;
    std::size_t numScvf = 0;
    for (const auto& element : elements(gridView))
    {
        const auto fvGeometry = localView(gridGeometry).bindElement(element);
        for ([[maybe_unused]] const auto& scv : scvs(fvGeometry))
            numScv++;
        for ([[maybe_unused]] const auto& scvf : scvfs(fvGeometry))
            numScvf++;
    }

    if (gridGeometry.numDofs() != gridView.size(0))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected numDofs()");
    if (numScv != gridView.size(0))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected numScv()");
    if (numScvf != expectedNumScvf)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected numScvf()");

    const auto tol = tolerance(gridGeometry);
    const auto testLocalView = [&] <typename LocalView> (const auto& element,
                                                         const LocalView& localView,
                                                         bool withStencil) {
        static_assert(Dumux::Concepts::CCGridGeometryLocalView<LocalView>);

        int count = 0;
        int scvDofIdx = -1;
        for (const auto& scv : scvs(localView))
        {
            scvDofIdx = static_cast<std::size_t>(scv.dofIndex());
            if (!equal(scv.center(), element.geometry().center(), tol))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scv center");
            if (!equal(scv.center(), localView.element().geometry().center(), tol))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected localView-element center");
            if (!equal(element.geometry().center(), localView.geometry(scv).center(), tol))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scv geometry center");
            if (!Dune::FloatCmp::eq(scv.volume(), element.geometry().volume(), tol))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scv volume");
            if (!Dune::FloatCmp::eq(localView.geometry(scv).volume(), element.geometry().volume(), tol))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scv geometry volume");
            if (localView.indexInElement(scv) != count)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scv indexin element");
            count++;
        }

        if (count != 1)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvs");

        if (withStencil)
        {
            std::vector<std::size_t> neighborScvIndices;
            for (const auto& scv : neighborScvs(localView))
                if (neighborScvIndices.push_back(scv.dofIndex()); scv.dofIndex() == scvDofIdx)
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected neighbor scv dof");

            std::vector<std::size_t> neighborInsideScvIndices;
            for (const auto& scvf : neighborScvfs(localView))
                neighborInsideScvIndices.push_back(localView.insideScv(scvf).dofIndex());

            if (neighborScvIndices.size() < GridView::dimension)
                DUNE_THROW(Dune::InvalidStateException,
                           "Expected at least " << int(GridView::dimension)
                                                << " neighbors, found: " << neighborScvIndices.size());

            if (neighborInsideScvIndices.size() != neighborScvIndices.size())
                DUNE_THROW(Dune::InvalidStateException, "Outside scv/scvf size mismatch");

            std::ranges::sort(neighborScvIndices);
            std::ranges::sort(neighborInsideScvIndices);
            if (!std::ranges::equal(neighborInsideScvIndices, neighborScvIndices))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected neighbor scv indices");

            const auto numNeighbors = neighborScvIndices.size();
            std::ranges::unique(neighborScvIndices);
            if (neighborScvIndices.size() != numNeighbors)
                DUNE_THROW(Dune::InvalidStateException, "Duplicate neighbor scvs");
        }

        count = 0;
        for (const auto& scvf : scvfs(localView))
        {
            const auto c = scvf.center();
            auto d = c - element.geometry().center();
            d /= d.two_norm();

            if (!equal(scvf.center(), localView.geometry(scvf).center(), tol))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected face center");
            if (!Dune::FloatCmp::eq(scvf.area(), localView.geometry(scvf).volume(), tol))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected face area " << localView.geometry(scvf).volume() << " vs " << scvf.area());
            if (!Dune::FloatCmp::eq(d*scvf.unitOuterNormal(), 1.0, tol))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf normal");
            if (localView.insideScv(scvf).dofIndex() != scvDofIdx)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected inside scv");
            if (localView.indexInElement(scvf) != count)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf index in element");
            if (!localView.onBoundary(scvf))
            {

                d = scvf.center()
                    - localView.insideScv(scvf).center();
                d /= d.two_norm();
                if (!Dune::FloatCmp::eq(d*scvf.unitOuterNormal(), 1.0, tol))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf normal");
                if (withStencil)
                {
                    d = localView.outsideScv(scvf).center()
                        - localView.insideScv(scvf).center();
                    d /= d.two_norm();
                    if (!Dune::FloatCmp::eq(d*scvf.unitOuterNormal(), 1.0, tol))
                        DUNE_THROW(Dune::InvalidStateException, "Unexpected outside scv center");
                    if (!equal(scvf.center(), localView.flipScvf(scvf).center(), tol))
                        DUNE_THROW(Dune::InvalidStateException, "Unexpected Flip Scvf center");
                    if (localView.insideScv(localView.flipScvf(scvf)).dofIndex()
                        != localView.outsideScv(scvf).dofIndex())
                        DUNE_THROW(Dune::InvalidStateException, "insiceScv(flipScvf) != outsideScv(scvf)");
                }
            }
            else if (!isOnDomainBoundary(gridGeometry, c))
                DUNE_THROW(Dune::InvalidStateException, "Unexpected boundary face position: " << c);

            count++;
        }

        if (count != element.subEntities(1))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvfs: " << count);
    };

    for (const auto& element : elements(gridView))
    {
        auto fvGeometry = localView(gridGeometry);

        fvGeometry.bindElement(element);
        testLocalView(element, fvGeometry, false);
        testLocalView(element, localView(gridGeometry).bindElement(element), false);

        fvGeometry.bind(element);
        testLocalView(element, fvGeometry, true);
        testLocalView(element, localView(gridGeometry).bind(element), true);
    }
}

template<typename GV>
struct DynamicallySizedTraits : public Dumux::DefaultMapperTraits<GV>
{
    static constexpr auto maxNumScvfsPerElement = Dumux::Size::dynamic;
    static constexpr auto maxNumBranchesPerScvf = Dumux::Size::dynamic;
};

void test(const Dune::YaspGrid<2>& grid)
{
    using GV = std::decay_t<decltype(grid.leafGridView())>;
    using DynamicTraits = DynamicallySizedTraits<GV>;

    const std::size_t numScvf = grid.leafGridView().size(0)*4;
    testCCTpfaGridGeometry<GV, true>(grid.leafGridView(), numScvf);
    testCCTpfaGridGeometry<GV, false>(grid.leafGridView(), numScvf);
    testCCTpfaGridGeometry<GV, true, DynamicTraits>(grid.leafGridView(), numScvf);
    testCCTpfaGridGeometry<GV, false, DynamicTraits>(grid.leafGridView(), numScvf);
}

void test(const Dune::YaspGrid<3>& grid)
{
    const std::size_t numScvf = grid.leafGridView().size(0)*6;
    testCCTpfaGridGeometry(grid.leafGridView(), numScvf);
}

int main (int argc, char *argv[])
{
    test(Dune::YaspGrid<2>{{1.0, 1.0}, {10, 10}});
    test(Dune::YaspGrid<3>{{1.0, 1.0, 1.0}, {10, 10, 10}});

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
