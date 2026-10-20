// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the two-ring connectivity map of the cell-centered tpfa grid geometry
 */
#include <config.h>

#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

namespace Dumux {

// expected map computed from the grid intersections: element i maps to the elements j
// and their scvfs whose stencil (face neighbors of j and of the outside element) contains i
template<class GridGeometry>
std::vector<std::map<std::size_t, std::set<std::size_t>>> expectedMap(const GridGeometry& gridGeometry)
{
    const auto& gridView = gridGeometry.gridView();
    const auto& mapper = gridGeometry.elementMapper();

    std::vector<std::set<std::size_t>> neighbors(gridView.size(0));
    for (const auto& element : elements(gridView))
        for (const auto& intersection : intersections(gridView, element))
            if (intersection.neighbor())
                neighbors[mapper.index(element)].insert(mapper.index(intersection.outside()));

    std::vector<std::map<std::size_t, std::set<std::size_t>>> expected(gridView.size(0));
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridView))
    {
        const auto j = mapper.index(element);
        fvGeometry.bindElement(element);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            auto stencil = neighbors[j];
            if (!scvf.boundary())
                stencil.insert(neighbors[scvf.outsideScvIdx()].begin(), neighbors[scvf.outsideScvIdx()].end());

            for (const auto i : stencil)
                if (i != j)
                    expected[i][j].insert(scvf.index());
        }
    }

    return expected;
}

template<class GridView>
void checkTwoRingMap(const GridView& gridView, const std::string& name)
{
    using Traits = CCTpfaTwoRingGridGeometryTraits<GridView>;
    using GridGeometry = CCTpfaFVGridGeometry<GridView, true, Traits>;
    GridGeometry gridGeometry(gridView);

    const auto expected = expectedMap(gridGeometry);
    const auto& map = gridGeometry.connectivityMap();
    for (const auto& element : elements(gridView))
    {
        const auto i = gridGeometry.elementMapper().index(element);
        std::map<std::size_t, std::set<std::size_t>> actual;
        for (const auto& dataJ : map[i])
        {
            if (actual.count(dataJ.globalJ))
                DUNE_THROW(Dune::Exception, name << ": element " << dataJ.globalJ << " listed twice for element " << i);
            actual[dataJ.globalJ].insert(dataJ.scvfsJ.begin(), dataJ.scvfsJ.end());
            if (actual[dataJ.globalJ].size() != dataJ.scvfsJ.size())
                DUNE_THROW(Dune::Exception, name << ": duplicate scvf for element " << i);
        }

        if (actual != expected[i])
            DUNE_THROW(Dune::Exception, name << ": connectivity of element " << i << " differs from the expected stencil");
    }

    std::cout << name << ": two-ring connectivity map matches the stencil derived from the intersections" << std::endl;
}

} // end namespace Dumux

int main (int argc, char *argv[])
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using GridView = typename Grid::LeafGridView;

    // the default traits keep the stencil size of a face-neighbor stencil
    static_assert(CCTpfaFVGridGeometry<GridView, true>::maxElementStencilSize == 4*2 + 1);
    static_assert(CCTpfaFVGridGeometry<GridView, false>::maxElementStencilSize == 4*2 + 1);
    static_assert(CCTpfaFVGridGeometry<GridView, true, CCTpfaTwoRingGridGeometryTraits<GridView>>::maxElementStencilSize == 8*8 + 1);

    const Dune::FieldVector<double, 2> lower(0.0), upper(1.0);
    const auto yasp = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, std::array<unsigned int, 2>{{6, 6}});
    checkTwoRingMap(yasp->leafGridView(), "YaspGrid 6x6");

    // on a structured quadrilateral grid, an element at least two layers away from the boundary depends on
    // 4 face neighbors (all 4 of their scvfs), 4 elements two steps along the axes (1 scvf each)
    // and 4 diagonal elements (2 scvfs each)
    {
        using GridGeometry = CCTpfaFVGridGeometry<GridView, true, CCTpfaTwoRingGridGeometryTraits<GridView>>;
        GridGeometry gridGeometry(yasp->leafGridView());
        for (const auto& element : elements(yasp->leafGridView()))
        {
            const auto center = element.geometry().center();
            if (center[0] < 2.0/6.0 || center[0] > 4.0/6.0 || center[1] < 2.0/6.0 || center[1] > 4.0/6.0)
                continue;

            const auto& dataJForI = gridGeometry.connectivityMap()[gridGeometry.elementMapper().index(element)];
            std::size_t numScvfs = 0;
            for (const auto& dataJ : dataJForI)
                numScvfs += dataJ.scvfsJ.size();
            if (dataJForI.size() != 12 || numScvfs != 28)
                DUNE_THROW(Dune::Exception, "Interior element depends on " << dataJForI.size() << " elements and "
                                             << numScvfs << " scvfs, expected 12 and 28");
        }
    }

#if HAVE_DUNE_UGGRID
    using SimplexGrid = Dune::UGGrid<2>;
    const auto simplex = Dune::StructuredGridFactory<SimplexGrid>::createSimplexGrid(lower, upper, std::array<unsigned int, 2>{{5, 4}});
    checkTwoRingMap(simplex->leafGridView(), "UGGrid simplices 5x4");
#endif

    return 0;
}
