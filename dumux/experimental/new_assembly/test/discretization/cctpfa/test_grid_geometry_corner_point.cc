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
#include <config.h>
#include <iostream>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/cpgridmanager.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometry.hh>


static constexpr double tolerance = 1e-6;

struct ElementSurface
{
    std::size_t numFaces = 0;
    double area = 0.0;
};

template<typename GridView>
ElementSurface getElementSurface(const GridView& gv, const typename GridView::template Codim<0>::Entity& element)
{
    ElementSurface result;
    for (const auto& is : intersections(gv, element))
    {
        result.numFaces++;
        result.area += is.geometry().volume();
    }
    return result;
}

template<typename LocalView>
ElementSurface getElementSurface(const LocalView& localView)
{
    ElementSurface result;
    for (const auto& scvf : scvfs(localView))
    {
        result.numFaces++;
        result.area += scvf.area();
    }
    return result;
}

template<typename GV>
struct DynamicallySizedTraits : public Dumux::DefaultMapperTraits<GV>
{
    static constexpr auto maxNumScvfsPerElement = Dumux::Size::dynamic;
    static constexpr auto maxNumBranchesPerScvf = Dumux::Size::dynamic;
};

template<bool useCaching, typename GridView>
void testGridGeometry(const GridView& gridView)
{
    Dumux::CCTpfaGridGeometry<GridView, useCaching, DynamicallySizedTraits<GridView>> gridGeometry{gridView};
    for (const auto& element : elements(gridView))
    {
        auto fvGeometry = localView(gridGeometry).bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
            if (Dune::FloatCmp::ne(element.geometry().volume(), scv.volume(), tolerance))
                DUNE_THROW(Dune::InvalidStateException, "Scv volume mismatch");

        const auto elemSurface = getElementSurface(gridView, element);
        const auto controlVolumeSurface = getElementSurface(fvGeometry);
        if (elemSurface.numFaces != controlVolumeSurface.numFaces)
            DUNE_THROW(Dune::InvalidStateException, "Mismatch in number of faces");
        if (Dune::FloatCmp::ne(elemSurface.area, controlVolumeSurface.area, tolerance))
            DUNE_THROW(
                Dune::InvalidStateException,
                "Cell surface area mismatch: " << elemSurface.area << " vs. " << controlVolumeSurface.area
            );
    }
}

int main (int argc, char *argv[])
{
    Dumux::initialize(argc, argv);
    Dune::MPIHelper::instance();
    Dumux::Parameters::init(argc, argv);

    Dumux::CpGridManager gridManager;
    gridManager.init();

    testGridGeometry<false>(gridManager.grid().leafGridView());
    testGridGeometry<true>(gridManager.grid().leafGridView());

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
