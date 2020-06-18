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
 * \brief Test for finite volume element geometry, sub control volume, and sub
          control volume faces
 */
#include <config.h>

#include <iostream>
#include <utility>
#include <iomanip>

#include <dune/common/test/iteratortest.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/intersectionmapper.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/facecentered/staggered/fvgridgeometry.hh>

#ifndef DOXYGEN
namespace Dumux {
namespace Detail {
template<class T>
class NoopFunctor {
public:
  NoopFunctor() {}
  void operator()(const T& t){}
};
} // end namespace Detail

} // end namespace Dumux
#endif

int main (int argc, char *argv[]) try
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    using Grid = Dune::YaspGrid<2>;

    constexpr int dim = Grid::dimension;
    constexpr int dimworld = Grid::dimensionworld;

    using GlobalPosition = typename Dune::FieldVector<Grid::ctype, dimworld>;
    using GridGeometry = FaceCenteredStaggeredFVGridGeometry<typename Grid::LeafGridView, /*enable caching=*/ true>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    // make a grid
    GlobalPosition lower = getParam<GlobalPosition>("Grid.LowerLeft", GlobalPosition(0.0));
    GlobalPosition upper = getParam<GlobalPosition>("Grid.UpperRight", GlobalPosition(1.0));
    const auto cells = getParam<std::array<unsigned int, dim>>("Grid.Cells", std::array<unsigned int, dim>{1,1});
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, cells);
    auto leafGridView = grid->leafGridView();

    GridGeometry gridGeometry(leafGridView);
    gridGeometry.update();

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        std::cout << std::endl << "Checking fvGeometry of element " << eIdx << std::endl;
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bind(element);

        auto range = scvs(fvGeometry);
        Detail::NoopFunctor<SubControlVolume> op;
        if(0 != testForwardIterator(range.begin(), range.end(), op))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        for (auto&& scv : scvs(fvGeometry))
        {
            std::cout << "-- scv " << scv.index() << " center at: " << scv.center() << ", dofPosition " << scv.dofPosition() << " , volume: " << scv.volume()  << std::endl;
        }

        auto range2 = scvfs(fvGeometry);
        Detail::NoopFunctor<SubControlVolumeFace> op2;
        if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        std::size_t boundaryCount = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            std::cout << "-- scvf " << scvf.index() << " ip at: " << scvf.ipGlobal() << " normal: " << scvf.unitOuterNormal();

            if (scvf.isFrontal())
                std::cout << ", frontal ";
            else
                std::cout << ", lateral ";

            if (scvf.boundary())
            {
                ++boundaryCount;
                std::cout << " (on boundary).";
            }

            std::cout << std::endl;
        }

        if ((boundaryCount > 0) != fvGeometry.hasBoundaryScvf())
            DUNE_THROW(Dune::InvalidStateException, "fvGeometry.hasBoundaryScvf() reports " << fvGeometry.hasBoundaryScvf()
                            << " but the number of boundary scvfs is " << boundaryCount);
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception& e) {

    std::cout << e << std::endl;
    return 1;
}
