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

#include <dumux/common/intersectionmapper.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/staggered/fvelementgeometry.hh>
#include <dumux/discretization/staggered/fvgridgeometry.hh>
#include <dumux/discretization/staggered/freeflow/subcontrolvolumeface.hh>
#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>
#include <dumux/discretization/staggered/freeflow/connectivitymap.hh>

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

//! the fv grid geometry traits for this test
template<class GridView>
struct TestFVGGTraits : public DefaultMapperTraits<GridView>
{
    using SubControlVolume = CCSubControlVolume<GridView>;
    using SubControlVolumeFace = FreeFlowStaggeredSubControlVolumeFace<GridView>;
    using IntersectionMapper = ConformingGridIntersectionMapper<GridView>;
    using GeometryHelper = FreeFlowStaggeredGeometryHelper<GridView>;

    struct DofTypeIndices
    {
        using CellCenterIdx = Dune::index_constant<0>;
        using FaceIdx = Dune::index_constant<1>;
    };

    template<class FVGridGeometry>
    using ConnectivityMap = StaggeredFreeFlowConnectivityMap<FVGridGeometry>;

    template<class FVGridGeometry, bool enableCache>
    using LocalView = StaggeredFVElementGeometry<FVGridGeometry, enableCache>;
};

} // end namespace Dumux
#endif

int main (int argc, char *argv[]) try
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    using Grid = Dune::YaspGrid<2>;

    constexpr int dim = Grid::dimension;

    using FVGridGeometry = StaggeredFVGridGeometry<typename Grid::LeafGridView, /*enable caching=*/ true,
                                                   TestFVGGTraits<typename Grid::LeafGridView> >;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, dim> els{{2, 4}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);
    auto leafGridView = grid->leafGridView();

    FVGridGeometry fvGridGeometry(leafGridView);
    fvGridGeometry.update();

    std::cout << "Abbreviatons:\n"
              << "pos - postition of face center\n"
              << "fIdx - face index\n"
              << "dofIdx (self/oppo.) - dofIdx on face (self/opposite)\n"
              << "norm in/out - dofIdx on face normal to own face (within own element / in adjacent element)\n"
              << "par - dofIdx on face parallel to own face\n"
              << "normDist - distance bewteen the dofs on the faces normal the own face\n"
              << "parDist - distance bewteen the dof on the parallel face and the one on the own face\n"
              << "norm in/out - dofIdx on face normal to own face (within own element / in adjacent element)" << std::endl;

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = fvGridGeometry.elementMapper().index(element);
        std::cout << std::endl << "Checking fvGeometry of element " << eIdx << std::endl;
        auto fvGeometry = localView(fvGridGeometry);
        fvGeometry.bind(element);

        auto range = scvs(fvGeometry);
        Detail::NoopFunctor<SubControlVolume> op;
        if(0 != testForwardIterator(range.begin(), range.end(), op))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        for (auto&& scv : scvs(fvGeometry))
        {
            std::cout << "-- scv " << scv.dofIndex() << " center at: " << scv.center() << std::endl;
        }

        auto range2 = scvfs(fvGeometry);
        Detail::NoopFunctor<SubControlVolumeFace> op2;
        if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");


        for (auto&& scvf : scvfs(fvGeometry))
        {
            std::cout <<  std::fixed << std::left << std::setprecision(2)
            << "pos "<< scvf.ipGlobal()
            << "; fIdx "  << std::setw(3)  << scvf.index()
            << "; dofIdx (self/oppo.) " << std::setw(3) << scvf.dofIndex() << "/" << std::setw(3) <<scvf.dofIndexOpposingFace()
            << "; dist (self/oppo.) " << std::setw(3) << scvf.selfToOppositeDistance()
            << ", norm1 in/out " << std::setw(3) << scvf.pairData(0).normalPair.first << "/" << std::setw(3) << scvf.pairData(0).normalPair.second
            << ", norm2 in/out " << std::setw(3) << scvf.pairData(1).normalPair.first << "/" << std::setw(3) << scvf.pairData(1).normalPair.second
            << ", par1 " << std::setw(3)  << std::setw(3) << scvf.pairData(0).outerParallelFaceDofIdx
            << ", par2 " << std::setw(3)  << std::setw(3) << scvf.pairData(1).outerParallelFaceDofIdx
            << ", normDist1 " << std::setw(3) << scvf.pairData(0).normalDistance
            << ", normDist2 " << std::setw(3) << scvf.pairData(1).normalDistance
            << ", parDist1 " << std::setw(3) << scvf.pairData(0).parallelDistance
            << ", parDist2 " << std::setw(3) << scvf.pairData(1).parallelDistance;
            if (scvf.boundary()) std::cout << " (on boundary)";
            std::cout << std::endl;
        }
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception &e) {

    std::cout << e << std::endl;
    return 1;
}
