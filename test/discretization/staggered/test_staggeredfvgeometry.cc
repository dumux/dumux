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
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/staggered/fvelementgeometry.hh>
#include <dumux/discretization/staggered/fvgridgeometry.hh>
#include <dumux/discretization/staggered/subcontrolvolumeface.hh>

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
    using SubControlVolumeFace = StaggeredSubControlVolumeFace<GridView>;
    using IntersectionMapper = ConformingGridIntersectionMapper<GridView>;
    using GeometryHelper = BaseStaggeredGeometryHelper<GridView>;

    struct DofTypeIndices
    {
        using CellCenterIdx = Dune::index_constant<0>;
        using FaceIdx = Dune::index_constant<1>;
    };

    //! Dummy connectivity map, required by GridGeometry
    template<class GridGeometry>
    struct MockConnectivityMap
    {
        void update(const GridGeometry& gridGeometry) {}
        void setStencilOrder(const int order) {}
    };

    template<class GridGeometry>
    using ConnectivityMap = MockConnectivityMap<GridGeometry>;

    template<class GridGeometry, bool enableCache>
    using LocalView = StaggeredFVElementGeometry<GridGeometry, enableCache>;

    struct PublicTraits
    {
        using CellSubControlVolume = SubControlVolume;
        using CellSubControlVolumeFace = SubControlVolumeFace;
        using FaceSubControlVolume = SubControlVolume;
        using FaceLateralSubControlVolumeFace = SubControlVolumeFace;
        using FaceFrontalSubControlVolumeFace = SubControlVolumeFace;
    };
};

} // end namespace Dumux
#endif

int main (int argc, char *argv[])
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
    using GridGeometry = StaggeredFVGridGeometry<typename Grid::LeafGridView, /*enable caching=*/ true,
                                                   TestFVGGTraits<typename Grid::LeafGridView> >;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, dim> els{{2, 4}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);
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
            std::cout << "-- scv " << scv.localDofIndex() << " center at: " << scv.center() << " , volume: " << scv.volume()  << std::endl;
        }

        auto range2 = scvfs(fvGeometry);
        Detail::NoopFunctor<SubControlVolumeFace> op2;
        if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        std::size_t boundaryCount = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            std::cout << "-- scvf " << scvf.index() << " ip at: " << scvf.ipGlobal() << " normal: " << scvf.unitOuterNormal();
            if (scvf.boundary())
            {
                ++boundaryCount;
                std::cout << " (on boundary).";
            }
            std::cout << std::endl;
        }

        if ((boundaryCount>0) != fvGeometry.hasBoundaryScvf())
            DUNE_THROW(Dune::InvalidStateException, "fvGeometry.hasBoundaryScvf() reports " << fvGeometry.hasBoundaryScvf()
                            << " but the number of boundary scvfs is " << boundaryCount);
    }
}
