// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for finite volume element geometry, sub control volume, and sub
          control volume faces
 */
#include <config.h>

#include <iostream>
#include <utility>

#include <dune/common/test/iteratortest.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>

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

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    using Grid = Dune::YaspGrid<3>;

    constexpr int dim = Grid::dimension;

    using GridGeometry = BoxFVGridGeometry<double, typename Grid::LeafGridView, ENABLE_CACHING>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, dim> els{{1, 1, 1}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);
    auto leafGridView = grid->leafGridView();

    GridGeometry gridGeometry(leafGridView);

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        std::cout << std::endl << "Checking fvGeometry of element " << eIdx << std::endl;
        auto fvGeometry = localView(gridGeometry);

        // bind the local view to the element
        if (fvGeometry.isBound()) DUNE_THROW(Dune::Exception, "Local view should not be bound at this point");
        fvGeometry.bind(element);
        if (!fvGeometry.isBound()) DUNE_THROW(Dune::Exception, "Local view should be bound at this point");

        // make sure the bound element fits
        auto eIdxBound = gridGeometry.elementMapper().index(fvGeometry.element());
        if (eIdx != eIdxBound)
            DUNE_THROW(Dune::Exception, "Bound element index does not match");

        auto range = scvs(fvGeometry);
        Detail::NoopFunctor<SubControlVolume> op;
        if(0 != testForwardIterator(range.begin(), range.end(), op))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        for (auto&& scv : scvs(fvGeometry))
        {
            std::cout << "-- scv " << scv.localDofIndex() << " center at: " << scv.center() << " , volume: " << scv.volume()  << std::endl;

            const auto scvGeo = fvGeometry.geometry(scv);

            if ((scv.center() - scvGeo.center()).two_norm2() > 1e-15)
                DUNE_THROW(Dune::Exception,
                    "scv.center() and fvGeometry.geometry(scv).center() do not match! "
                    << scv.center() << " :: " << scvGeo.center()
                );

            if (std::abs(scv.volume() - scvGeo.volume()) > 1e-7)
                DUNE_THROW(Dune::Exception,
                    "scv.volume() and fvGeometry.geometry(scv).volume() do not match! "
                    << scv.volume() << " :: " << scvGeo.volume()
                );
        }

        auto range2 = scvfs(fvGeometry);
        Detail::NoopFunctor<SubControlVolumeFace> op2;
        if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        std::size_t boundaryCount = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            std::cout << "-- scvf " << scvf.index() << " ip at: " << scvf.ipGlobal() << " normal: " << scvf.unitOuterNormal();

            const auto scvfGeo = fvGeometry.geometry(scvf);

            if ((scvf.ipGlobal() - scvfGeo.center()).two_norm2() > 1e-15)
                DUNE_THROW(Dune::Exception,
                    "scvf.ipGlobal() and fvGeometry.geometry(scvf).center() do not match! "
                    << scvf.ipGlobal() << " :: " << scvfGeo.center()
                );

            if (std::abs(scvf.area() - scvfGeo.volume()) > 1e-7)
                DUNE_THROW(Dune::Exception,
                    "scvf.area() and fvGeometry.geometry(scvf).volume() do not match! "
                    << scvf.area() << " :: " << scvfGeo.volume()
                );

            if (scvf.boundary())
            {
                ++boundaryCount;
                std::cout << " (on boundary).";
            }

            // verify that boundary faces have no neighbor
            if (scvf.boundary() && scvf.numOutsideScvs() != 0)
                DUNE_THROW(Dune::Exception, "Boundary face states that it has a neighbor");

            // verify that non-boundary faces have a single neighbor
            if (!scvf.boundary() && scvf.numOutsideScvs() != 1)
                DUNE_THROW(Dune::Exception, "Expected non-boundary face to have a single neighbor");

            std::cout << std::endl;
        }

        if ((boundaryCount>0) != fvGeometry.hasBoundaryScvf())
            DUNE_THROW(Dune::InvalidStateException, "fvGeometry.hasBoundaryScvf() reports " << fvGeometry.hasBoundaryScvf()
                            << " but the number of boundary scvfs is " << boundaryCount);
    }
}
