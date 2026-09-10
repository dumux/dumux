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

#include <cmath>
#include <iostream>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/test/iteratortest.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

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

    using Grid = Dune::YaspGrid<2>;

    constexpr int dim = Grid::dimension;

    using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, ENABLE_CACHING>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, dim> els{{2, 2}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);

    // obtain leaf and make GridGeometry
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
            std::cout << "-- scv " << scv.dofIndex() << " center at: " << scv.center() << std::endl;
        }

        auto range2 = scvfs(fvGeometry);
        Detail::NoopFunctor<SubControlVolumeFace> op2;
        if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        std::size_t boundaryCount = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            std::cout << "-- scvf " << scvf.index() << " ip at: " << scvf.ipGlobal();
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

        // there is one boundary face per boundary intersection, and it carries the geometry of it
        std::size_t faceCount = 0;
        for (const auto& face : boundaryFaces(fvGeometry))
        {
            ++faceCount;
            bool matched = false;
            for (const auto& intersection : intersections(leafGridView, element))
                if (intersection.indexInInside() == face.intersectionIndex())
                {
                    matched = true;
                    if (!intersection.boundary())
                        DUNE_THROW(Dune::InvalidStateException, "Boundary face on an interior intersection");
                    if ((intersection.geometry().center() - face.center()).two_norm() > 1e-12)
                        DUNE_THROW(Dune::InvalidStateException, "Boundary face center does not match its intersection");
                    if ((intersection.centerUnitOuterNormal() - face.unitOuterNormal()).two_norm() > 1e-12)
                        DUNE_THROW(Dune::InvalidStateException, "Boundary face normal does not match its intersection");
                    if (std::abs(intersection.geometry().volume() - face.area()) > 1e-12*face.area())
                        DUNE_THROW(Dune::InvalidStateException, "Boundary face area does not match its intersection");
                }
            if (!matched)
                DUNE_THROW(Dune::InvalidStateException, "Boundary face has no intersection with its index");

            // a cell-centered scheme has exactly one sub-control volume face per intersection
            std::size_t subFaceCount = 0;
            for (const auto& scvf : scvfs(fvGeometry, face))
            {
                ++subFaceCount;
                if (!scvf.boundary())
                    DUNE_THROW(Dune::InvalidStateException, "Sub-control volume face on a boundary face is not a boundary face");
                if ((scvf.center() - face.center()).two_norm() > 1e-12)
                    DUNE_THROW(Dune::InvalidStateException, "Sub-control volume face does not lie on its boundary face");
            }
            if (subFaceCount != 1)
                DUNE_THROW(Dune::InvalidStateException, "Expected one sub-control volume face per boundary face, got " << subFaceCount);
        }

        if (faceCount != boundaryCount)
            DUNE_THROW(Dune::InvalidStateException, "Found " << faceCount << " boundary faces but " << boundaryCount << " boundary scvfs");
    }
}
