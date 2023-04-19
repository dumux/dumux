// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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
#include <iomanip>

#include <dune/common/test/iteratortest.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/intersectionmapper.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/facecentered/staggered/fvgridgeometry.hh>
#include <dumux/io/vtk/intersectionwriter.hh>

#include "drawgrid.hh"

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
    constexpr bool useHigherOrder = true; // do this properly with setProp

    // make a grid
    GlobalPosition lower = getParam<GlobalPosition>("Grid.LowerLeft", GlobalPosition(0.0));
    GlobalPosition upper = getParam<GlobalPosition>("Grid.UpperRight", GlobalPosition(1.0));
    const auto cells = getParam<std::array<unsigned int, dim>>("Grid.Cells", std::array<unsigned int, dim>{1,1});
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, cells);
    auto leafGridView = grid->leafGridView();

    GridGeometry gridGeometry(leafGridView);

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        if (eIdx == 0 || eIdx == 12)
        {
            std::cout << std::endl << "Checking fvGeometry of element " << eIdx << std::endl;
            fvGeometry.bind(element);

            auto range = scvs(fvGeometry);
            Detail::NoopFunctor<SubControlVolume> op;
            if(0 != testForwardIterator(range.begin(), range.end(), op))
                DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

            auto range2 = scvfs(fvGeometry);
            Detail::NoopFunctor<SubControlVolumeFace> op2;
            if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
                DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

            const auto& firstScv = (*scvs(fvGeometry).begin());
            auto range3 = scvfs(fvGeometry, firstScv);
            if(0 != testForwardIterator(range3.begin(), range3.end(), op2))
                DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

            for (auto&& scv : scvs(fvGeometry))
            {
                std::cout << "- scv " << scv.index() << " center at: " << scv.center()
                          << ", dofPosition " << scv.dofPosition() << " , volume: " << scv.volume()
                          << " , normal: "<< int(scv.dofAxis())  << " , in direction: " << int(scv.directionSign());
                std::cout << " contains: \n";

                std::size_t boundaryCount = 0;
                for (auto&& scvf : scvfs(fvGeometry, scv))
                {
                    std::cout << "    scvf (globalIdx) " << scvf.index() << " ip at: " << scvf.ipGlobal() << " normal: " << scvf.unitOuterNormal();
                    if (scvf.isLateral())
                    {
                        const auto& lateralScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                        std::cout << "    lateral scvf (globalIdx ): " << lateralScvf.index() << " ip at: " << lateralScvf.ipGlobal() << " normal: " << lateralScvf.unitOuterNormal();
                    }
                    std::cout << ", inside and outside SCVs:" << scvf.insideScvIdx() << ", " << scvf.outsideScvIdx();

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
        else
            continue;
    }

    // draw the grid with indices
    if constexpr (dim == 2 && dimworld == 2)
        Dumux::drawGridGeometry(gridGeometry, "FaceCenteredIndicies", 4000, true, false);

    // write face dof indices to vtk
    auto intersectionWriter = Dumux::ConformingIntersectionWriter(leafGridView);
    std::vector<std::size_t> dofIndices(gridGeometry.numDofs());
    for (const auto& element : elements(leafGridView))
    {
        fvGeometry.bindElement(element);
        for (auto&& scv : scvs(fvGeometry))
            dofIndices[scv.dofIndex()] = scv.dofIndex();
    }
    intersectionWriter.addField(dofIndices, "dofIdx");
    intersectionWriter.write("staggered");

    // Check indices for all dofs called in the stencil. Center element, scv local index 1 (global index 49)
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        if (eIdx == 12)
        {
            fvGeometry.bind(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto scvIdx = scv.index();
                if (scvIdx == 49)
                {
                    std::cout << "\nEvaluating the dofs associated with the scv with global index "<< scvIdx << " of element " << eIdx << ": \n";
                    for (auto&& scvf : scvfs(fvGeometry, scv))
                    {
                        if(!scvf.isFrontal())
                            continue;

                        std::cout << "-- When evaluating the advective flux across the FRONTAL face the following scvfs and dofs are involved: \n";
                        std::cout << "-- -- The frontal scvf has the index " << scvf.index() << "\n";
                        std::cout << "-- -- The \"self\" velocity is located at scv " << scvf.insideScvIdx() << "\n";
                        std::cout << "-- -- The \"opposite\" velocity is located at scv " << scvf.outsideScvIdx() << "\n";
                        if constexpr (useHigherOrder)
                        {
                            // do higher order stuff
//    PSEUDO                if(scvf.hasForwardNeighbor())
//    CODE                      std::cout << "-- -- The \"forward\" velocity is located at scv " << scvf.forwardScvIdx() << "\n";
//    PSEUDO                if(scvf.hasBackwardNeighbor())
//    CODE                      std::cout << "-- -- The \"backward\" velocity is located at scv " << scvf.backwardScvIdx() << "\n";
                        }
                    }
                    int latCount = 0;
                    std::cout << "-- When evaluating the advective flux across the LATERAL face(s) the following scvfs and dofs are involved: \n";
                    for (auto&& scvf : scvfs(fvGeometry, scv))
                    {
                        if(!scvf.isLateral())
                            continue;

                        // Write out the indices for the transported velocity calculation
                        std::string firstOrSecond = (latCount < 1) ? "first" : "second";
                        std::cout << "-- -- The " << firstOrSecond << "lateral scvf has the index " << scvf.index() << "\n";
                        std::cout << "-- -- - The \"inner\" transported velocity is located at scv " << scvf.insideScvIdx() << "\n";
                        std::cout << "-- -- - The \"outer\" transported velocity is located at scv " << scvf.outsideScvIdx() << "\n";
                        if constexpr (useHigherOrder)
                        {
                            // do higher order stuff for transported velocity
//  PSEUDO                  if (scvf.hasSecondParallelNeighbor())
//  CODE                        std::cout << "-- -- - The \"second outer\" transported velocity is located at scv " << scvf.secondOuterScvIdx() << "\n";

                            // const auto& firstParallelScv = fvGeometry.scv(scvf.outsideScvIdx());
                        }
                        // Write out the indices for the transporting velocity calculation
                        const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                        std::cout << "-- -- - The orthogonal scvf has the index " << orthogonalScvf.index() << "\n";
                        std::cout << "-- -- - The \"inner transporting velocity is located at scv " << orthogonalScvf.insideScvIdx() << "\n";
                        std::cout << "-- -- - The \"outer transporting velocity is located at scv " << orthogonalScvf.outsideScvIdx() << "\n";
                        latCount++;
                    }

                }
            }
        }
    }
}
