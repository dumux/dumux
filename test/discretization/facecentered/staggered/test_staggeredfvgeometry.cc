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
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/common/intersectionmapper.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/facecentered/staggered/fvgridgeometry.hh>

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

int main (int argc, char *argv[]) try
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    using Grid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2> >;
    GridManager<Grid> gridManager;
    gridManager.init();

    constexpr int dim = Grid::dimension;
    constexpr int dimworld = Grid::dimensionworld;

    using GlobalPosition = typename Dune::FieldVector<Grid::ctype, dimworld>;
    using GridGeometry = FaceCenteredStaggeredFVGridGeometry<typename Grid::LeafGridView,
                                                             /*enable caching=*/ true,
                                                             /*upwindSchemeOrder*/ 2>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    constexpr bool useHigherOrder = GridGeometry::useHigherOrder;

    // make a grid
    const auto& leafGridView = gridManager.grid().leafGridView();
    GridGeometry gridGeometry(leafGridView);
    gridGeometry.update();

    if constexpr (dim == 2 && dimworld == 2)
        Dumux::drawGridGeometry(gridGeometry, "FaceCenteredIndicies", 4000, true, true);

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        if (eIdx == 0 || eIdx == 12)
        {
            std::cout << std::endl << "Checking fvGeometry of element " << eIdx << std::endl;
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bind(element);

            auto range = scvs(fvGeometry);
            Detail::NoopFunctor<SubControlVolume> op;
            if(0 != testForwardIterator(range.begin(), range.end(), op))
                DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

            auto range2 = scvfs(fvGeometry);
            Detail::NoopFunctor<SubControlVolumeFace> op2;
            if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
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
                    std::cout << " scvfIdxWithCommonEntity: " << scvf.scvfIdxWithCommonEntity()
                              << ", inside and outside SCVs:" << scvf.insideScvIdx() << ", " << scvf.outsideScvIdx();

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

    // Check indicies for all dofs called in the stencil. Center element, scv local index 1 (global index 49)
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        if (eIdx == 12)
        {
            auto fvGeometry = localView(gridGeometry);
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
                        std::cout << "-- -- These dofs are located " << fvGeometry.selfToOppositeDistance(scvf) << " units away from each other \n";
                        if constexpr (useHigherOrder)
                        {
                            // Find the dofs located at the forward and backward axial scvs
                            if (fvGeometry.hasForwardNeighbor(scvf))
                            {
                                std::cout << "-- -- [Higher Order] The \"forward\" velocity is located at scv "
                                          << fvGeometry.forwardScvIdx(scvf) << "\n";
                                std::cout << "-- -- [Higher Order] The \"forward\" dof is located " << fvGeometry.selfToForwardDistance(scvf)
                                          << " units from the self dof \n";
                            }

                            if (fvGeometry.hasBackwardNeighbor(scvf))
                            {
                                std::cout << "-- -- [Higher Order] The \"backward\" velocity is located at scv "
                                          << fvGeometry.backwardScvIdx(scvf) << "\n";
                                std::cout << "-- -- [Higher Order] The \"backward\" dof is located " << fvGeometry.oppositeToBackwardDistance(scvf)
                                          << " units from the opposite dof \n";
                            }
                        }
                    }

                    int latCount = 0;
                    std::cout << "-- When evaluating the advective flux across the LATERAL face(s) the following scvfs and dofs are involved: \n";
                    for (auto&& scvf : scvfs(fvGeometry, scv))
                    {
                        if(!scvf.isLateral())
                            continue;

                        // Write out the indicies for the transported velocity calculation
                        std::string firstOrSecond = (latCount < 1) ? "first " : "second ";
                        std::cout << "-- -- The " << firstOrSecond << "lateral scvf has the index " << scvf.index() << "\n";
                        std::cout << "-- -- - The \"inner parallel\" transported velocity is located at scv " << scvf.insideScvIdx() << "\n";
                        std::cout << "-- -- - The \"inner parallel\" scv is " << fvGeometry.insideScvLateralLength(scvf)
                                  << " units in length "<< "\n";
                        if (fvGeometry.hasParallelNeighbor(scvf))
                        {
                            std::cout << "-- -- - The \"1st outer parallel\" transported velocity is located at scv " << scvf.outsideScvIdx() << "\n";
                            std::cout << "-- -- - The \"1st outer parallel\" dof is located " << fvGeometry.selfToParallelDistance(scvf)
                                      << " units away from the self dof \n";
                            std::cout << "-- -- - The \"1st outer parallel\" Scv is " << fvGeometry.outsideScvLateralLength(scvf)
                                      << " units in length "<< "\n";
                        }
                        if constexpr (useHigherOrder)
                        {
                            // Find the dofs related to the second parallel scvs
                            if (fvGeometry.hasSecondParallelNeighbor(scvf))
                            {
                                std::cout << "-- -- - [Higher Order] The \"2nd outer parallel\" transported velocity is located at scv "
                                          << fvGeometry.secondParallelScvIdx(scvf) << "\n";
                                std::cout << "-- -- - [Higher Order] The \"2nd outer parallel\" dof is located "
                                          << fvGeometry.paralleltoSecondParallelDistance(scvf)
                                          << " units away from the first parallel dof \n";
                            }
                        }
                        // Write out the indicies for the transporting velocity calculation
                        const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                        std::cout << "-- -- - The orthagonal scvf has the index " << orthogonalScvf.index() << "\n";
                        std::cout << "-- -- - The \"inner transporting velocity is located at scv " << orthogonalScvf.insideScvIdx() << "\n";
                        std::cout << "-- -- - The \"outer transporting velocity is located at scv " << orthogonalScvf.outsideScvIdx() << "\n";
                        latCount++;
                    }
                }
            }
        }
        else if (eIdx == 22)
        {
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bind(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto scvIdx = scv.index();
                if (scvIdx == 89)
                {
                    for (auto&& scvf : scvfs(fvGeometry, scv))
                    {
                        if (scvf.index() == 284)
                        {
                            std::cout << "\n \n Checking the outer parallel lateral face \n";
                            std::cout << "scvf index is : " << scvf.index()
                                      << " and the outer parallel lateral face index is: " << fvGeometry.outerParallelLateralScvf(scvf).index() << "\n";
                        }
                    }
                }
            }
        }
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception& e) {

    std::cout << e << std::endl;
    return 1;
}
