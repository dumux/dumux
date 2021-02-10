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
 * \ingroup RANSModel
 * \copydoc Dumux::RANSGridVariables
 */
#ifndef DUMUX_RANS_GRIDVARIABLES_HH
#define DUMUX_RANS_GRIDVARIABLES_HH

#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include "model.hh"

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Reynolds-Averaged Navier-Stokes grid variables base class.
 *
 * This implements some base functionality for RANS models.
 * Especially vectors containing all wall-relevant properties, which are accessed
 * by the volumevariables.
 */
template<class NavierStokesGridVariables>
class RANSGridVariables : public NavierStokesGridVariables
{
    using ParentType = NavierStokesGridVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

    static constexpr auto dim = GridView::dimension;
    static constexpr int numCorners = SubControlVolumeFace::numCornersPerFace;
    using DimVector = GlobalPosition;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    struct WallElementInformation
    {
        // store the element indicies for all elements with an intersection on the wall
        unsigned int wallElementIdx;
        // for each wall element, store the faces normal axis
        unsigned int wallFaceNormalAxis;
        // for each wall element, store the location of the face center and each corner.
        GlobalPosition wallFaceCenter;
        std::array<GlobalPosition, numCorners> wallFaceCorners;
    };

public:
    /*!
     * \brief The constructor
     * \param ?
     */
    template<class Problem>
    RANSGridVariables(std::shared_ptr<Problem> problem,
                      std::shared_ptr<FVGridGeometry> fvGridGeometry)
    {
        if (!wallInitialized_)
            initializeWallDistance(fvGridGeometry);
        else if (!neighborsInitialized_)
            storeNeighborIndices();
    }

    /*!
     * \brief Update the (solution independent) relations to the walls
     *
     * This function finds all  wall intersection, and calculates the
     * distance to the nearest wall intersection for each element..
     */
    template<class GridGeometry>
    void initializeWallDistance(const GridGeometry& gridGeometry)
    {
        std::cout << "Determining wall distances. ";
        const auto numElements =  this->fvGridGeometry().gridView().size(0);

        // update size and initial values of the global vectors
        isWallBound_.resize(numElements, false);
        wallDistance_.resize(numElements, std::numeric_limits<Scalar>::max());

        // collect all wall elements and fill wall related vectors
        fillWallDistances_(getWallElements_(gridGeometry));
        wallInitialized_ = true;
    }

    /*!
     * \brief Update the static neighbor indicies for each element
     *
     * This function collects stores the indicies for all
     * neighbor elements for each element
     */
    void storeNeighborIndices()
    {
        std::cout << "Determining cell neighbor indicies. ";
        const auto numElements =  this->fvGridGeometry().gridView().size(0);

        // Store the cell center and neighbor Idx
        neighborIdx_.resize(numElements);

        // search for neighbor Idxs
        fillNeighborIndices_();
        neighborsInitialized_ = true;
    }

    Scalar wallDistance(const int elementIdx) const
    { return wallDistance_[elementIdx]; }

    bool isWallBound(const int elementIdx)
    { return isWallBound_[elementIdx]; }

    unsigned int neighborIndex(const int elementIdx, const int dimIdx, const int sideIdx) const
    { return neighborIdx_[elementIdx][dimIdx][sideIdx];}

private:

    /*!
     * \brief Collects a vector of all elements adjacent to a wall
     *
     * \param scvf The sub control volume face.
     */
    std::vector<WallElementInformation> getWallElements_(const GridGeometry& gridGeometry) const
    {
        std::vector<WallElementInformation> wallElements;

        const auto gridView = gridGeometry.gridView();
        auto fvGeometry = localView(gridGeometry);

        for (const auto& element : elements(gridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                // only search for walls at a global boundary
                if (!scv.boundary())
                    continue;

                for (const auto& scvf : scvfs(fvGeometry, scv))
                {
                    if (!scvf.isBoundary() || scvf.isLateral())
                        continue;

                    if (isOnWall_(scvf))
                    {
                        isWallBound_[scv.elementIndex()] = true;

                        // store the location of the wall adjacent face's center and all corners
                        WallElementInformation wallElementInformation;
                        wallElementInformation.wallFaceCenter = scvf.center();
                        wallElementInformation.wallFaceCorners = fvGeometry.faceCornerPositions(scvf, element);

                        // Store the wall element and face index and face's normal direction (used only with isFlatWallBounded on)
                        wallElementInformation.wallElementIdx = scv.elementIndex();
                        wallElementInformation.wallElementFaceIdx = scvf.index();
                        wallElementInformation.wallFaceNormalAxis = scvf.directionIndex();

                        // add the wall element to the vector of wall elements
                        wallElements.push_back(wallElementInformation);
                    }
                }
            }
        }
        // output the number of wall adjacent faces. Check that this is non-zero.
        std::cout << "NumWallIntersections=" << wallElements.size() << std::endl;
        if (wallElements.size() == 0)
            DUNE_THROW(Dune::InvalidStateException,
                       "No wall intersections have been found. Make sure that the isOnWall(globalPos) is working properly.");

        return wallElements;
    }

    void fillWallDistances_(const std::vector<WallElementInformation>& wallElements)
    {
        const auto gridView = gridGeometry.gridView();
        // search for shortest distance to the wall for each element
        for (const auto& element : elements(gridView))
        {
            // Store the cell center position for each element
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            GlobalPosition cellCenter = element.geometry().center();

            for (unsigned int i = 0; i < wallElements.size(); ++i)
            {
                // Find the minimum distance from the cell center to the wall face (center and corners)
                std::array<Scalar,numCorners+1> cellToWallDistances;
                for (unsigned int j = 0; j < numCorners; j++)
                    cellToWallDistances[j] = (cellCenter - wallElements[i].wallFaceCorners[j]).two_norm();
                cellToWallDistances[numCorners] = (cellCenter - wallElements[i].wallFaceCenter).two_norm();
                Scalar distanceToWall = *std::min_element(cellToWallDistances.begin(), cellToWallDistances.end());

                if (distanceToWall < wallDistance_[elementIdx])
                    wallDistance_[elementIdx] = distanceToWall;
            }
        }
    }

    void fillNeighborIndices_() const
    {
        const auto gridView = gridGeometry.gridView();
        for (const auto& element : elements(gridView))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                neighborIdx_[elementIdx][dimIdx][0] = elementIdx;
                neighborIdx_[elementIdx][dimIdx][1] = elementIdx;
            }

            for (const auto& intersection : intersections(gridView, element))
            {
                if (intersection.boundary())
                    continue;

                unsigned int neighborIdx = this->gridGeometry().elementMapper().index(intersection.outside());
                for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(cellCenter(elementIdx)[dimIdx] - cellCenter(neighborIdx)[dimIdx]) > 1e-8)
                    {
                        if (cellCenter(elementIdx)[dimIdx] > cellCenter(neighborIdx)[dimIdx])
                            neighborIdx_[elementIdx][dimIdx][0] = neighborIdx;

                        if (cellCenter(elementIdx)[dimIdx] < cellCenter(neighborIdx)[dimIdx])
                            neighborIdx_[elementIdx][dimIdx][1] = neighborIdx;
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns whether a given sub control volume face is on a wall
     * \param scvf The sub control volume face.
     */
    bool isOnWall_(const SubControlVolumeFace& scvf) const
    {
        // return ( ( scvf.isBoundary()) && boundaryTypes(scvf).isWall());

        // {
        // boundaryTypes boundaryTypes; // ??
        // return boundaryTypes().isWall(scvf); // How to get boundary types
        // }

        if (scvf.isBoundary())
            return true; // TODO: SOMEHOW EVALUATE THE BOUNDARY TYPES?
        else
            return false;
    }

    std::vector<Scalar> wallDistance_;
    std::vector<bool> isWallBound_;
    std::vector<std::array<std::array<unsigned int, 2>, dim>> neighborIdx_;

    bool neighborsInitialized_= false;
    bool wallInitialized_ = false;
};

} // end namespace Dumux

#endif
