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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridVariables
 */
#ifndef DUMUX_RANS_GRID_VARIABLES_HH
#define DUMUX_RANS_GRID_VARIABLES_HH

#include <dumux/discretization/staggered/gridvariables.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class storing data associated to scvs and scvfs
 * \tparam GG the type of the grid geometry
 * \tparam GVV the type of the grid volume variables
 * \tparam GFVC the type of the grid flux variables cache
 * \tparam GFV the type of the grid face variables
 */
template<class NavierStokesGridVariables>
class RANSGridVariables : public NavierStokesGridVariables
{
    using ParentType = NavierStokesGridVariables;
    // using FVGridGeometry = GG;
    // using ThisType = StaggeredGridVariables<GG, GVV, GFVC, GFV>;
    // friend class StaggeredGridVariablesView<ThisType>;

    // static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
    // static constexpr auto faceIdx = FVGridGeometry::faceIdx();

public:

    using ParentType::ParentType;

    /*!
     * \brief Update the static (solution independent) relations to the walls
     *
     * This function determines all element with a wall intersection,
     * the wall distances and the relation to the neighboring elements.
     */
    void updateStaticWallProperties()
    {
        using std::abs;
        std::cout << "Update static wall properties. ";
        calledUpdateStaticWallProperties = true;

        // update size and initial values of the global vectors
        wallElementIdx_.resize(this->fvGridGeometry().elementMapper().size());
        wallDistance_.resize(this->fvGridGeometry().elementMapper().size(), std::numeric_limits<Scalar>::max());
        neighborIdx_.resize(this->fvGridGeometry().elementMapper().size());
        cellCenter_.resize(this->fvGridGeometry().elementMapper().size(), GlobalPosition(0.0));
        velocity_.resize(this->fvGridGeometry().elementMapper().size(), DimVector(0.0));
        velocityMaximum_.resize(this->fvGridGeometry().elementMapper().size(), DimVector(0.0));
        velocityGradients_.resize(this->fvGridGeometry().elementMapper().size(), DimMatrix(0.0));
        stressTensorScalarProduct_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        vorticityTensorScalarProduct_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        flowNormalAxis_.resize(this->fvGridGeometry().elementMapper().size(), 0);
        wallNormalAxis_.resize(this->fvGridGeometry().elementMapper().size(), 1);
        kinematicViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        sandGrainRoughness_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);

        // retrieve all wall intersections and corresponding elements
        std::vector<unsigned int> wallElements;
        std::vector<GlobalPosition> wallPositions;
        std::vector<unsigned int> wallNormalAxisTemp;

        const auto gridView = this->fvGridGeometry().gridView();
        auto fvGeometry = localView(this->fvGridGeometry());

        for (const auto& element : elements(gridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                // only search for walls at a global boundary
                if (!scvf.boundary())
                    continue;

                if (asImp_().isOnWall(scvf))
                {
                    wallElements.push_back(this->fvGridGeometry().elementMapper().index(element));
                    wallPositions.push_back(scvf.center());
                    wallNormalAxisTemp.push_back(scvf.directionIndex());
                }
            }
        }
        std::cout << "NumWallIntersections=" << wallPositions.size() << std::endl;
        if (wallPositions.size() == 0)
            DUNE_THROW(Dune::InvalidStateException,
                       "No wall intersections have been found. Make sure that the isOnWall(globalPos) is working properly.");


        // search for shortest distance to wall for each element
        for (const auto& element : elements(gridView))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            cellCenter_[elementIdx] = element.geometry().center();
            for (unsigned int i = 0; i < wallPositions.size(); ++i)
            {
                static const int problemWallNormalAxis
                    = getParamFromGroup<int>(this->paramGroup(), "RANS.WallNormalAxis", -1);
                int searchAxis = problemWallNormalAxis;

                // search along wall normal axis of the intersection
                if (problemWallNormalAxis < 0 || problemWallNormalAxis >= dim)
                {
                    searchAxis = wallNormalAxisTemp[i];
                }

                GlobalPosition global = element.geometry().center();
                global -= wallPositions[i];
                // second and argument ensures to use only aligned elements
                if (abs(global[searchAxis]) < wallDistance_[elementIdx]
                    && abs(global[searchAxis]) < global.two_norm() + 1e-8
                    && abs(global[searchAxis]) > global.two_norm() - 1e-8)
                {
                    wallDistance_[elementIdx] = abs(global[searchAxis]);
                    wallElementIdx_[elementIdx] = wallElements[i];
                    wallNormalAxis_[elementIdx] = searchAxis;
                    sandGrainRoughness_[elementIdx] = asImp_().sandGrainRoughnessAtPos(wallPositions[i]);
                }
            }
        }

        // search for neighbor Idxs
        for (const auto& element : elements(gridView))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                neighborIdx_[elementIdx][dimIdx][0] = elementIdx;
                neighborIdx_[elementIdx][dimIdx][1] = elementIdx;
            }

            for (const auto& intersection : intersections(gridView, element))
            {
                if (intersection.boundary())
                    continue;

                unsigned int neighborIdx = this->fvGridGeometry().elementMapper().index(intersection.outside());
                for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(cellCenter_[elementIdx][dimIdx] - cellCenter_[neighborIdx][dimIdx]) > 1e-8)
                    {
                        if (cellCenter_[elementIdx][dimIdx] > cellCenter_[neighborIdx][dimIdx])
                        {
                            neighborIdx_[elementIdx][dimIdx][0] = neighborIdx;
                        }
                        if (cellCenter_[elementIdx][dimIdx] < cellCenter_[neighborIdx][dimIdx])
                        {
                            neighborIdx_[elementIdx][dimIdx][1] = neighborIdx;
                        }
                    }
                }
            }
        }
    }

private:
    bool calledUpdateStaticWallProperties = false;
    std::vector<unsigned int> wallElementIdx_;
    std::vector<Scalar> wallDistance_;
    std::vector<std::array<std::array<unsigned int, 2>, dim>> neighborIdx_;
    std::vector<GlobalPosition> cellCenter_;
    std::vector<DimVector> velocity_;
    std::vector<DimVector> velocityMaximum_;
    std::vector<DimVector> velocityMinimum_;
    std::vector<DimMatrix> velocityGradients_;
    std::vector<Scalar> stressTensorScalarProduct_;
    std::vector<Scalar> vorticityTensorScalarProduct_;
    std::vector<unsigned int> wallNormalAxis_;
    std::vector<unsigned int> flowNormalAxis_;
    std::vector<Scalar> kinematicViscosity_;
    std::vector<Scalar> sandGrainRoughness_;

};

} // end namespace Dumux

#endif
