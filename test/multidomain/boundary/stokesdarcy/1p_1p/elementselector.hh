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
 *
 * \brief Auxiliary class used to select the elements for the final grid
 */
#ifndef DUMUX_MICROMODEL_ELEMENT_SELECTOR_HH
#define DUMUX_MICROMODEL_ELEMENT_SELECTOR_HH

#include <string>
#include <dune/common/fvector.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dumux/geometry/intersectspointgeometry.hh>

namespace Dumux
{

/*!
 * \brief Auxiliary class used to select the elements for the final grid
 */
template <class GridView, bool isDarcy = false>
class ElementSelector
{
    using Scalar = typename GridView::ctype;
    static constexpr int coordDim = 2;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, coordDim>;
    using Rectangle = Dune::AxisAlignedCubeGeometry<Scalar, coordDim, coordDim>;
    using Triangle = Dune::MultiLinearGeometry<Scalar, coordDim, coordDim>;

public:
    ElementSelector(const std::string& modelParamGroup = "")
    {
        const Scalar lengthCavity = getParamFromGroup<Scalar>(modelParamGroup, "Grid.CavityLength");
        const Scalar heightCavity = getParamFromGroup<Scalar>(modelParamGroup, "Grid.CavityHeight");


        // std::cout << std::setprecision(15) << "lengthCavity " << lengthCavity << std::endl;
        // std::cout << std::setprecision(15) << "numChannelsX " << numChannelsX << std::endl;
        // std::cout << std::setprecision(15) << "channelWidthX " << channelWidthX << std::endl;
        // std::cout << std::setprecision(15) << "numPillarsX " << numPillarsX << std::endl;
        // std::cout << std::setprecision(15) << "pillarWidthX " << pillarWidthX << std::endl;

        const Scalar inletLength = getParamFromGroup<Scalar>(modelParamGroup, "Grid.InletLength");

        const auto positions0 = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.Positions0");
        const auto positions1 = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.Positions1");

        const GlobalPosition lowerLeft = GlobalPosition{positions0[0], positions1[0]};
        const GlobalPosition upperRight = GlobalPosition{positions0.back(), positions1.back()};

        Rectangle belowInlet(lowerLeft, GlobalPosition{lowerLeft[0]+inletLength, heightCavity});
        Rectangle belowOutlet(GlobalPosition{inletLength+lengthCavity+lowerLeft[0], lowerLeft[1]}, GlobalPosition{upperRight[0], heightCavity});

        const Scalar xMin = lowerLeft[0]+inletLength;
        const Scalar xMax = inletLength+lengthCavity+lowerLeft[0];
        const Scalar yMin = positions1[positions1.size()-3];
        const Scalar yMax = positions1[positions1.size()-2];
        Rectangle porousMedium(GlobalPosition{xMin, yMin}, GlobalPosition{xMax, yMax});
        rectangles_.push_back(porousMedium);

        if constexpr (!isDarcy)
        {
            rectangles_.push_back(belowInlet);
            rectangles_.push_back(belowOutlet);
        }

        const bool triangularReservoir = getParamFromGroup<bool>(modelParamGroup, "Grid.TriangularReservoir", false);
        if (triangularReservoir)
        {
            // account for slope at bottom reservoir
            const Scalar reservoirY = getParamFromGroup<Scalar>(modelParamGroup, "Grid.TriangleHeight");
            const Scalar triangleTipWidth = getParamFromGroup<Scalar>(modelParamGroup, "Grid.TriangleTipWidth");

            std::vector<GlobalPosition> pointsLeftTriangle({lowerLeft + GlobalPosition{inletLength, 0.0},
                                                            lowerLeft + GlobalPosition{inletLength + 0.5*lengthCavity - 0.5*triangleTipWidth, 0.0},
                                                            lowerLeft + GlobalPosition{inletLength, reservoirY}});

            std::vector<GlobalPosition> pointsRightTriangle({lowerLeft + GlobalPosition{inletLength + 0.5*lengthCavity + 0.5*triangleTipWidth, 0.0},
                                                             lowerLeft + GlobalPosition{inletLength + lengthCavity, 0.0},
                                                             lowerLeft + GlobalPosition{inletLength + lengthCavity, reservoirY}});

            triangles_.emplace_back(Dune::GeometryTypes::simplex(2), pointsLeftTriangle);
            triangles_.emplace_back(Dune::GeometryTypes::simplex(2), pointsRightTriangle);
        }
    }

    //! Select all elements that are not cut-out by the recangles or triangles
    bool operator() (const Element& element) const
    {
        // make the selector work for 2D and 3D
        const auto center = [&]()
        {
            // for 3D, only consider the x and y component
            const auto tmpCenter = element.geometry().center();
            return GlobalPosition{tmpCenter[0], tmpCenter[1]};
        }();

        if constexpr (!isDarcy)
        {
            for (auto&& triangle : triangles_)
            {
                if (intersectsPointGeometry(center, triangle))
                    return false;
            }

            for (auto&& rectangle : rectangles_)
            {
                if (intersectsPointGeometry(center, rectangle))
                    return false;
            }

            return true;
        }
        else
        {
            for (auto&& rectangle : rectangles_)
            {
                if (intersectsPointGeometry(center, rectangle))
                    return true;
            }

            return false;
        }


    }

private:
    std::vector<Rectangle> rectangles_;
    std::vector<Triangle> triangles_;
};

} //end namespace

#endif
