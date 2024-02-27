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
 * \ingroup BoubdaryCoupling
 * \ingroup BoxModel
 * \copydoc Dumux::PNMStokesCouplingManager
 */

#ifndef DUMUX_DROP_INTERSECTION_HH
#define DUMUX_DROP_INTERSECTION_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dumux/mixeddimension/boundary/pnmstokes/geometry.hh>

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionBoundary
 * \brief Coupling manager for drop coupled at the boundary to the bulk and low-dim
 *        domain.
 */
template <class Scalar, class GlobalPosition>
class DropIntersection
{

public:

    static auto distancePointToLine(const GlobalPosition& point1, const GlobalPosition& point2, const GlobalPosition& dropCenter)
    {
        static constexpr auto dimWorld = GlobalPosition::dimension;
        auto firstCornerPosition = point1;
        auto secondCornerPosition = point2;

        GlobalPosition directionVectorAB;
        GlobalPosition directionVectorAD;
        GlobalPosition directionVectorBD;

        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            directionVectorAB[dimIdx] = secondCornerPosition[dimIdx] - firstCornerPosition[dimIdx];
            directionVectorAD[dimIdx] = dropCenter[dimIdx] - firstCornerPosition[dimIdx];
            directionVectorBD[dimIdx] = dropCenter[dimIdx] - secondCornerPosition[dimIdx];
        }

        Scalar runningParamNom = 0.0;
        Scalar runningParamDnom = 0.0;
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            runningParamNom += directionVectorAB[dimIdx] * directionVectorAD[dimIdx];
            runningParamDnom += directionVectorAB[dimIdx] * directionVectorAB[dimIdx];
        }

        Scalar runningParam = runningParamNom/runningParamDnom;
        Scalar distance = 0.0;

        if (runningParam == 0 || runningParam <0)
        {
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                distance += directionVectorAD[dimIdx] * directionVectorAD[dimIdx];
            }
        }
        else if (runningParam == 1 || runningParam >1)
        {
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                distance += directionVectorBD[dimIdx] * directionVectorBD[dimIdx];
            }
        }
        else
        {
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                Scalar temp = 0.0;
                temp += directionVectorAD[dimIdx];
                temp -= runningParam * directionVectorAB[dimIdx];
                distance += temp * temp;
            }
        }

        distance = std::sqrt(distance);

        return distance;
    }

    template<class SubControlVolumeFace>
    static auto distancePointToSurface(const SubControlVolumeFace& scvf, const GlobalPosition& dropCenter)
    {
        static constexpr auto dimWorld = 3;
        const auto surfaceUnitNormal = scvf.unitOuterNormal();
        Scalar distance = 0.0;

        int count = 0;
        std::vector<Scalar> maxValues(dimWorld);
        std::vector<Scalar> minValues(dimWorld);

        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            maxValues[dimIdx] = scvf.corner(0)[dimIdx];
            minValues[dimIdx] = scvf.corner(0)[dimIdx];

            for (int i = 0; i<4; i++)
            {
                maxValues[dimIdx] = std::max(maxValues[dimIdx], scvf.corner(i)[dimIdx]);
                minValues[dimIdx] = std::min(minValues[dimIdx], scvf.corner(i)[dimIdx]);
            }
            if(surfaceUnitNormal[dimIdx])
            {
                Scalar temp = std::abs(dropCenter[dimIdx] - scvf.center()[dimIdx]);
                temp *= temp;
                distance += temp;
            }
            else
            {
                Scalar temp = 0.0;
                if(dropCenter[dimIdx] > maxValues[dimIdx])
                    temp = std::abs(dropCenter[dimIdx] - maxValues[dimIdx]);
                else if(dropCenter[dimIdx] < minValues[dimIdx])
                    temp = std::abs(dropCenter[dimIdx] - minValues[dimIdx]);

                temp *= temp;
                distance += temp;
            }
        }
        distance = std::sqrt(distance);
        return distance;
    }

    static auto distancePointToPoint(const GlobalPosition& point, const GlobalPosition& dropCenter)
    {
        static constexpr auto dimWorld = GlobalPosition::dimension;

        Scalar distance = 0.0;

        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            Scalar temp = 0.0;
            temp += point[dimIdx] - dropCenter[dimIdx];
            distance += temp*temp;
        }
        distance = std::sqrt(distance);

        return distance;
    }

    template<class SubControlVolumeFace>
    static auto dropToBulkFaceProjectedArea(const GlobalPosition& dropCenter, const Scalar& dropRadius, const SubControlVolumeFace& scvf)
    {
        Scalar area = 0.0;
        static constexpr auto dimWorld = GlobalPosition::dimension;

        if(dimWorld == 3)
            area = dropToBulkFaceProjectedAreaThreeDim_(dropCenter, dropRadius, scvf);
        else
            area = dropToBulkFaceProjectedAreaTwoDim_(dropCenter, dropRadius, scvf);

        return area;
    }

    template<class SubControlVolumeFace>
    static auto dropToBulkFaceProjectedAreaExtrusionFactor(const GlobalPosition& dropCenter, const Scalar& dropRadius, const SubControlVolumeFace& scvf)
    {
        auto circleRadius = dropRadius;//projectedCircleRadius(dropCenter, dropRadius, scvf);
        auto circleCenter = dropCenter;//projectedCircleCenter(dropCenter, scvf);
        if (circleRadius == 0.0)
            return 0.0;
        auto scvfCenter = scvf.center();
        auto onDropNormalVector = dropCenter;
        onDropNormalVector -= scvfCenter; // it is not unit normal


        static constexpr auto dimWorld = GlobalPosition::dimension;
        const auto& faceUnitOuterNormal = scvf.unitOuterNormal();
        const auto& faceCenter = scvf.center();



        auto distance0 = distancePointToPoint(scvf.corner(0), circleCenter);
        auto distance1 = distancePointToPoint(scvf.corner(1), circleCenter);
        auto intersectingPoint0 = scvf.corner(0);
        auto intersectingPoint1 = scvf.corner(0);

        if (distance0 < dropRadius && distance1 < dropRadius)
        {
            intersectingPoint0 = scvf.corner(0);
            intersectingPoint1 = scvf.corner(1);
        }
        else
        {

            auto maxEntities = scvf.corner(0);
            auto minEntities = scvf.corner(0);
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                maxEntities[dimIdx] = std::max(scvf.corner(0)[dimIdx], scvf.corner(1)[dimIdx]);
                minEntities[dimIdx] = std::min(scvf.corner(0)[dimIdx], scvf.corner(1)[dimIdx]);
            }



            if (faceUnitOuterNormal[0])
            {

                Scalar temp = circleRadius * circleRadius;
                temp -= (maxEntities[0] - circleCenter[0]) * (maxEntities[0] - circleCenter[0]) ;
                temp = std::sqrt(temp);

                intersectingPoint0[1] = (circleCenter[1] + temp);
                intersectingPoint1[1] = (circleCenter[1] - temp);
            }
            else
            {
                Scalar temp = circleRadius * circleRadius;
                temp -= (maxEntities[1] - circleCenter[1]) * (maxEntities[1] - circleCenter[1]) ;
                temp = std::sqrt(temp);

                intersectingPoint0[0] = (circleCenter[0] + temp);
                intersectingPoint1[0] = (circleCenter[0] - temp);
            }
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {

                // if(onDropNormalVector[0]>0)
                // {
                intersectingPoint0[dimIdx] = intersectingPoint0[dimIdx] < maxEntities[dimIdx] ? intersectingPoint0[dimIdx]:maxEntities[dimIdx];
                intersectingPoint0[dimIdx] = intersectingPoint0[dimIdx] > minEntities[dimIdx] ? intersectingPoint0[dimIdx]:minEntities[dimIdx];

                intersectingPoint1[dimIdx] = intersectingPoint1[dimIdx] < maxEntities[dimIdx] ? intersectingPoint1[dimIdx]:maxEntities[dimIdx];
                intersectingPoint1[dimIdx] = intersectingPoint1[dimIdx] > minEntities[dimIdx] ? intersectingPoint1[dimIdx]:minEntities[dimIdx];

            }
        }

        auto h0 = (intersectingPoint0[1] - circleCenter[1]);
        auto h1 = (intersectingPoint1[1] - circleCenter[1]);
        if ((h0 - h1 == 0.0))
            return 0.0;
        else if (h0*h1 > 0)
        {
            h0 = std::abs(h0);
            h1 = std::abs(h1);
            auto tetha0 = h0/circleRadius<1.0 ? 2*std::acos(h0/circleRadius): 0.0;
            auto tetha1 = h1/circleRadius<1.0 ? 2*std::acos(h1/circleRadius): 0.0;

            auto A0 = 0.5 * (tetha0 - std::sin(tetha0)) * dropRadius * dropRadius;
            auto A1 = 0.5 * (tetha1 - std::sin(tetha1)) * dropRadius * dropRadius;

            return std::abs(A0 - A1);
        }
        else
        {
            h0 = std::abs(h0);
            h1 = std::abs(h1);
            auto tetha0 = h0/circleRadius<1.0 ? 2*std::acos(h0/circleRadius): 0.0;
            auto tetha1 = h1/circleRadius<1.0 ? 2*std::acos(h1/circleRadius): 0.0;
            auto A0 = 0.5 *(M_PI -  (tetha0 - std::sin(tetha0))) * dropRadius * dropRadius;
            auto A1 = 0.5 * (M_PI - (tetha1 - std::sin(tetha1))) * dropRadius * dropRadius;

            return std::abs(A0 + A1);

        }
    }
private:

    template<class SubControlVolumeFace>
    static auto dropToBulkFaceProjectedAreaThreeDim_(const GlobalPosition& dropCenter, const Scalar& dropRadius, const SubControlVolumeFace& scvf)
    {
        auto circleRadius = projectedCircleRadius_(dropCenter, dropRadius, scvf);
        auto circleCenter = projectedCircleCenter_(dropCenter, scvf);


        static constexpr auto dimWorld = GlobalPosition::dimension;
        const auto& faceUnitOuterNormal = scvf.unitOuterNormal();
        const auto& faceCenter = scvf.center();

        std::vector<Scalar> distances(4);
        int count = 0;
        for (int i = 0; i < 4; i++)
        {
            distances[i] = distancePointToPoint(scvf.corner(i), circleCenter);
            if (distances[i] <= circleRadius)
                count +=1;
        }

        if (count == 4)
            return scvf.area();

        auto maxEntities = scvf.corner(0);
        auto minEntities = scvf.corner(0);
        for (size_t i = 0; i<4; i++)
        {
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                maxEntities[dimIdx] = std::max(maxEntities[dimIdx], scvf.corner(i)[dimIdx]);
                minEntities[dimIdx] = std::min(minEntities[dimIdx], scvf.corner(i)[dimIdx]);
            }
        }

        std::vector<GlobalPosition> intersectingPoints;
        for (int dimIdx = 0; dimIdx < dimWorld; dimIdx++)
        {
            if (faceUnitOuterNormal[dimIdx])
                continue;

            GlobalPosition intersectingPointMax;
            for  (int i = 0; i < dimWorld; i++)
            {
                if (faceUnitOuterNormal[i])
                    intersectingPointMax[i] = faceCenter[i];
                else
                    intersectingPointMax[i] = -1.0;
            }

            intersectingPointMax[dimIdx] = maxEntities[dimIdx];

            Scalar temp = circleRadius * circleRadius;
            int missingIndex;
            for (int i = 0; i < dimWorld; i++)
            {
                if (intersectingPointMax[i] == -1.0)
                {
                    missingIndex = i;
                    continue;
                }

                temp -= (intersectingPointMax[i] - circleCenter[i]) * (intersectingPointMax[i] - circleCenter[i]);
            }
            if (temp < 0)
                continue;
            temp = std::sqrt(temp);
            if (faceCenter[missingIndex] < circleCenter[missingIndex])
                temp *= -1.0;

            intersectingPointMax[missingIndex] = circleCenter[missingIndex];
            intersectingPointMax[missingIndex] += temp;

            intersectingPoints.push_back(intersectingPointMax);
        }


        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            if (faceUnitOuterNormal[dimIdx])
                continue;
            GlobalPosition intersectingPointMin;
            for  (int i = 0; i < dimWorld; i++)
            {
                if (faceUnitOuterNormal[i])
                    intersectingPointMin[i] = faceCenter[i];
                else
                    intersectingPointMin[i] = -1.0;
            }

            intersectingPointMin[dimIdx] = minEntities[dimIdx];

            Scalar temp = circleRadius * circleRadius;
            int missingIndex;
            for (int i = 0; i < dimWorld; ++i)
            {
                if (intersectingPointMin[i] == -1.0)
                {
                    missingIndex = i;
                    continue;
                }

                temp -= (intersectingPointMin[i] - circleCenter[i]) * (intersectingPointMin[i] - circleCenter[i]);
            }
            if (temp < 0)
                continue;
            temp = std::sqrt(temp);
            if (faceCenter[missingIndex] < circleCenter[missingIndex])
                temp *= -1.0;

            intersectingPointMin[missingIndex] = circleCenter[missingIndex];
            intersectingPointMin[missingIndex] += temp;

            intersectingPoints.push_back(intersectingPointMin);
        }

        std::vector<GlobalPosition> inBoundIntersectingPoints;
        for (auto& intersectingPoint:intersectingPoints)
        {
            bool isInBox = true;
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                if (intersectingPoint[dimIdx] > maxEntities[dimIdx] || intersectingPoint[dimIdx] < minEntities[dimIdx])
                {
                    isInBox = false;
                }
            }
            if(isInBox)
                inBoundIntersectingPoints.push_back(intersectingPoint);
        }

        if (inBoundIntersectingPoints.size() > 2)
            return scvf.area();
        else if (inBoundIntersectingPoints.size() < 2)
            return 0.0;

        auto maxEntitiesPoints = inBoundIntersectingPoints[0];
        auto minEntitiesPoints = inBoundIntersectingPoints[0];

        for (auto &intersectingPoint : inBoundIntersectingPoints)
        {
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                if (!faceUnitOuterNormal[dimIdx])
                {
                    maxEntitiesPoints[dimIdx] = std::max(maxEntitiesPoints[dimIdx], intersectingPoint[dimIdx]);
                    minEntitiesPoints[dimIdx] = std::min(minEntitiesPoints[dimIdx], intersectingPoint[dimIdx]);
                }
            }
        }

        std::vector<Scalar> PointToFaceBound(dimWorld, 0.0);
        std::vector<Scalar> PointToPointBound(dimWorld, 0.0);
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            PointToFaceBound[dimIdx] = faceCenter[dimIdx] > circleCenter[dimIdx] ? (minEntitiesPoints[dimIdx] - minEntities[dimIdx]) : (maxEntitiesPoints[dimIdx] - maxEntities[dimIdx]);
            PointToPointBound[dimIdx] = maxEntitiesPoints[dimIdx] - minEntitiesPoints[dimIdx];
        }

        Scalar area = 0.0;
        Scalar temp0 = 1.0;
        for (int i = 0; i < dimWorld; i++)
        {
            if (faceUnitOuterNormal[i])
                continue;
            Scalar temp1 = 1.0;
            for (int j = 0; j < dimWorld; j++)
            {
                if (faceUnitOuterNormal[j])
                    continue;
                else if (j == i)
                    temp1 *= PointToPointBound[j];
                else
                    temp1 *= PointToFaceBound[j];
            }
            temp0 *= PointToFaceBound[i];

            area += std::abs(temp1);
        }
        area += std::abs(temp0);

        int firstIntegralIndex = -1;
        int secondIntegralIndex = -1;
        for (int i = 0; i < dimWorld; i++)
        {
            if (faceUnitOuterNormal[i])
                continue;
            if (maxEntitiesPoints[i] == minEntitiesPoints[i])
                secondIntegralIndex = i;
            else
            {
                if (firstIntegralIndex == -1)
                    firstIntegralIndex = i;
                else
                    secondIntegralIndex = i;
            }
        }

        auto integralFirstPart = [&](Scalar value) {
        Scalar integralResult = 0.0;

        if (value > circleRadius)
            value = circleRadius;
        else if(value<(-circleRadius))
            value = -circleRadius;

        integralResult += circleRadius * circleRadius * std::asin(value/ circleRadius) / 2;
        integralResult += value * std::sqrt(circleRadius * circleRadius - value * value) / 2;

        return integralResult;

        };

        Scalar segmentArea = 0.0;
        Scalar segmentArea2 = 0.0;
        segmentArea += integralFirstPart(maxEntitiesPoints[firstIntegralIndex] - circleCenter[firstIntegralIndex]);
        segmentArea -= integralFirstPart(minEntitiesPoints[firstIntegralIndex] - circleCenter[firstIntegralIndex]);

        if (faceCenter[secondIntegralIndex] > circleCenter[secondIntegralIndex])
        {
            segmentArea2 += PointToPointBound[firstIntegralIndex] * (minEntitiesPoints[secondIntegralIndex] - circleCenter[secondIntegralIndex]);
        }
        else
        {
            segmentArea2 += PointToPointBound[firstIntegralIndex] * (maxEntitiesPoints[secondIntegralIndex] - circleCenter[secondIntegralIndex]);
        }

        area += std::abs(segmentArea);
        area -= std::abs(segmentArea2);

        return area;
    }

    template<class SubControlVolumeFace>
    static auto dropToBulkFaceProjectedAreaTwoDim_(const GlobalPosition& dropCenter, const Scalar& dropRadius, const SubControlVolumeFace& scvf)
    {
        auto circleRadius = dropRadius;//projectedCircleRadius(dropCenter, dropRadius, scvf);
        auto circleCenter = dropCenter;//projectedCircleCenter(dropCenter, scvf);
        auto scvfCenter = scvf.center();
        auto onDropNormalVector = dropCenter;
        onDropNormalVector -= scvfCenter; // it is not unit normal


        static constexpr auto dimWorld = GlobalPosition::dimension;
        const auto& faceUnitOuterNormal = scvf.unitOuterNormal();
        const auto& faceCenter = scvf.center();



            auto distance0 = distancePointToPoint(scvf.corner(0), circleCenter);
            auto distance1 = distancePointToPoint(scvf.corner(1), circleCenter);

        Scalar epsilon = dropRadius * 1e-5;
        if ((distance0 < dropRadius + epsilon)  && (distance1 < dropRadius + epsilon))
            return scvf.area();

        auto maxEntities = scvf.corner(0);
        auto minEntities = scvf.corner(0);
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            maxEntities[dimIdx] = std::max(scvf.corner(0)[dimIdx], scvf.corner(1)[dimIdx]);
            minEntities[dimIdx] = std::min(scvf.corner(0)[dimIdx], scvf.corner(1)[dimIdx]);
        }

        auto intersectingPoint0 = scvf.corner(0);
        auto intersectingPoint1 = scvf.corner(0);

        if (faceUnitOuterNormal[0])
        {

            Scalar temp = circleRadius * circleRadius;
            temp -= (maxEntities[0] - circleCenter[0]) * (maxEntities[0] - circleCenter[0]) ;
            temp = std::sqrt(temp);

            intersectingPoint0[1] = (circleCenter[1] + temp);
            intersectingPoint1[1] = (circleCenter[1] - temp);
        }
        else
        {
            Scalar temp = circleRadius * circleRadius;
            temp -= (maxEntities[1] - circleCenter[1]) * (maxEntities[1] - circleCenter[1]) ;
            temp = std::sqrt(temp);

            intersectingPoint0[0] = (circleCenter[0] + temp);
            intersectingPoint1[0] = (circleCenter[0] - temp);
        }
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {

            // if(onDropNormalVector[0]>0)
            // {
                intersectingPoint0[dimIdx] = intersectingPoint0[dimIdx] < maxEntities[dimIdx] ? intersectingPoint0[dimIdx]:maxEntities[dimIdx];
                intersectingPoint0[dimIdx] = intersectingPoint0[dimIdx] > minEntities[dimIdx] ? intersectingPoint0[dimIdx]:minEntities[dimIdx];

                intersectingPoint1[dimIdx] = intersectingPoint1[dimIdx] < maxEntities[dimIdx] ? intersectingPoint1[dimIdx]:maxEntities[dimIdx];
                intersectingPoint1[dimIdx] = intersectingPoint1[dimIdx] > minEntities[dimIdx] ? intersectingPoint1[dimIdx]:minEntities[dimIdx];

        }
        return (intersectingPoint1 - intersectingPoint0).two_norm();
    }

    template <class SubControlVolumeFace>
    static auto projectedCircleRadius_(const GlobalPosition &dropCenter, const Scalar &dropRadius, const SubControlVolumeFace &scvf)
    {
        static constexpr auto dimWorld = GlobalPosition::dimension;
        Scalar distanceToCenter;
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            if (scvf.unitOuterNormal()[dimIdx])
                distanceToCenter = dropCenter[dimIdx] - scvf.center()[dimIdx];
        }

        Scalar circleRadius = dropRadius * dropRadius;
        circleRadius -= distanceToCenter * distanceToCenter;

        return std::sqrt(circleRadius);
    }

    template <class SubControlVolumeFace>
    static auto projectedCircleCenter_(const GlobalPosition &dropCenter, const SubControlVolumeFace &scvf)
    {
        static constexpr auto dimWorld = GlobalPosition::dimension;
        GlobalPosition circleCenter = dropCenter;
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            if (scvf.unitOuterNormal()[dimIdx])
            {
                Scalar distanceToCenter = dropCenter[dimIdx] - scvf.center()[dimIdx];
                circleCenter[dimIdx] -= distanceToCenter;
            }
        }
        return circleCenter;
    }
};

} // end namespace Dumux

#endif
