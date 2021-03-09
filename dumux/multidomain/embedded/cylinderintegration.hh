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
 * \ingroup EmbeddedCoupling
 * \brief Integration over cylindrical and elliptic cylindrical domains
 * Lowest order integration formulas that are mostly useful to evenly distribute
 * mass in a cylindrical domain.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_CYLINDERINTEGRATOR_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_CYLINDERINTEGRATOR_HH

#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <utility>
#include <tuple>
#include <limits>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>

namespace Dumux::EmbeddedCoupling {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Helper class to integrate over a cylinder domain
 * \note This is mostly useful if the integral is known and the
 *       integral mass has to be locally distributed to some non-matching domain
 * \note The algorithm creates (almost) evenly spaced area elements on a circle
 *       and extrudes this pattern into a cylinder. Each sub-cell is an integration
 *       point with constant ansatz on the sub-cell (lowest order integration formula).
 *       Hence to improve the quality of the integral, increase the sample size.
 * \note The volume integral of a constant function is exact
 * See Beckers & Beckers (2012) doi:10.1016/j.comgeo.2012.01.011 for circle distribution
 * and T. Koch PhD thesis (2020), Section 7.2.4.
 */
template<class Scalar>
class CylinderIntegration
{
    using GlobalPosition = Dune::FieldVector<Scalar, 3>;

public:
    /*!
     * \brief Constructor
     * \param rStep characteristic relative integration length (positive real number between 0 and 1)
     * \note half the characteristic length means 2*2*2=8 times more integration points
     */
    CylinderIntegration(const Scalar rStep)
    {
        rSamples_ = characteristicLengthToNumSamples_(rStep);
        initialize_();
    }

    /*!
     * \brief Constructor
     * \param rStep characteristic relative integration length in r-direction (positive real number between 0 and 1)
     * \param zStep characteristic relative integration length in z-direction (positive real number between 0 and 1)
     * \note Use this constructor to achieve a non-balanced (away from 1) aspect ratio between r and z-direction
     */
    CylinderIntegration(const Scalar rStep, const Scalar zStep)
    : zStepFixed_(true)
    {
        using std::ceil; using std::clamp;
        rSamples_ = characteristicLengthToNumSamples_(rStep);
        zSamples_ = characteristicLengthToNumSamples_(zStep);
        initialize_();
    }

    /*!
     * \brief Set the geometry of the cylinder
     * \param bottomCenter bottom center position
     * \param topCenter top center position
     * \param radius cylinder radius
     * \param verbosity the verbosity level (default: 0 -> no terminal output)
     */
    void setGeometry(const GlobalPosition& bottomCenter,
                     const GlobalPosition& topCenter,
                     const Scalar radius,
                     int verbosity = 0)
    {
        const auto zAxis = (topCenter-bottomCenter);
        const auto height = zAxis.two_norm();
        const auto zAxisUnitVector = zAxis/height;

        // compute step size in r-direction
        const auto rStep = radius/Scalar(rSamples_);

        // compute zSamples from r samples if not specified by the user
        using std::ceil;
        if (!zStepFixed_)
            zSamples_ = std::max<std::size_t>(1, ceil(height/rStep));
        const auto zStep = height/Scalar(zSamples_);

        // compute offsets for index calculation
        auto kOffset = numRingSamples_;
        std::partial_sum(kOffset.begin(), kOffset.end(), kOffset.begin());

        // compute total number of samples
        numPoints_ = zSamples_*circleSamples_;
        if (verbosity > 0)
            std::cout << "CylinderIntegration -- number of integration points: " << numPoints_ << std::endl;

        // resize integration points and integration elements
        integrationElement_.resize(numPoints_);
        points_.resize(numPoints_);

        // reserve enough memory for the points on each ring
        std::vector<GlobalPosition> ringPoints;
        ringPoints.reserve(numRingSamples_.back());

        // compute integration points and integration elements
        for (std::size_t i = 0; i < zSamples_; ++i)
        {
            // for each cylinder slice i
            auto sliceCenter = bottomCenter;
            sliceCenter.axpy((Scalar(i)+0.5)*zStep, zAxisUnitVector);

            // generate circle points for each ring j
            for (std::size_t j = 0; j < rSamples_; ++j)
            {
                const auto r = (Scalar(j)+0.5)*rStep;
                circlePoints(ringPoints, sincos_[j], sliceCenter, zAxisUnitVector, r);
                const auto ringOffest = j > 0 ? kOffset[j-1] : 0;
                const auto ringSamples = numRingSamples_[j];
                // volume of each element in ring
                // total ring volume is given by M_PI*rStep*rStep*zStep*((j+1)^2 - j^2)
                const auto integrationElement = M_PI*rStep*rStep*zStep*(1.0 + 2.0*Scalar(j))/ringSamples;
                for (std::size_t k = 0; k < ringSamples; ++k)
                {
                    const std::size_t idx = k + ringOffest + circleSamples_*i;
                    points_[idx] = ringPoints[k];
                    integrationElement_[idx] = integrationElement;
                }
            }
        }
    }

    //! The integration element of the ith integration point
    Scalar integrationElement(std::size_t i) const
    { return integrationElement_[i]; }

    //! The ith integration point
    const GlobalPosition& integrationPoint(std::size_t i) const
    { return points_[i]; }

    //! The number of integration points
    std::size_t size() const
    { return numPoints_; }

private:
    void initialize_()
    {
        // precompute number of cells in ring and total circle samples
        circleSamples_ = 0;
        numRingSamples_.resize(rSamples_);
        for (int j = 0; j < rSamples_; ++j)
        {
            // number of cells in ring j
            using std::floor;
            numRingSamples_[j] = floor(2*M_PI*(Scalar(j)+0.5));
            circleSamples_ += numRingSamples_[j];
        }

        // optimization because calling too many sin/cos can be really expensive
        sincos_.resize(rSamples_);
        for (int j = 0; j < rSamples_; ++j)
        {
            const auto numPoints = numRingSamples_[j];
            sincos_[j].resize(2*numPoints);
            Scalar t = 0 + 0.1;
            for (std::size_t i = 0; i < numPoints; ++i)
            {
                using std::sin; using std::cos;
                sincos_[j][2*i] = sin(t);
                sincos_[j][2*i + 1] = cos(t);
                t += 2*M_PI/numPoints;
                if (t > 2*M_PI) t -= 2*M_PI;
            }
        }
    }

    // invert characteristic length (between 0 and 1)
    // and make sure there is no integer wrap around for small ratios
    std::size_t characteristicLengthToNumSamples_(const Scalar cl) const
    {
        using std::clamp; using std::min; using std::max; using std::ceil; using std::nexttoward;
        const Scalar clampedCl = clamp(cl, std::numeric_limits<Scalar>::min(), 1.0);
        const Scalar floatNumSamples = ceil(1.0/clampedCl);
        const Scalar largestValidFloat = nexttoward(std::numeric_limits<std::size_t>::max(), 0.0);
        return static_cast<std::size_t>(min(floatNumSamples, largestValidFloat));
    }

    bool zStepFixed_ = false;
    std::size_t zSamples_, rSamples_, numPoints_, circleSamples_;
    std::vector<Scalar> integrationElement_;
    std::vector<GlobalPosition> points_;
    std::vector<std::size_t> numRingSamples_;
    std::vector<std::vector<Scalar>> sincos_;
};

namespace Detail {
//! check if a point is in an ellipse
template<class GlobalPosition>
inline bool pointInEllipse(const GlobalPosition& p,
                           const GlobalPosition& center,
                           const GlobalPosition& firstAxis,
                           const GlobalPosition& secondAxis,
                           const GlobalPosition& normal,
                           const typename GlobalPosition::value_type a,
                           const typename GlobalPosition::value_type b)
{
    const auto d = p-center;
    // early return if point is not on the ellipse plane
    if (d*normal > 1e-7*a)
        return false;

    // check if it's actually within the ellipse
    const auto da = (d*firstAxis);
    const auto db = (d*secondAxis);

    return (da*da/(a*a) + db*db/(b*b) < 1.0);
}

//! construct evenly distributed integration points on an ellipse
template<class GlobalPosition>
inline std::pair<std::vector<GlobalPosition>, typename GlobalPosition::value_type>
ellipseIntegrationPoints(const GlobalPosition& center,
                         const GlobalPosition& firstUnitAxis,
                         const GlobalPosition& secondUnitAxis,
                         typename GlobalPosition::value_type a,
                         typename GlobalPosition::value_type b,
                         const GlobalPosition& normal,
                         typename GlobalPosition::value_type characteristicLength)
{
    // choose step size in a-direction
    const auto aStep = characteristicLength;

    using std::floor;
    const std::size_t aSamples = std::max<std::size_t>(floor(2*a/aStep), 1);
    const std::size_t bSamples = std::max<std::size_t>(floor(2*b/aStep), 1);
    const auto bStep = 2*b / bSamples;

    // reserve upper limit estimate for memory needed
    std::vector<GlobalPosition> points;
    points.reserve(aSamples*bSamples);

    // walk over lattice grid and reject points outside the ellipse
    auto startAB = center;
    startAB.axpy(-a + 0.5*aStep, firstUnitAxis);
    startAB.axpy(-b + 0.5*bStep, secondUnitAxis);
    for (std::size_t as = 0; as < aSamples; ++as)
    {
        auto posA = startAB;
        posA.axpy(as*aStep, firstUnitAxis);
        for (std::size_t bs = 0; bs < bSamples; ++bs)
        {
            auto pos = posA;
            pos.axpy(bs*bStep, secondUnitAxis);
            if (pointInEllipse(pos, center, firstUnitAxis, secondUnitAxis, normal, a, b))
                points.emplace_back(std::move(pos));
        }
    }
    // return points and integration element
    return {points, aStep*bStep};
}
} // end namespace Detail

/*!
 * \ingroup EmbeddedCoupling
 * \brief Helper class to integrate over an elliptic cylinder domain
 * \note The cylinder caps do not have to be orthogonal to the centerline axis but top and bottom cap are parallel planes
 * \note This is mostly useful if the integral is known and the
 *       integral mass has to be locally distributed to some non-matching domain
 * \note The algorithm is based on evenly spaced area elements with rejection sampling.
 *       To improve the quality of the integral, increase the sample size.
 * \note The volume integral of a constant function is exact
 * See Koch et al. (2020) https://doi.org/10.1016/j.jcp.2020.109369
 */
template<class Scalar>
class EllipticCylinderIntegration
{
    using GlobalPosition = Dune::FieldVector<Scalar, 3>;

public:
    /*!
     * \brief Constructor
     * \param relCharLength characteristic relative integration length (real number between 0 and 1)
     * \note half the characteristic length means 2*2*2=8 times more integration points
     */
    explicit EllipticCylinderIntegration(const Scalar relCharLength)
    : relCharLength_(relCharLength)
    {}

    /*!
     * \brief Set the geometry of elliptic cylinder
     * \param bottomCenter bottom center position
     * \param topCenter top center position
     * \param firstAxis first ellipse axis (length corresponding to axis length)
     * \param secondAxis second ellipse axis (length corresponding to axis length)
     * \param verbosity the verbosity level (default: 0 -> no terminal output)
     */
    void setGeometry(const GlobalPosition& bottomCenter,
                     const GlobalPosition& topCenter,
                     const GlobalPosition& firstAxis,
                     const GlobalPosition& secondAxis,
                     int verbosity = 0)
    {
        const auto a = firstAxis.two_norm();
        const auto b = secondAxis.two_norm();
        const auto characteristicLength = relCharLength_*a;

        auto firstUnitAxis = firstAxis; firstUnitAxis /= a;
        auto secondUnitAxis = secondAxis; secondUnitAxis /= b;
        const auto normal = crossProduct(firstUnitAxis, secondUnitAxis);

        auto zAxis = (topCenter-bottomCenter);
        const auto length = zAxis.two_norm();
        zAxis /= length;
        const auto height = (zAxis*normal)*length;
        using std::floor;
        const std::size_t zSamples = std::max<std::size_t>(floor(height/characteristicLength), 1);
        const auto zStep = length / Scalar(zSamples);
        const auto zStepFactor = length/height;

        // walk over lattice grid and reject points outside the ellipse
        auto startZ = bottomCenter;
        startZ.axpy(0.5*zStep, zAxis);
        auto [ellipseIPs, ellipseIntegrationElement]
            = Detail::ellipseIntegrationPoints(startZ, firstUnitAxis, secondUnitAxis, a, b, normal, characteristicLength);

        // compute total number of points
        const auto abPoints = ellipseIPs.size();
        numPoints_ = abPoints*zSamples;

        // first move in the first layer of ellipse integration points
        points_.clear();
        points_.reserve(numPoints_);
        std::move(ellipseIPs.begin(), ellipseIPs.end(), std::back_inserter(points_));

        // then extrude along the z axis
        auto zStepVector = zAxis; zStepVector *= zStep;
        for (std::size_t zs = 1; zs < zSamples; ++zs)
            std::transform(points_.end() - abPoints, points_.end(), std::back_inserter(points_),
                           [&](const auto& p){ return p + zStepVector; });

        // computation and error correction for integral element to obtain exact integral of constant functions
        integrationElement_ = ellipseIntegrationElement*zStep/zStepFactor;
        const auto meanLocalError = (height*a*b*M_PI - integrationElement_*numPoints_)/Scalar(numPoints_);
        integrationElement_ += meanLocalError;

        if (verbosity > 0)
        {
            std::cout << "EllipticCylinderIntegration -- number of integration points: " << numPoints_ <<"\n";
            if (verbosity > 1)
                std::cout << "  -- a: " << a << ", b: " << b << ", h: " << height << "\n"
                          << "  -- volume: " << integrationElement_*numPoints_ << " (expected " << height*a*b*M_PI << ")\n"
                          << "  -- corrected a volume error of: " << meanLocalError/integrationElement_*100 << "%\n";
        }
    }

    //! obtain ith integration element
    Scalar integrationElement(std::size_t i) const
    { return integrationElement_; }

    //! obtain ith integration point
    const GlobalPosition& integrationPoint(std::size_t i) const
    { return points_[i]; }

    //! number of samples points
    std::size_t size() const
    { return numPoints_; }

private:
    const Scalar relCharLength_;
    std::size_t numPoints_;
    Scalar integrationElement_;
    std::vector<GlobalPosition> points_;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief Helper class to integrate over an elliptic domain
 * \note This is mostly useful if the integral is known and the
 *       integral mass has to be locally distributed to some non-matching domain
 * \note The algorithm is based on evenly spaced area elements with rejection sampling.
 *       To improve the quality of the integral, increase the sample size.
 * \note The area integral of a constant function is exact
 * See Koch et al. (2020) https://doi.org/10.1016/j.jcp.2020.109369
 */
template<class Scalar>
class EllipseIntegration
{
    using GlobalPosition = Dune::FieldVector<Scalar, 3>;

public:
    /*!
     * \brief Constructor
     * \param relCharLength characteristic relative integration length (real number between 0 and 1)
     * \note half the characteristic length means 2*2=4 times more integration points
     */
    explicit EllipseIntegration(const Scalar relCharLength)
    : relCharLength_(relCharLength)
    {}

    /*!
     * \brief set geometry of an ellipse
     * \param center the center position
     * \param firstAxis first ellipse axis (length corresponding to axis length)
     * \param secondAxis second ellipse axis (length corresponding to axis length)
     * \param verbosity the verbosity level (default: 0 -> no terminal output)
     */
    void setGeometry(const GlobalPosition& center,
                     const GlobalPosition& firstAxis,
                     const GlobalPosition& secondAxis,
                     int verbosity = 0)
    {
        const auto a = firstAxis.two_norm();
        const auto b = secondAxis.two_norm();
        const auto characteristicLength = relCharLength_*a;

        auto firstUnitAxis = firstAxis; firstUnitAxis /= a;
        auto secondUnitAxis = secondAxis; secondUnitAxis /= b;
        const auto normal = crossProduct(firstUnitAxis, secondUnitAxis);

        // generate integration points
        std::tie(points_, integrationElement_)
            = Detail::ellipseIntegrationPoints(center, firstUnitAxis, secondUnitAxis, a, b, normal, characteristicLength);

        // store number of sample points
        const auto abPoints = points_.size();
        numPoints_ = abPoints;

        // computation and error correction for integral element to obtain exact integral of constant functions
        const auto meanLocalError = (a*b*M_PI - integrationElement_*numPoints_)/Scalar(numPoints_);
        integrationElement_ += meanLocalError;

        if (verbosity > 0)
        {
            std::cout << "EllipseIntegration -- number of integration points: " << numPoints_ << "\n";
            if (verbosity > 1)
                std::cout << "  -- a: " << a << ", b: " << b << "\n"
                          << "  -- area: " << integrationElement_*numPoints_ << " (expected " << a*b*M_PI << ")\n"
                          << "  -- corrected an area error of: " << meanLocalError/integrationElement_*100 << "%\n";
        }
    }

    //! obtain ith integration element
    Scalar integrationElement(std::size_t i) const
    { return integrationElement_; }

    //! obtain ith integration point
    const GlobalPosition& integrationPoint(std::size_t i) const
    { return points_[i]; }

    //! number of integration points
    std::size_t size() const
    { return numPoints_; }

private:
    const Scalar relCharLength_;
    std::size_t numPoints_;
    Scalar integrationElement_;
    std::vector<GlobalPosition> points_;
};

} // end namespace Dumux::EmbeddedCoupling

#endif
