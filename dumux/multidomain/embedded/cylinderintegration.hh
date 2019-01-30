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
 * \brief Integration over a cylindrical non-matching domain
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_CYLINDERINTEGRATOR_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_CYLINDERINTEGRATOR_HH

#include <cmath>
#include <random>
#include <algorithm>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>

namespace Dumux {

enum class CylinderIntegrationMethod
{ random, spaced };

template<class Scalar,
         CylinderIntegrationMethod m = CylinderIntegrationMethod::random>
class CylinderIntegration;


template<class Scalar>
class CylinderIntegration<Scalar, CylinderIntegrationMethod::random>
{
    using GlobalPosition = Dune::FieldVector<Scalar, 3>;

public:
    CylinderIntegration(std::size_t samples)
    : samples_(samples)
    {
        // generate random numbers
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<Scalar> dist0(0.0, 1.0), dist1(-1.0, 1.0);
        randomVectors_.resize(samples);
        randomNumbers_.resize(2*samples);
        std::generate(randomVectors_.begin(), randomVectors_.end(), [&](){ return GlobalPosition({dist1(gen), dist1(gen), dist1(gen)}); });
        std::generate(randomNumbers_.begin(), randomNumbers_.end(), [&](){ return dist0(gen); });
    }

    void setGeometry(const GlobalPosition& a, const GlobalPosition& b, Scalar radius)
    {
        a_ = a;
        ab_ = (b-a);
        radius_ = radius;
        height_ = ab_.two_norm();
        integrationElement_ = M_PI*radius*radius*height_/samples_;
    }

    Scalar integrationElement(std::size_t i) const
    { return integrationElement_; }

    GlobalPosition getIntegrationPoint(std::size_t i) const
    {
        using std::sqrt;

        auto point = a_;
        point.axpy(randomNumbers_[i], ab_); // walk random step on centerline
        const auto& p = randomVectors_[i];
        auto r = crossProduct(p, ab_); r /= r.two_norm();
        const auto rho = radius_*sqrt(randomNumbers_[samples_ + i]);
        point.axpy(rho, r);

        return point;
    }

    std::size_t size() const
    { return samples_; }

private:
    std::vector<GlobalPosition> randomVectors_;
    std::vector<Scalar> randomNumbers_;
    GlobalPosition a_, ab_;
    Scalar radius_, height_;
    Scalar integrationElement_;
    std::size_t samples_;
};

template<class Scalar>
class CylinderIntegration<Scalar, CylinderIntegrationMethod::spaced>
{
    using GlobalPosition = Dune::FieldVector<Scalar, 3>;

public:
    //! samples: number of samples along z axis
    CylinderIntegration(std::size_t zSamples)
    : zSamples_(zSamples)
    {}

    void setGeometry(const GlobalPosition& a, const GlobalPosition& b, Scalar radius)
    {
        const auto ab = (b-a);
        const auto height = ab.two_norm();
        std::size_t zSamples = 1;

        const auto zStep = height/Scalar(zSamples);
        const auto rSamples = std::max(1, int(std::ceil(radius/height*Scalar(zSamples_))));
        const auto rStep = radius/Scalar(rSamples);
        std::size_t circleSamples = 0;
        std::vector<std::size_t> ks(rSamples);
        for (int j = 0; j < rSamples; ++j)
        {
            // number of cells in ring j
            ks[j] = std::floor(2*M_PI*(Scalar(j)+0.5));
            circleSamples += ks[j];
        }

        // compute offsets for index calculation
        auto kOffset = ks;
        std::partial_sum(kOffset.begin(), kOffset.end(), kOffset.begin());

        // compute total number of samples
        samples_ = zSamples*circleSamples;
        std::cout << "CylinderIntegration: number of sample points: " << samples_ << std::endl;

        integrationElement_.resize(samples_);
        points_.resize(samples_);

        // generate the points on the disc
        Scalar integral = 0.0;
        for (std::size_t i = 0; i < zSamples; ++i)
        {
            // for each cylinder slice i
            auto center = a;
            center.axpy((Scalar(i)+0.5)*zStep, ab);

            // generate circle point for each ring j
            for (std::size_t j = 0; j < rSamples; ++j)
            {
                auto r = (Scalar(j)+0.5)*rStep;
                const auto ringSamples = ks[j];
                auto cPoints = EmbeddedCoupling::circlePoints(center, ab, r, ringSamples);
                const auto ringOffest = j > 0 ? kOffset[j-1] : 0;
                const auto integrationElement = M_PI*rStep*rStep*zStep*(1.0 + 2.0*Scalar(j))/ringSamples;
                for (std::size_t k = 0; k < ringSamples; ++k)
                {
                    std::size_t idx = k + ringOffest + circleSamples*i;
                    points_[idx] = cPoints[k];
                    integrationElement_[idx] = integrationElement;
                    integral += integrationElement_[idx];
                }
            }
        }

        const auto exactIntegral = M_PI*radius*radius*height;
        std::cout << "Integration error in %: " << std::abs(exactIntegral-integral)/exactIntegral*100.0 << std::endl;
    }

    Scalar integrationElement(std::size_t i) const
    { return integrationElement_[i]; }

    GlobalPosition getIntegrationPoint(std::size_t i) const
    { return points_[i]; }

    std::size_t size() const
    { return samples_; }

private:
    std::size_t zSamples_, samples_;
    std::vector<Scalar> integrationElement_;
    std::vector<GlobalPosition> points_;
};

} // end namespace Dumux

#endif
