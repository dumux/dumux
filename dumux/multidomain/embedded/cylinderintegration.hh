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

    const GlobalPosition& getIntegrationPoint(std::size_t i) const
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
    CylinderIntegration(std::size_t rSamples, std::size_t zSamples = 0)
    : zSamples_(zSamples), rSamples_(rSamples), circleSamples_(0)
    {
        // precompute sin/cos for the circle points
        ks_.resize(rSamples);
        for (int j = 0; j < rSamples; ++j)
        {
            // number of cells in ring j
            ks_[j] = std::floor(2*M_PI*(Scalar(j)+0.5));
            circleSamples_ += ks_[j];
        }

        // optimization because calling too many sin/cos can be really expensive
        sincos_.resize(rSamples);
        for (int j = 0; j < rSamples; ++j)
        {
            const auto numPoints = ks_[j];
            sincos_[j].resize(2*numPoints);
            Scalar t = 0 + 0.1;
            for (std::size_t i = 0; i < numPoints; ++i)
            {
                sincos_[j][2*i] = sin(t);
                sincos_[j][2*i + 1] = cos(t);
                t += 2*M_PI/numPoints;
                if(t > 2*M_PI) t -= 2*M_PI;
            }
        }
    }

    void setGeometry(const GlobalPosition& a, const GlobalPosition& b, Scalar radius)
    {
        const auto ab = (b-a);
        const auto height = ab.two_norm();

        // compute number of samples in z and r direction
        const std::size_t rSamples = rSamples_;
        // compute zSamples from r samples if not specified by the user
        const std::size_t zSamples = zSamples_ == 0 ? std::max(1, int(std::ceil(height/radius*Scalar(rSamples)))) : zSamples_;
        const auto zStep = height/Scalar(zSamples);
        const auto rStep = radius/Scalar(rSamples);

        // compute offsets for index calculation
        auto kOffset = ks_;
        std::partial_sum(kOffset.begin(), kOffset.end(), kOffset.begin());

        // compute total number of samples
        samples_ = zSamples*circleSamples_;
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
            std::vector<GlobalPosition> ringPoints(ks_[rSamples-1]);
            for (std::size_t j = 0; j < rSamples; ++j)
            {
                auto r = (Scalar(j)+0.5)*rStep;
                const auto ringSamples = ks_[j];
                EmbeddedCoupling::circlePoints(ringPoints, sincos_[j], center, ab, r, ringSamples);
                const auto ringOffest = j > 0 ? kOffset[j-1] : 0;
                const auto integrationElement = M_PI*rStep*rStep*zStep*(1.0 + 2.0*Scalar(j))/ringSamples;
                for (std::size_t k = 0; k < ringSamples; ++k)
                {
                    std::size_t idx = k + ringOffest + circleSamples_*i;
                    points_[idx] = ringPoints[k];
                    integrationElement_[idx] = integrationElement;
                    integral += integrationElement_[idx];
                }
            }
        }
    }

    Scalar integrationElement(unsigned int i) const
    { return integrationElement_[i]; }

    const GlobalPosition& getIntegrationPoint(unsigned int i) const
    { return points_[i]; }

    std::size_t size() const
    { return samples_; }

private:
    std::size_t zSamples_, rSamples_, samples_, circleSamples_;
    std::vector<Scalar> integrationElement_;
    std::vector<GlobalPosition> points_;
    std::vector<std::size_t> ks_;
    std::vector<std::vector<Scalar>> sincos_;

};

template<class Scalar,
         CylinderIntegrationMethod m = CylinderIntegrationMethod::spaced>
class EllipticCylinderIntegration;

template<class Scalar>
class EllipticCylinderIntegration<Scalar, CylinderIntegrationMethod::spaced>
{
    using GlobalPosition = Dune::FieldVector<Scalar, 3>;

public:
    //! samples: number of samples along z axis
    explicit EllipticCylinderIntegration(Scalar deltaA)
    : deltaA_(deltaA)
    {}

    // elliptic cylinder with centerline points p, q and ellipse axis a, b
    void setGeometry(const GlobalPosition& p, const GlobalPosition& q, GlobalPosition aVec, GlobalPosition bVec)
    {
        auto zAxis = (q-p);
        const auto length = zAxis.two_norm();
        zAxis /= length;
        const auto a = aVec.two_norm();
        const auto b = bVec.two_norm();
        aVec /= a;
        bVec /= b;
        const auto normal = crossProduct(aVec, bVec);
        const auto height = (zAxis*normal)*length;

        const auto aStep = deltaA_;
        std::size_t aSamples = std::max(int(std::floor(2*a/aStep)), 1);
        std::size_t bSamples = std::max(int(std::floor(2*b/aStep)), 1);
        const auto bStep = 2*b / Scalar(bSamples);
        std::size_t zSamples = std::max(int(std::floor(height/aStep)), 1);
        const auto zStep = length / Scalar(zSamples);
        auto zStepFactor = length/height;

        integrationElement_ = aStep*bStep*zStep/zStepFactor;
        points_.clear();
        points_.reserve(aSamples*bSamples*zSamples);

        // walk over lattice grid an reject point outside the ellipse
        auto startZ = p; startZ.axpy(0.5*zStep, zAxis);
        auto startAB = startZ;
        startAB.axpy(-a + 0.5*aStep, aVec);
        startAB.axpy(-b + 0.5*bStep, bVec);
        for (std::size_t as = 0; as < aSamples; ++as)
        {
            for (std::size_t bs = 0; bs < bSamples; ++bs)
            {
                auto pos = startAB;
                pos.axpy(as*aStep, aVec);
                pos.axpy(bs*bStep, bVec);
                if (pointInEllipse(pos, startZ, aVec, bVec, normal, a, b))
                    points_.emplace_back(std::move(pos));
            }
        }

        const auto abPoints = points_.size();
        samples_ = abPoints*zSamples;
        points_.resize(samples_);

        for (std::size_t zs = 1; zs < zSamples; ++zs)
        {
            auto add = zAxis; add *= zs*zStep;
            std::transform(points_.begin(), points_.begin() + abPoints, points_.begin() + zs*abPoints,
                           [&add](const auto& p){ return p + add; });
        }

        const auto meanLocalError = (height*a*b*M_PI - integrationElement_*samples_)/Scalar(samples_);
        integrationElement_ += meanLocalError;

        std::cout << "Elliptic cylinder integration -- samples: " << samples_ <<"\n";
        // std::cout << "  -- samples: " << samples_ << ", a: " << a << ", b: " << b << ", h: " << height
        //           << "  -- a-samples: " << aSamples << ", b-samples: " << bSamples << ", z-samples: " << zSamples
        //           << ", volume: " << integrationElement_*samples_ << " (expected " << height*a*b*M_PI << ")\n";
        // std::cout << "  -- volume error was: " << meanLocalError/integrationElement_*100 << "%\n";
    }

    Scalar integrationElement(unsigned int i) const
    { return integrationElement_; }

    const GlobalPosition& getIntegrationPoint(unsigned int i) const
    { return points_[i]; }

    std::size_t size() const
    { return samples_; }

private:
    bool pointInEllipse(const GlobalPosition& p, const GlobalPosition& center,
                        const GlobalPosition& aVec, const GlobalPosition& bVec,
                        const GlobalPosition& normal,
                        const Scalar a, const Scalar b)
    {
        const auto d = p-center;
        // check if point is in ellipse plane
        if (d*normal > 1e-7*a)
            return false;

        const auto da = (d*aVec);
        const auto db = (d*bVec);

        return (da*da/(a*a) + db*db/(b*b) < 1.0);
    }

    Scalar deltaA_;
    std::size_t samples_;
    Scalar integrationElement_;
    std::vector<GlobalPosition> points_;
};


template<class Scalar>
class EllipseIntegration
{
    using GlobalPosition = Dune::FieldVector<Scalar, 3>;

public:
    EllipseIntegration(Scalar deltaA)
    : deltaA_(deltaA)
    {}

    // ellipse cylinder with center point p and ellipse axes a, b
    void setGeometry(const GlobalPosition& p, GlobalPosition aVec, GlobalPosition bVec)
    {
        const auto a = aVec.two_norm();
        const auto b = bVec.two_norm();
        aVec /= a;
        bVec /= b;
        const auto normal = crossProduct(aVec, bVec);

        const auto aStep = deltaA_;
        std::size_t aSamples = std::max(int(std::floor(2*a/aStep)), 1);
        std::size_t bSamples = std::max(int(std::floor(2*b/aStep)), 1);
        const auto bStep = 2*b / Scalar(bSamples);

        integrationElement_ = aStep*bStep;
        points_.reserve(aSamples*bSamples);

        // walk over lattice grid an reject point outside the ellipse
        auto startZ = p;
        auto startAB = startZ;
        startAB.axpy(-a + 0.5*aStep, aVec);
        startAB.axpy(-b + 0.5*bStep, bVec);
        for (std::size_t as = 0; as < aSamples; ++as)
        {
            for (std::size_t bs = 0; bs < bSamples; ++bs)
            {
                auto pos = startAB;
                pos.axpy(as*aStep, aVec);
                pos.axpy(bs*bStep, bVec);
                if (pointInEllipse(pos, startZ, aVec, bVec, normal, a, b))
                    points_.emplace_back(std::move(pos));
            }
        }

        const auto abPoints = points_.size();
        samples_ = abPoints;

        const auto meanLocalError = (a*b*M_PI - integrationElement_*samples_)/Scalar(samples_);
        integrationElement_ += meanLocalError;

        std::cout << "Ellipse integration:\n";
        std::cout << "  -- samples: " << samples_ << ", a: " << a << ", b: " << b
                  << ", area: " << integrationElement_*samples_ << " (expected " << a*b*M_PI << ")\n";
        std::cout << "  -- area error was: " << meanLocalError/integrationElement_*100 << "%\n";
    }

    Scalar integrationElement(unsigned int i) const
    { return integrationElement_; }

    const GlobalPosition& getIntegrationPoint(unsigned int i) const
    { return points_[i]; }

    std::size_t size() const
    { return samples_; }

private:
    bool pointInEllipse(const GlobalPosition& p, const GlobalPosition& center,
                        const GlobalPosition& aVec, const GlobalPosition& bVec,
                        const GlobalPosition& normal,
                        const Scalar a, const Scalar b)
    {
        const auto d = p-center;
        // check if point is in ellipse plane
        if (d*normal > 1e-7*a)
            return false;

        const auto da = (d*aVec);
        const auto db = (d*bVec);

        return (da*da/(a*a) + db*db/(b*b) < 1.0);
    }

    Scalar deltaA_;
    std::size_t samples_;
    Scalar integrationElement_;
    std::vector<GlobalPosition> points_;
};

} // end namespace Dumux

#endif
