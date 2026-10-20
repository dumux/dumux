// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterFlux
 * \brief Limited linear reconstruction of the shallow water state at the faces of a cell
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_RECONSTRUCTION_HH
#define DUMUX_FLUX_SHALLOW_WATER_RECONSTRUCTION_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>

namespace Dumux::ShallowWater {

/*!
 * \ingroup ShallowWaterFlux
 * \brief The state of a cell at one of its faces
 */
template<class Scalar>
struct FaceState
{
    Scalar waterDepth;
    Scalar velocityX;
    Scalar velocityY;
    Scalar bedSurface;
};

/*!
 * \ingroup ShallowWaterFlux
 * \brief Limited average of two estimates a and b of the same increment (van Albada average)
 *
 * Returns \f$ \frac{a + b}{2} \frac{2 p(ab) + 2\epsilon^2}{a^2 + b^2 + 2\epsilon^2} \f$ with a continuously
 * differentiable positive part \f$ p(x) = x \f$ for \f$ x \geq 0 \f$ and
 * \f$ p(x) = \epsilon^2 x / (\epsilon^2 - x) \f$ for \f$ x < 0 \f$. The average equals \f$ b \f$ for
 * \f$ a = b \f$ and vanishes for estimates of opposite signs that are large compared to \f$ \epsilon \f$,
 * so the increment does not change sign. Increments small compared to \f$ \epsilon \f$ are averaged
 * without limiting, which makes the average smooth where the field is nearly constant.
 * For \f$ \epsilon = 0 \f$ it is the clipped van Albada average, between zero and
 * \f$ 2\min(|a|, |b|) \f$ for estimates of the same sign.
 *
 * \param a first estimate
 * \param b second estimate
 * \param eps2 the square of the smoothing scale \f$ \epsilon \f$
 */
template<class Scalar>
Scalar vanAlbadaAverage(const Scalar a, const Scalar b, const Scalar eps2)
{
    const Scalar ab = a*b;
    const Scalar positivePart = ab >= 0.0 ? ab : eps2*ab/(eps2 - ab);
    const Scalar numerator = 2.0*(positivePart + eps2);

    // a^2 + b^2 + 2 eps^2 written such that it equals the numerator exactly for a = b
    const Scalar difference = a - b;
    const Scalar denominator = numerator + difference*difference + 2.0*(ab - positivePart);
    if (denominator == 0.0)
        return 0.0;
    return 0.5*(a + b)*(numerator/denominator);
}

/*!
 * \ingroup ShallowWaterFlux
 * \brief Parameters of the linear reconstruction (see CellReconstruction)
 */
template<class Scalar>
struct ReconstructionParameters
{
    bool limit;
    Scalar relativeEpsilon;
    Scalar lowerWaterDepth;
    Scalar upperWaterDepth;
};

/*!
 * \ingroup ShallowWaterFlux
 * \brief Read the parameters of the linear reconstruction
 */
template<class Scalar>
ReconstructionParameters<Scalar> readReconstructionParameters(const std::string& paramGroup = "")
{
    ReconstructionParameters<Scalar> params;
    const auto limiter = getParamFromGroup<std::string>(paramGroup, "ShallowWater.Limiter", "VanAlbada");
    if (limiter == "VanAlbada")
        params.limit = true;
    else if (limiter == "None")
        params.limit = false;
    else
        DUNE_THROW(Dune::InvalidStateException, "Unknown ShallowWater.Limiter " << limiter << " (use VanAlbada or None)");

    params.relativeEpsilon = getParamFromGroup<Scalar>(paramGroup, "ShallowWater.LimiterRelativeEpsilon", 1e-3);

    params.lowerWaterDepth = getParamFromGroup<Scalar>(paramGroup, "ShallowWater.ReconstructionLowerWaterDepth", 1e-5);
    params.upperWaterDepth = getParamFromGroup<Scalar>(paramGroup, "ShallowWater.ReconstructionUpperWaterDepth", 1e-3);
    if (!(params.upperWaterDepth > params.lowerWaterDepth))
        DUNE_THROW(Dune::InvalidStateException, "ShallowWater.ReconstructionUpperWaterDepth has to be larger than "
                                                << "ShallowWater.ReconstructionLowerWaterDepth");
    return params;
}

namespace Detail {

//! Smooth step from zero at the lower to one at the upper water depth, continuously differentiable
template<class Scalar>
Scalar wetWeight(const Scalar waterDepth, const ReconstructionParameters<Scalar>& params)
{
    using std::clamp;
    const Scalar t = clamp((waterDepth - params.lowerWaterDepth)/(params.upperWaterDepth - params.lowerWaterDepth),
                           Scalar(0.0), Scalar(1.0));
    return t*t*(3.0 - 2.0*t);
}

} // end namespace Detail

/*!
 * \ingroup ShallowWaterFlux
 * \brief Limited linear reconstruction of water depth, free surface elevation and velocity in a cell
 *
 * The gradients are least-squares fits over the face neighbors of the cell (weighted with the inverse
 * squared distance), exact for linear fields. In directions without neighbor information (e.g. across
 * a channel one cell wide) the gradient component vanishes.
 *
 * At a face with center \f$ x_f \f$ and neighbor \f$ j \f$, the increment from the cell center predicted
 * by the cell gradient is \f$ b = \nabla\phi_i \cdot (x_f - x_i) \f$, and the neighbor value deviates from
 * the value predicted by the cell gradient by \f$ e = (\phi_j - \phi_i) - \nabla\phi_i \cdot (x_j - x_i) \f$.
 * The face increment is vanAlbadaAverage(b - e/2, b + e/2), or \f$ b \f$ without limiting (parameter
 * `ShallowWater.Limiter`: `VanAlbada` (default) or `None`). Both estimates coincide with \f$ b \f$ for linear
 * fields; on a uniform Cartesian grid they are \f$ \Delta^-/2 \f$ and \f$ \Delta^+/2 \f$ with the backward and
 * forward differences \f$ \Delta^- \f$ and \f$ \Delta^+ \f$, so the face increment is that of the
 * one-dimensional van Albada limiter. It vanishes at extrema and tends to \f$ b \f$ for increments small
 * compared to the smoothing scale of the limiter.
 * At boundary faces, water depth, free surface elevation and the velocity component normal to the face are
 * extrapolated with the increment \f$ b \f$, and the tangential velocity is that of the cell: boundary
 * conditions take incoming quantities such as the tangential velocity at an inflow from the inside state,
 * and extrapolating these can make the scheme unstable. The velocity is limited in its components normal and
 * tangential to the face, so the reconstruction does not depend on the coordinate system. The smoothing scale
 * of the limiter is the parameter `ShallowWater.LimiterRelativeEpsilon` (default 1e-3) times the water depth
 * of the cell for water depth and free surface elevation, and times \f$ \sqrt{|\mathbf{u}|^2 + g h} \f$ for
 * the velocity.
 *
 * Reconstructed are the water depth \f$ h \f$, the free surface elevation \f$ \eta = h + z \f$ and the
 * velocity. The bed at the face follows as \f$ z_f = \eta_f - h_f \f$, so a lake at rest is
 * reconstructed exactly (Audusse et al. 2004, \cite Audusse2004).
 *
 * Near dry cells all increments are multiplied by a smooth weight, the product over the cell and its face
 * neighbors of a step function of the water depth that is zero below
 * `ShallowWater.ReconstructionLowerWaterDepth` (default 1e-5 m) and one above
 * `ShallowWater.ReconstructionUpperWaterDepth` (default 1e-3 m). Negative reconstructed depths are set to zero.
 *
 * \note Requires caching of the grid geometry and the grid volume variables, since the cell and its
 *       neighbors need not be in the stencil bound by the element-local views.
 */
template<class Scalar, class GlobalPosition>
class CellReconstruction
{
    static_assert(GlobalPosition::dimension == 2, "The shallow water reconstruction is implemented for two dimensions");
    using Matrix = Dune::FieldMatrix<Scalar, 2, 2>;

public:
    template<class FVElementGeometry, class ElementVolumeVariables>
    CellReconstruction(const ReconstructionParameters<Scalar>& params,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const std::size_t cellIdx,
                       const Scalar gravity)
    {
        static_assert(ElementVolumeVariables::GridVolumeVariables::cachingEnabled,
                      "The shallow water reconstruction requires caching of the grid volume variables");

        limit_ = params.limit;

        const auto& gridGeometry = fvGeometry.gridGeometry();
        center_ = fvGeometry.scv(cellIdx).center();
        const auto& volVars = elemVolVars[cellIdx];
        waterDepth_ = volVars.waterDepth();
        bedSurface_ = volVars.bedSurface();
        velocity_ = GlobalPosition{volVars.velocity(0), volVars.velocity(1)};
        weight_ = Detail::wetWeight(waterDepth_, params);

        using std::max;
        const Scalar kappa2 = params.relativeEpsilon*params.relativeEpsilon;
        const Scalar depth = max(waterDepth_, Scalar(0.0));
        epsDepth2_ = kappa2*depth*depth;
        epsVelocity2_ = kappa2*(velocity_.two_norm2() + gravity*depth);

        Matrix normalMatrix(0.0);
        GlobalPosition rhsDepth(0.0), rhsSurface(0.0);
        Matrix rhsVelocity(0.0);
        for (const auto scvfIdx : gridGeometry.scvfIndicesOfScv(cellIdx))
        {
            const auto& scvf = fvGeometry.scvf(scvfIdx);
            if (scvf.boundary())
                continue;

            const auto& neighborVolVars = elemVolVars[scvf.outsideScvIdx()];
            auto d = fvGeometry.scv(scvf.outsideScvIdx()).center();
            d -= center_;
            const Scalar w = 1.0/d.two_norm2();

            for (int k = 0; k < 2; ++k)
                for (int l = 0; l < 2; ++l)
                    normalMatrix[k][l] += w*d[k]*d[l];

            const Scalar deltaDepth = neighborVolVars.waterDepth() - waterDepth_;
            const Scalar deltaSurface = deltaDepth + neighborVolVars.bedSurface() - bedSurface_;
            rhsDepth.axpy(w*deltaDepth, d);
            rhsSurface.axpy(w*deltaSurface, d);
            for (int k = 0; k < 2; ++k)
                rhsVelocity[k].axpy(w*(neighborVolVars.velocity(k) - velocity_[k]), d);

            weight_ *= Detail::wetWeight(neighborVolVars.waterDepth(), params);
        }

        const auto pinv = pseudoInverse_(normalMatrix);
        pinv.mv(rhsDepth, gradDepth_);
        pinv.mv(rhsSurface, gradSurface_);
        for (int k = 0; k < 2; ++k)
            pinv.mv(rhsVelocity[k], gradVelocity_[k]);
    }

    //! The state at a face with a neighbor whose volume variables and cell center are given
    template<class VolumeVariables>
    FaceState<Scalar> faceState(const GlobalPosition& faceCenter,
                                const GlobalPosition& unitNormal,
                                const VolumeVariables& neighborVolVars,
                                const GlobalPosition& neighborCenter) const
    {
        auto r = faceCenter; r -= center_;
        auto d = neighborCenter; d -= center_;

        auto increment = [&](const Scalar delta, const GlobalPosition& grad)
        {
            const Scalar b = grad*r;
            if (!limit_)
                return b;
            const Scalar halfDeviation = 0.5*(delta - grad*d);
            return vanAlbadaAverage(b - halfDeviation, b + halfDeviation, epsDepth2_);
        };

        const Scalar deltaDepth = neighborVolVars.waterDepth() - waterDepth_;
        const Scalar deltaSurface = deltaDepth + neighborVolVars.bedSurface() - bedSurface_;

        GlobalPosition b, halfDeviation;
        for (int k = 0; k < 2; ++k)
        {
            b[k] = gradVelocity_[k]*r;
            halfDeviation[k] = 0.5*((neighborVolVars.velocity(k) - velocity_[k]) - gradVelocity_[k]*d);
        }

        GlobalPosition velocityIncrement = b;
        if (limit_)
        {
            const GlobalPosition tangent{-unitNormal[1], unitNormal[0]};
            velocityIncrement = unitNormal;
            const Scalar bn = b*unitNormal, en = halfDeviation*unitNormal;
            const Scalar bt = b*tangent, et = halfDeviation*tangent;
            velocityIncrement *= vanAlbadaAverage(bn - en, bn + en, epsVelocity2_);
            velocityIncrement.axpy(vanAlbadaAverage(bt - et, bt + et, epsVelocity2_), tangent);
        }

        return makeFaceState_(increment(deltaDepth, gradDepth_), increment(deltaSurface, gradSurface_), velocityIncrement);
    }

    //! The state at a boundary face with the given unit outer normal
    FaceState<Scalar> faceState(const GlobalPosition& faceCenter, const GlobalPosition& unitNormal) const
    {
        auto r = faceCenter; r -= center_;
        GlobalPosition b;
        for (int k = 0; k < 2; ++k)
            b[k] = gradVelocity_[k]*r;

        auto velocityIncrement = unitNormal;
        velocityIncrement *= b*unitNormal;
        return makeFaceState_(gradDepth_*r, gradSurface_*r, velocityIncrement);
    }

    Scalar waterDepth() const { return waterDepth_; }
    Scalar bedSurface() const { return bedSurface_; }

private:
    FaceState<Scalar> makeFaceState_(const Scalar depthIncrement,
                                     const Scalar surfaceIncrement,
                                     const GlobalPosition& velocityIncrement) const
    {
        FaceState<Scalar> state;
        state.waterDepth = waterDepth_ + weight_*depthIncrement;
        if (state.waterDepth >= 0.0)
            state.bedSurface = bedSurface_ + weight_*(surfaceIncrement - depthIncrement);
        else
        {
            state.waterDepth = 0.0;
            state.bedSurface = waterDepth_ + bedSurface_ + weight_*surfaceIncrement;
        }
        state.velocityX = velocity_[0] + weight_*velocityIncrement[0];
        state.velocityY = velocity_[1] + weight_*velocityIncrement[1];
        return state;
    }

    // Moore-Penrose inverse of the symmetric positive semi-definite normal matrix
    static Matrix pseudoInverse_(const Matrix& m)
    {
        const Scalar trace = m[0][0] + m[1][1];
        const Scalar det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
        if (!(trace > 0.0))
            return Matrix(0.0);
        if (det > 1e-10*trace*trace)
        {
            auto inv = m;
            inv.invert();
            return inv;
        }

        // rank one: m = trace v v^T with a unit vector v
        auto pinv = m;
        pinv /= trace*trace;
        return pinv;
    }

    bool limit_;
    GlobalPosition center_;
    Scalar waterDepth_;
    Scalar bedSurface_;
    GlobalPosition velocity_;
    Scalar weight_;
    Scalar epsDepth2_;
    Scalar epsVelocity2_;
    GlobalPosition gradDepth_;
    GlobalPosition gradSurface_;
    std::array<GlobalPosition, 2> gradVelocity_;
};

} // end namespace Dumux::ShallowWater

#endif
