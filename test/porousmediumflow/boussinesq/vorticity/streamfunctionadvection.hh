// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief CVFE advection type for the vector-potential Boussinesq formulation (2D/3D).
 *
 * The Darcy velocity is recovered as u = ∇ × A:
 *
 *   2D (single potential ψ ≡ A_z):
 *     u_x =  ∂ψ/∂y,   u_y = -∂ψ/∂x
 *
 *   3D (three components A_x, A_y, A_z):
 *     u   = ∇ × (A_x, A_y, A_z)  (standard curl)
 *
 * This guarantees ∇·u = 0 identically and encodes buoyancy through the
 * vector-potential equations (handled in the local residual).
 *
 * Note: this class is primarily used by external callers (e.g. compositional
 * parent residuals).  The standalone BoussinesqVorticityLocalResidual computes
 * the curl directly without going through this helper.
 */
#ifndef DUMUX_BOUSSINESQ_VECTOR_POTENTIAL_ADVECTION_HH
#define DUMUX_BOUSSINESQ_VECTOR_POTENTIAL_ADVECTION_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \brief CVFE advection type computing  (∇ × A)·n·|σ|  as the volumetric face flux.
 *
 * \tparam Scalar       floating-point type
 * \tparam GridGeometry grid geometry type (2D or 3D)
 */
template<class Scalar, class GridGeometry>
class VectorPotentialAdvection
{
    using FVElementGeometry    = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Extrusion            = Extrusion_t<GridGeometry>;
    using GridView             = typename GridGeometry::GridView;
    using Element              = typename GridView::template Codim<0>::Entity;
    using GlobalPosition       = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int nPot     = (dimWorld == 2) ? 1 : dimWorld;

public:
    /*!
     * \brief Returns (∇ × A)·n·|σ|.
     *
     * In 2D: ψ lives in vectorPotential(0); curl gives (∂ψ/∂y, -∂ψ/∂x).
     * In 3D: full curl of the three potential components.
     */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem&,
                       const Element&,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int /*phaseIdx*/,
                       const ElementFluxVarsCache& elemFluxVarCache)
    {
        const auto& fluxVarCache = elemFluxVarCache[scvf];

        std::array<GlobalPosition, nPot> gradA;
        for (auto& g : gradA) g = GlobalPosition(0.0);

        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& gradN = fluxVarCache.gradN(scv.indexInElement());
            for (int k = 0; k < nPot; ++k)
                gradA[k].axpy(elemVolVars[scv].vectorPotential(k), gradN);
        }

        const auto& n    = scvf.unitOuterNormal();
        const Scalar area = Extrusion::area(fvGeometry, scvf);

        if constexpr (dimWorld == 2)
            return (gradA[0][1]*n[0] - gradA[0][0]*n[1]) * area;
        else
        {
            const Scalar ux = gradA[2][1] - gradA[1][2];
            const Scalar uy = gradA[0][2] - gradA[2][0];
            const Scalar uz = gradA[1][0] - gradA[0][1];
            return (ux*n[0] + uy*n[1] + uz*n[2]) * area;
        }
    }

    /*!
     * \brief Per-node transmissibilities for the 2D analytical Jacobian.
     *
     * flux = Σ_i ti[i] · ψ_i,  so  ti[i] = (∂N_i/∂y·n_x - ∂N_i/∂x·n_y)·|σ|.
     *
     * For 3D the transmissibilities depend on all three potential components and
     * a single vector is insufficient; use numerical differentiation instead.
     */
    template<class Problem, class ElementVolumeVariables, class FluxVarCache>
    static std::vector<Scalar>
    calculateTransmissibilities(const Problem&,
                                const Element&,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables&,
                                const SubControlVolumeFace& scvf,
                                const FluxVarCache& fluxVarCache)
    {
        static_assert(dimWorld == 2,
            "calculateTransmissibilities is only implemented for 2D. "
            "Use numerical differentiation for 3D.");

        const auto& n    = scvf.unitOuterNormal();
        const Scalar area = Extrusion::area(fvGeometry, scvf);

        std::vector<Scalar> ti(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& gradN = fluxVarCache.gradN(scv.indexInElement());
            ti[scv.indexInElement()] = (gradN[1]*n[0] - gradN[0]*n[1]) * area;
        }
        return ti;
    }
};

// Backward-compat alias
template<class Scalar, class GridGeometry>
using StreamfunctionAdvection = VectorPotentialAdvection<Scalar, GridGeometry>;

} // namespace Dumux

#endif
