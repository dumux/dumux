
// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief CVFE advection type for the streamfunction-vorticity Boussinesq formulation.
 *
 * Primary variable 0 is the streamfunction ψ (stored in the pressure slot).
 * The Darcy velocity is recovered as u = curl(ψ):
 *
 *   u_x =  ∂ψ/∂y,   u_y = -∂ψ/∂x
 *
 * This automatically satisfies ∇·u = 0 and encodes buoyancy through the
 * streamfunction Poisson equation (handled in the local residual).
 */
#ifndef DUMUX_BOUSSINESQ_STREAMFUNCTION_ADVECTION_HH
#define DUMUX_BOUSSINESQ_STREAMFUNCTION_ADVECTION_HH

#include <dune/common/fvector.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \brief CVFE advection type that computes face fluxes from u = curl(ψ).
 *
 * \tparam Scalar      scalar type
 * \tparam GridGeometry grid geometry type (must be 2D)
 */
template<class Scalar, class GridGeometry>
class StreamfunctionAdvection
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dimWorld == 2, "StreamfunctionAdvection is only implemented for 2D");

public:
    /*!
     * \brief Returns u·n·|σ| where u = curl(ψ) and ψ is primary variable 0.
     *
     * The volumetric face flux is u·n·|σ| = (∂ψ/∂y·nₓ - ∂ψ/∂x·nᵧ)·|σ|.
     * The compositional local residual multiplies this by the upwind mass or
     * mole fraction to get component fluxes.
     */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarCache)
    {
        const auto& fluxVarCache = elemFluxVarCache[scvf];

        // ∇ψ = Σ_i ψ_i · ∇N_i   (ψ lives in the pressure slot)
        Dune::FieldVector<Scalar, dimWorld> gradPsi(0.0);
        for (auto&& scv : scvs(fvGeometry))
            gradPsi.axpy(elemVolVars[scv].pressure(phaseIdx),
                         fluxVarCache.gradN(scv.indexInElement()));

        // u = curl(ψ):  u_x = ∂ψ/∂y,  u_y = -∂ψ/∂x
        const Dune::FieldVector<Scalar, dimWorld> u{ gradPsi[1], -gradPsi[0] };

        return (u * scvf.unitOuterNormal()) * Extrusion::area(fvGeometry, scvf);
    }

    /*!
     * \brief Per-node transmissibilities for the analytical Jacobian / VTK velocity output.
     *
     * flux = Σ_i ti[i] · ψ_i,  so  ti[i] = (∂N_i/∂y·nₓ - ∂N_i/∂x·nᵧ)·|σ|.
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
        const auto& n = scvf.unitOuterNormal();
        const Scalar area = Extrusion::area(fvGeometry, scvf);

        std::vector<Scalar> ti(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& gradN = fluxVarCache.gradN(scv.indexInElement());
            // ti = (∂N/∂y · nₓ  -  ∂N/∂x · nᵧ) · area
            ti[scv.indexInElement()] = (gradN[1]*n[0] - gradN[0]*n[1]) * area;
        }
        return ti;
    }
};

} // namespace Dumux

#endif
