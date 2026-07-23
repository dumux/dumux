// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief CVFE Darcy law for the dimensionless Boussinesq approximation.
 *
 * Identical to CVFEDarcysLaw except that the density used in the gravity
 * (buoyancy) term is
 *
 *   ρ_buoy = 1 + C
 *
 * where C = massFraction(phase0, soluteIdx) is the dimensionless solute
 * concentration.  The fluid-system density is kept at 1 (Boussinesq), so
 * the continuity and transport storage terms remain incompressible.
 *
 * With g = (0, -1) (dimensionless unit gravity pointing down) this gives
 *
 *   q = -(K/μ)( ∇p + (1+C)·g ) = -(K/μ)( ∇p - (1+C)·ẑ ) = -(∇P - C·ẑ)
 *
 * where P = p - z is the dynamic pressure, matching eq. (2) of the
 * dimensionless Boussinesq system.
 */
#ifndef DUMUX_FLUX_CVFE_BOUSSINESQ_DARCYS_LAW_HH
#define DUMUX_FLUX_CVFE_BOUSSINESQ_DARCYS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/facetensoraverage.hh>
#include <dumux/flux/cvfe/darcyslaw.hh>

namespace Dumux {

/*!
 * \ingroup CVFEFlux
 * \brief CVFE Darcy law with Boussinesq buoyancy.
 *
 * The buoyancy density is assembled as ρ·(1 + problem.boussinesq_term(volVars)),
 * i.e. the fluid-system density scaled by one plus a problem-supplied dimensionless
 * driver (typically the transported solute concentration) rather than a density that
 * itself varies with composition -- the defining simplification of the Boussinesq
 * approximation.
 *
 * \tparam Scalar      scalar type
 * \tparam GridGeometry grid geometry type
 */
template<class Scalar, class GridGeometry>
class BoussinesqCVFEDarcyLaw : public CVFEDarcysLaw<Scalar, GridGeometry>
{
    using ParentType = CVFEDarcysLaw<Scalar, GridGeometry>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView  = typename GridGeometry::GridView;
    using Element   = typename GridView::template Codim<0>::Entity;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    // transmissibilities are gravity-independent, so the base class's implementation
    // applies unchanged
    using ParentType::calculateTransmissibilities;

    /*!
     * \brief Advective (Darcy) flux with Boussinesq buoyancy.
     *
     * The buoyancy density is (1 + C) rather than the fluid-system density
     * (which equals 1 for Boussinesq).
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

        const auto& insideVolVars  = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
        const auto& outsideVolVars = elemVolVars[fvGeometry.scv(scvf.outsideScvIdx())];

        auto insideK  = insideVolVars.permeability();
        auto outsideK = outsideVolVars.permeability();
        insideK  *= insideVolVars.extrusionFactor();
        outsideK *= outsideVolVars.extrusionFactor();
        const auto K = faceTensorAverage(insideK, outsideK, scvf.unitOuterNormal());

        static const bool enableGravity =
            getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        const auto& shapeValues = fluxVarCache.shapeValues();

        Dune::FieldVector<Scalar, dimWorld> gradP(0.0);
        Scalar rho(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            const Scalar N = shapeValues[scv.indexInElement()][0];

            if (enableGravity)
                rho += N * volVars.density(phaseIdx) * (1.0 + problem.boussinesq_term(volVars));

            gradP.axpy(volVars.pressure(phaseIdx), fluxVarCache.gradN(scv.indexInElement()));
        }

        if (enableGravity)
            gradP.axpy(-rho, problem.spatialParams().gravity(scvf.center()));

        return -1.0 * vtmv(scvf.unitOuterNormal(), K, gradP)
               * Extrusion::area(fvGeometry, scvf);
    }
};

} // namespace Dumux

#endif