// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief CVFE Darcy law for the Boussinesq approximation.
 *
 * Identical to CVFEDarcysLaw except that the density used in the gravity
 * (buoyancy) term is
 *
 *   ρ_buoy = ρ + Δρ
 *
 * where ρ = volVars.density(phaseIdx) is the (Boussinesq: constant reference)
 * fluid-system density, and Δρ = problem.buoyantDensity(volVars) is a
 * problem-supplied density deviation, in the same units as ρ (kg/m^3), that
 * drives buoyancy.
 *
 * \note Interface contract for `Problem::buoyantDensity(volVars)`:
 * - It must return a *real* density deviation Δρ, in kg/m^3 -- not a
 *   dimensionless or otherwise normalized quantity. There is no implicit
 *   scaling applied by this class.
 * - Δρ must be defined consistently with the reference density ρ returned
 *   by the fluid system in use (`volVars.density(phaseIdx)`), since the two
 *   are summed directly in the gravity term below. Typically ρ is a constant
 *   reference density (e.g. that of the solvent) and Δρ is a linear function
 *   of whatever drives the density change (e.g. solute mass fraction or
 *   temperature), such as Δρ = ρ_ref * β * C for a solute with a linear
 *   density law of slope β.
 * - Everywhere else (storage, transport, mobility), the fluid-system density
 *   ρ is used unmodified and is expected to stay at its (Boussinesq)
 *   reference value -- this class does not touch those terms.
 *
 * For a fully worked example (a fluid system with a linear density law, used
 * both in a "full" compressible-density model and in a Boussinesq-approximated
 * one sharing the same problem), see
 * `test/porousmediumflow/1pnc/1p2c/isothermal/boussinesqintrusion/`.
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
 * The buoyancy density is assembled as ρ + problem.buoyantDensity(volVars), i.e.
 * the (Boussinesq: constant reference) fluid-system density plus a problem-supplied
 * real density deviation Δρ [kg/m^3], rather than a density that itself varies with
 * composition -- the defining simplification of the Boussinesq approximation.
 * See the file-level documentation above for the full interface contract of
 * `buoyantDensity()`.
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
     * The buoyancy density is ρ + Δρ, where ρ is the fluid-system density
     * (expected to stay at its Boussinesq reference value) and
     * Δρ = problem.buoyantDensity(volVars) is the problem-supplied density
     * deviation that drives buoyancy.
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
                rho += N * (volVars.density(phaseIdx) + problem.buoyantDensity(volVars));

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
