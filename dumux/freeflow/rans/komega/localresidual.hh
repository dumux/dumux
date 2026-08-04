// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KOmegaMassLocalResidual
 */
#ifndef DUMUX_RANS_KOMEGA_MASS_LOCAL_RESIDUAL_HH
#define DUMUX_RANS_KOMEGA_MASS_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/navierstokes/mass/1p/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Element-wise storage/flux for the single-phase Navier-Stokes mass balance fused with
 *        the two-equation k-omega turbulence transport equations (k, omega).
 *
 * Extends Dumux::NavierStokesMassOnePLocalResidual's storage/flux with the k/omega transport
 * terms, ported from the deleted releases/3.10:dumux/freeflow/rans/twoeq/komega/staggered/
 * fluxvariables.hh's computeMassFlux() (advective + diffusive terms), mirroring
 * dumux/freeflow/rans/oneeq/localresidual.hh's structure exactly (two equations instead of
 * one). See that file's documentation for why boundary "outside"/ghost volVars must be
 * fetched via elemVars[scvf.outsideScvIdx()] directly rather than fvGeometry.scv(...) - the
 * exact bug that caused the one-equation port's singular Jacobian, repeated here to avoid it.
 *
 * Production/destruction/cross-diffusion source terms are *not* added here: they go through
 * the standard problem.source(...) extension point instead (Dumux::KOmegaMassProblem).
 */
template<class TypeTag>
class KOmegaMassLocalResidual
: public NavierStokesMassOnePLocalResidual<TypeTag>
{
    using ParentType = NavierStokesMassOnePLocalResidual<TypeTag>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = Concept::GridVariablesCache_t<GridVariables>;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using Variables = Concept::Variables_t<GridVariables>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    using Extrusion = Extrusion_t<GridGeometry>;

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const Variables& vars) const
    {
        NumEqVector storage = ParentType::computeStorage(problem, scv, vars);
        storage[Indices::turbulentKineticEnergyEqIdx] = vars.density()*vars.turbulentKineticEnergy();
        storage[Indices::dissipationEqIdx] = vars.density()*vars.dissipation();
        return storage;
    }

    NumEqVector storageIntegral(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        NumEqVector storage = ParentType::storageIntegral(fvGeometry, elemVars, scv, isPreviousTimeLevel);
        const auto& vars = elemVars[scv];
        const auto volumeFactor = Extrusion::volume(fvGeometry, scv)*vars.extrusionFactor();
        storage[Indices::turbulentKineticEnergyEqIdx] = vars.density()*vars.turbulentKineticEnergy()*volumeFactor;
        storage[Indices::dissipationEqIdx] = vars.density()*vars.dissipation()*volumeFactor;
        return storage;
    }

    template<class ElementFluxVariablesCache>
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVariables& elemVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        NumEqVector flux = ParentType::computeFlux(problem, element, fvGeometry, elemVars, scvf, elemFluxVarsCache);
        addTurbulenceFlux_(flux, problem, element, fvGeometry, elemVars, scvf);
        return flux;
    }

    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux = ParentType::fluxIntegral(fvGeometry, elemVars, scvf);
        addTurbulenceFlux_(flux, this->asImp().problem(), fvGeometry.element(), fvGeometry, elemVars, scvf);
        return flux;
    }

private:
    //! Two-point-TPFA diffusive ((mu+sigma_k*mu_t), (mu+sigma_omega*mu_t)) terms - ported from
    //! releases/3.10:dumux/freeflow/rans/twoeq/komega/staggered/fluxvariables.hh::computeMassFlux().
    //! See class docs above re: boundary ghost-volVars handling.
    //!
    //! \note The advective (upwind) part of the k/omega flux is *not* added here - it is already
    //! included in ParentType::computeFlux()'s result via the KOmegaMassModelTraits-specific
    //! Dumux::AdvectiveFlux specialization (dumux/freeflow/rans/komega/advectiveflux.hh), which
    //! NavierStokesMassOnePFluxVariables::advectiveFlux() dispatches to automatically for every
    //! equation index, not just conti0EqIdx. Adding it again here would double-count k/omega
    //! advection (a bug found and fixed while implementing Phases 6/7).
    void addTurbulenceFlux_(NumEqVector& flux,
                            const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVariables& elemVars,
                            const SubControlVolumeFace& scvf) const
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVars = elemVars[insideScv];
        const auto& outsideVars = elemVars[scvf.outsideScvIdx()];

        const auto area = Extrusion::area(fvGeometry, scvf)*insideVars.extrusionFactor();

        const auto diffCoeffKInside = insideVars.viscosity() + insideVars.sigmaK()*insideVars.density()*insideVars.turbulentKineticEnergy()/insideVars.dissipation();
        const auto diffCoeffKOutside = outsideVars.viscosity() + outsideVars.sigmaK()*outsideVars.density()*outsideVars.turbulentKineticEnergy()/outsideVars.dissipation();
        const auto diffCoeffOmegaInside = insideVars.viscosity() + insideVars.sigmaOmega()*insideVars.density()*insideVars.turbulentKineticEnergy()/insideVars.dissipation();
        const auto diffCoeffOmegaOutside = outsideVars.viscosity() + outsideVars.sigmaOmega()*outsideVars.density()*outsideVars.turbulentKineticEnergy()/outsideVars.dissipation();

        // Distance-weighted (not plain-averaged) interpolation of the diffusion coefficient to
        // the face, and the corresponding cell-to-cell distance - reproducing
        // releases/3.10:dumux/freeflow/rans/twoeq/komega/staggered/fluxvariables.hh's
        // arithmeticMean()-based interior-face treatment verbatim. On a graded mesh, a plain
        // 0.5*(inside+outside) average (used here previously) deviates most exactly where the
        // diffusion coefficient itself jumps most sharply between neighboring cells - found to
        // matter right at the wall-adjacent cell, where the analytically-pinned omega there
        // implies a much larger eddy viscosity than its neighbor.
        Scalar distance, diffCoeffK, diffCoeffOmega;
        if (scvf.boundary())
        {
            distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            diffCoeffK = diffCoeffKInside;
            diffCoeffOmega = diffCoeffOmegaInside;
        }
        else
        {
            const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
            const auto distanceInside = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            const auto distanceOutside = (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            distance = (outsideScv.dofPosition() - insideScv.dofPosition()).two_norm();
            diffCoeffK = arithmeticMean(diffCoeffKInside, diffCoeffKOutside, distanceOutside, distanceInside);
            diffCoeffOmega = arithmeticMean(diffCoeffOmegaInside, diffCoeffOmegaOutside, distanceOutside, distanceInside);
        }

        flux[Indices::turbulentKineticEnergyEqIdx]
            += diffCoeffK*(insideVars.turbulentKineticEnergy() - outsideVars.turbulentKineticEnergy())/distance*area;
        flux[Indices::dissipationEqIdx]
            += diffCoeffOmega*(insideVars.dissipation() - outsideVars.dissipation())/distance*area;
    }
};

} // end namespace Dumux

#endif
