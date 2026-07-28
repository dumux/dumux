// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KEpsilonMassLocalResidual
 */
#ifndef DUMUX_RANS_KEPSILON_MASS_LOCAL_RESIDUAL_HH
#define DUMUX_RANS_KEPSILON_MASS_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/navierstokes/mass/1p/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Element-wise storage/flux for the single-phase Navier-Stokes mass balance fused with
 *        the two-equation, wall-function k-epsilon turbulence transport equations (k, epsilon).
 *
 * Extends Dumux::NavierStokesMassOnePLocalResidual's storage/flux with the k/epsilon transport
 * terms, ported from the deleted releases/3.10:dumux/freeflow/rans/twoeq/kepsilon/staggered/
 * fluxvariables.hh's computeMassFlux() (advective + diffusive terms, plus the k-diffusive-flux
 * suppression between matching-point/near-wall-region cells - "don't let the algebraic
 * wall-function region diffusively couple back into the resolved region"). See
 * dumux/freeflow/rans/oneeq/localresidual.hh for why boundary "outside"/ghost volVars must be
 * fetched via elemVars[scvf.outsideScvIdx()] directly rather than fvGeometry.scv(...).
 *
 * Production/destruction source terms (including the matching-point production/destruction
 * skip for k) are *not* added here: they go through the standard problem.source(...) extension
 * point instead (Dumux::KEpsilonMassProblem).
 */
template<class TypeTag>
class KEpsilonMassLocalResidual
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
    //! Upwinded advective (rho*k, rho*epsilon) plus two-point-TPFA diffusive
    //! ((mu_t/sigma_k+mu), (mu_t/sigma_epsilon+mu)) terms - ported from
    //! releases/3.10:dumux/freeflow/rans/twoeq/kepsilon/staggered/fluxvariables.hh::computeMassFlux().
    //! The k-diffusive term is additionally suppressed across a face where both sides, or one
    //! side, are matching-point cells adjacent to a near-wall-region cell, exactly as old code
    //! did - preventing the algebraic wall-function region from diffusively coupling back into
    //! the resolved log-layer solution.
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

        const auto velocity = problem.faceVelocity(element, fvGeometry, scvf);
        const auto vn = velocity*scvf.unitOuterNormal();
        const auto area = Extrusion::area(fvGeometry, scvf)*insideVars.extrusionFactor();

        const auto rhoKUpwind = vn > 0.0 ? insideVars.density()*insideVars.turbulentKineticEnergy()
                                         : outsideVars.density()*outsideVars.turbulentKineticEnergy();
        const auto rhoEpsilonUpwind = vn > 0.0 ? insideVars.density()*insideVars.dissipation()
                                               : outsideVars.density()*outsideVars.dissipation();

        flux[Indices::turbulentKineticEnergyEqIdx] += vn*rhoKUpwind*area;
        flux[Indices::dissipationEqIdx] += vn*rhoEpsilonUpwind*area;

        const auto distance = scvf.boundary()
            ? (insideScv.dofPosition() - scvf.ipGlobal()).two_norm()
            : (fvGeometry.scv(scvf.outsideScvIdx()).dofPosition() - insideScv.dofPosition()).two_norm();

        const auto diffCoeffKInside = insideVars.dynamicEddyViscosity()/insideVars.sigmaK() + insideVars.viscosity();
        const auto diffCoeffKOutside = outsideVars.dynamicEddyViscosity()/outsideVars.sigmaK() + outsideVars.viscosity();
        const auto diffCoeffEpsilonInside = insideVars.dynamicEddyViscosity()/insideVars.sigmaEpsilon() + insideVars.viscosity();
        const auto diffCoeffEpsilonOutside = outsideVars.dynamicEddyViscosity()/outsideVars.sigmaEpsilon() + outsideVars.viscosity();

        const auto diffCoeffK = scvf.boundary() ? diffCoeffKInside : 0.5*(diffCoeffKInside + diffCoeffKOutside);
        const auto diffCoeffEpsilon = scvf.boundary() ? diffCoeffEpsilonInside : 0.5*(diffCoeffEpsilonInside + diffCoeffEpsilonOutside);

        // Suppressed between matching-point/near-wall-region cells (avoids diffusively coupling
        // the algebraic wall-function region back into the resolved solution), and unconditionally
        // at a wall boundary face itself (k has no natural diffusive flux to add there - its value
        // is set entirely by the internal Dirichlet constraint instead, see massproblem.hh).
        const bool suppressKDiffusion = (insideVars.isMatchingPoint() && outsideVars.isMatchingPoint())
            || (insideVars.isMatchingPoint() && outsideVars.inNearWallRegion())
            || (insideVars.inNearWallRegion() && outsideVars.isMatchingPoint())
            || (scvf.boundary() && problem.isOnWallAtPos(scvf.center()));

        if (!suppressKDiffusion)
        {
            flux[Indices::turbulentKineticEnergyEqIdx]
                += diffCoeffK*(insideVars.turbulentKineticEnergy() - outsideVars.turbulentKineticEnergy())/distance*area;
        }

        flux[Indices::dissipationEqIdx]
            += diffCoeffEpsilon*(insideVars.dissipation() - outsideVars.dissipation())/distance*area;
    }
};

} // end namespace Dumux

#endif
