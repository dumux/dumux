// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::LowReKEpsilonMassLocalResidual
 */
#ifndef DUMUX_RANS_LOWREKEPSILON_MASS_LOCAL_RESIDUAL_HH
#define DUMUX_RANS_LOWREKEPSILON_MASS_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/navierstokes/mass/1p/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Element-wise storage/flux for the single-phase Navier-Stokes mass balance fused with
 *        the two-equation, low-Reynolds-number k-epsilon turbulence transport equations
 *        (k, epsilon-tilde).
 *
 * Extends Dumux::NavierStokesMassOnePLocalResidual's storage/flux with the k/epsilon-tilde
 * diffusive transport terms, ported from the deleted releases/3.10:dumux/freeflow/rans/
 * twoeq/lowrekepsilon/staggered/fluxvariables.hh's computeMassFlux(). See
 * dumux/freeflow/rans/oneeq/localresidual.hh for why boundary "outside"/ghost volVars must be
 * fetched via elemVars[scvf.outsideScvIdx()] directly rather than fvGeometry.scv(...).
 *
 * \note The advective part of the k/epsilon-tilde flux is *not* added here - it is already
 * included in ParentType::computeFlux()'s result via the LowReKEpsilonMassModelTraits-specific
 * Dumux::AdvectiveFlux specialization (dumux/freeflow/rans/lowrekepsilon/advectiveflux.hh),
 * which NavierStokesMassOnePFluxVariables::advectiveFlux() dispatches to automatically for
 * every equation index, not just conti0EqIdx - see komega/localresidual.hh's note for the
 * double-counting bug this avoids.
 *
 * Production/destruction/damping-function source terms are *not* added here: they go through
 * the standard problem.source(...) extension point instead (Dumux::LowReKEpsilonMassProblem).
 */
template<class TypeTag>
class LowReKEpsilonMassLocalResidual
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
        storage[Indices::dissipationEqIdx] = vars.density()*vars.dissipationTilde();
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
        storage[Indices::dissipationEqIdx] = vars.density()*vars.dissipationTilde()*volumeFactor;
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
        addTurbulenceFlux_(flux, fvGeometry, elemVars, scvf);
        return flux;
    }

    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux = ParentType::fluxIntegral(fvGeometry, elemVars, scvf);
        addTurbulenceFlux_(flux, fvGeometry, elemVars, scvf);
        return flux;
    }

private:
    //! Two-point-TPFA diffusive ((mu+sigma_k*mu_t), (mu+sigma_epsilon*mu_t)) terms - ported from
    //! releases/3.10:dumux/freeflow/rans/twoeq/lowrekepsilon/staggered/fluxvariables.hh::
    //! computeMassFlux(). See class docs above re: boundary ghost-volVars handling and why the
    //! advective part is not (re-)added here.
    void addTurbulenceFlux_(NumEqVector& flux,
                            const FVElementGeometry& fvGeometry,
                            const ElementVariables& elemVars,
                            const SubControlVolumeFace& scvf) const
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVars = elemVars[insideScv];
        const auto& outsideVars = elemVars[scvf.outsideScvIdx()];

        const auto area = Extrusion::area(fvGeometry, scvf)*insideVars.extrusionFactor();

        const auto distance = scvf.boundary()
            ? (insideScv.dofPosition() - scvf.ipGlobal()).two_norm()
            : (fvGeometry.scv(scvf.outsideScvIdx()).dofPosition() - insideScv.dofPosition()).two_norm();

        const auto diffCoeffKInside = insideVars.dynamicEddyViscosity()/insideVars.sigmaK() + insideVars.viscosity();
        const auto diffCoeffKOutside = outsideVars.dynamicEddyViscosity()/outsideVars.sigmaK() + outsideVars.viscosity();
        const auto diffCoeffEpsilonInside = insideVars.dynamicEddyViscosity()/insideVars.sigmaEpsilon() + insideVars.viscosity();
        const auto diffCoeffEpsilonOutside = outsideVars.dynamicEddyViscosity()/outsideVars.sigmaEpsilon() + outsideVars.viscosity();

        const auto diffCoeffK = scvf.boundary() ? diffCoeffKInside : 0.5*(diffCoeffKInside + diffCoeffKOutside);
        const auto diffCoeffEpsilon = scvf.boundary() ? diffCoeffEpsilonInside : 0.5*(diffCoeffEpsilonInside + diffCoeffEpsilonOutside);

        flux[Indices::turbulentKineticEnergyEqIdx]
            += diffCoeffK*(insideVars.turbulentKineticEnergy() - outsideVars.turbulentKineticEnergy())/distance*area;
        flux[Indices::dissipationEqIdx]
            += diffCoeffEpsilon*(insideVars.dissipationTilde() - outsideVars.dissipationTilde())/distance*area;
    }
};

} // end namespace Dumux

#endif
