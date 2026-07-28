// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::OneEqMassLocalResidual
 */
#ifndef DUMUX_RANS_ONEEQ_MASS_LOCAL_RESIDUAL_HH
#define DUMUX_RANS_ONEEQ_MASS_LOCAL_RESIDUAL_HH

#include <cmath>
#include <iostream>

#include <dumux/common/properties.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/navierstokes/mass/1p/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Element-wise storage/flux for the single-phase Navier-Stokes mass balance fused
 *        with the one-equation (Spalart-Allmaras) turbulence transport equation for ν̃.
 *
 * Extends Dumux::NavierStokesMassOnePLocalResidual's storage/flux with the ν̃ transport
 * terms, ported from the deleted releases/3.10:dumux/freeflow/rans/oneeq/staggered/
 * fluxvariables.hh's computeMassFlux() (advective + diffusive terms) - reusing this mass
 * model's own, already-correct, always-live velocity lookup (problem.velocity(...), the
 * same one the pressure equation's own advective flux already uses, see
 * whatisimplemented.md for why this avoids the stale-cached-velocity bug that hit the
 * previous (parked) 3-domain attempt).
 *
 * Production/destruction/cross-diffusion source terms are *not* added here: they go through
 * the standard DuMux problem.source(...) extension point instead (implemented in
 * Dumux::RANSMassOneEqProblem, dumux/freeflow/rans/oneeq/massproblem.hh), since
 * Dumux::DiscretizationDefaultLocalOperator's computeSource() already forwards there
 * generically - no need to override source handling here.
 */
template<class TypeTag>
class OneEqMassLocalResidual
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
        storage[Indices::viscosityTildeEqIdx] = vars.density()*vars.viscosityTilde();
        return storage;
    }

    NumEqVector storageIntegral(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        NumEqVector storage = ParentType::storageIntegral(fvGeometry, elemVars, scv, isPreviousTimeLevel);
        const auto& vars = elemVars[scv];
        storage[Indices::viscosityTildeEqIdx] = vars.density()*vars.viscosityTilde()
                                                * Extrusion::volume(fvGeometry, scv) * vars.extrusionFactor();
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
        flux[Indices::viscosityTildeEqIdx] = viscosityTildeFlux_(problem, element, fvGeometry, elemVars, scvf);
        return flux;
    }

    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux = ParentType::fluxIntegral(fvGeometry, elemVars, scvf);
        flux[Indices::viscosityTildeEqIdx] = viscosityTildeFlux_(this->asImp().problem(), fvGeometry.element(), fvGeometry, elemVars, scvf);
        return flux;
    }

private:
    //! Upwinded advective ρν̃(v·n) plus a two-point-TPFA diffusive ((μ+ρν̃)/σ) term - ported
    //! from releases/3.10:dumux/freeflow/rans/oneeq/staggered/fluxvariables.hh::computeMassFlux(),
    //! reusing problem.faceVelocity(element,fvGeometry,scvf) - the exact same, always-live
    //! velocity lookup dumux/freeflow/navierstokes/scalarfluxvariables.hh's getAdvectiveFlux()
    //! already uses for the pressure equation's own advective flux (see
    //! whatisimplemented.md/proposedimplementation.md).
    //!
    //! \note Called for interior faces *and* for boundary faces with a pure-Dirichlet BC (the
    //!       inlet here, see CCLocalResidual::evalFlux) - Neumann boundaries (walls/outlet) go
    //!       through Problem::neumann() instead. For a boundary face, scvf.outsideScvIdx() does
    //!       not index a real scv in fvGeometry (CCTpfa has no "ghost" scv objects) - the
    //!       corresponding volVars (populated from problem.dirichlet(...) by the framework) must
    //!       be looked up directly from elemVars by that raw index instead, and the "outside"
    //!       position is the face itself (one-sided distance), not a dofPosition().
    Scalar viscosityTildeFlux_(const Problem& problem,
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

        const auto rhoNuUpwind = vn > 0.0 ? insideVars.density()*insideVars.viscosityTilde()
                                          : outsideVars.density()*outsideVars.viscosityTilde();

        Scalar flux = vn*rhoNuUpwind*area;

        const auto distance = scvf.boundary()
            ? (insideScv.dofPosition() - scvf.ipGlobal()).two_norm()
            : (fvGeometry.scv(scvf.outsideScvIdx()).dofPosition() - insideScv.dofPosition()).two_norm();
        const auto diffCoeffInside = (insideVars.viscosity() + insideVars.density()*insideVars.viscosityTilde())/insideVars.sigma();
        const auto diffCoeffOutside = (outsideVars.viscosity() + outsideVars.density()*outsideVars.viscosityTilde())/outsideVars.sigma();
        const auto diffCoeff = scvf.boundary() ? diffCoeffInside : 0.5*(diffCoeffInside + diffCoeffOutside);

        flux -= diffCoeff*(outsideVars.viscosityTilde() - insideVars.viscosityTilde())/distance*area;

        return flux;
    }
};

} // end namespace Dumux

#endif
