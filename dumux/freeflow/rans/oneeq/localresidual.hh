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

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/navierstokes/mass/1p/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Element-wise storage/flux for the single-phase Navier-Stokes mass balance
 *        together with the one-equation (Spalart-Allmaras) turbulence transport equation for ν̃.
 *
 * Extends Dumux::NavierStokesMassOnePLocalResidual's storage/flux with the ν̃ transport
 * terms (advective + diffusive), using the same face velocity as the advective flux of the
 * pressure equation.
 *
 * The production/destruction/cross-diffusion source terms are provided by
 * problem.source(...), see Dumux::RANSMassOneEqProblem in
 * dumux/freeflow/rans/oneeq/massproblem.hh.
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
    //! Upwinded advective ρν̃(v·n) plus a two-point-TPFA diffusive ((μ+ρν̃)/σ) term. The face
    //! velocity is obtained from problem.faceVelocity(element,fvGeometry,scvf), as for the
    //! advective flux of the pressure equation.
    //!
    //! \note Called for interior faces and for boundary faces with a pure-Dirichlet BC.
    //!       Neumann boundaries go through Problem::neumann() instead. For a boundary face,
    //!       scvf.outsideScvIdx() does not index a scv in fvGeometry (CCTpfa has no ghost
    //!       scvs), so the corresponding volVars (set from problem.dirichlet(...)) are looked up
    //!       directly from elemVars by that index, and the one-sided distance to the face is used.
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
        // distance-weighted mean of the diffusion coefficient on interior faces
        // (the inside/outside distances to the face are the weights)
        const auto weightedMean = [&](const auto inside, const auto outside)
        {
            if (scvf.boundary())
                return inside;
            const auto distanceInside = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            const auto distanceOutside = (fvGeometry.scv(scvf.outsideScvIdx()).dofPosition() - scvf.ipGlobal()).two_norm();
            return arithmeticMean(inside, outside, distanceOutside, distanceInside);
        };

        const auto diffCoeff = weightedMean(diffCoeffInside, diffCoeffOutside);

        flux -= diffCoeff*(outsideVars.viscosityTilde() - insideVars.viscosityTilde())/distance*area;

        return flux;
    }
};

} // end namespace Dumux

#endif
