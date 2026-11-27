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

#include <dumux/common/math.hh>
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
 * diffusive transport terms. The volume variables on the outside of a boundary face are
 * fetched via elemVars[scvf.outsideScvIdx()] directly.
 *
 * \note The advective part of the k/epsilon-tilde flux is not added here: it is already
 * contained in the result of ParentType::computeFlux() through the Dumux::AdvectiveFlux
 * specialization for LowReKEpsilonMassModelTraits
 * (dumux/freeflow/rans/lowrekepsilon/advectiveflux.hh), which is applied to every equation index.
 *
 * Production/destruction/damping-function source terms are not added here: they are provided
 * by problem.source(...) (Dumux::LowReKEpsilonMassProblem).
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
    //! Adds the two-point diffusive terms with the diffusion coefficients (mu_t/sigma_k + mu) and
    //! (mu_t/sigma_epsilon + mu). The advective part is not added here, see the class
    //! documentation.
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

        // mean of the diffusion coefficient on interior faces, weighted by the distances of the
        // two cell centers to the face
        const auto weightedMean = [&](const auto inside, const auto outside)
        {
            if (scvf.boundary())
                return inside;
            const auto distanceInside = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            const auto distanceOutside = (fvGeometry.scv(scvf.outsideScvIdx()).dofPosition() - scvf.ipGlobal()).two_norm();
            return arithmeticMean(inside, outside, distanceOutside, distanceInside);
        };

        const auto diffCoeffK = weightedMean(diffCoeffKInside, diffCoeffKOutside);
        const auto diffCoeffEpsilon = weightedMean(diffCoeffEpsilonInside, diffCoeffEpsilonOutside);

        flux[Indices::turbulentKineticEnergyEqIdx]
            += diffCoeffK*(insideVars.turbulentKineticEnergy() - outsideVars.turbulentKineticEnergy())/distance*area;
        flux[Indices::dissipationEqIdx]
            += diffCoeffEpsilon*(insideVars.dissipationTilde() - outsideVars.dissipationTilde())/distance*area;
    }
};

} // end namespace Dumux

#endif
