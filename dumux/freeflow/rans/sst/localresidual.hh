// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::SSTMassLocalResidual
 */
#ifndef DUMUX_RANS_SST_MASS_LOCAL_RESIDUAL_HH
#define DUMUX_RANS_SST_MASS_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/navierstokes/mass/1p/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Element-wise storage/flux for the single-phase Navier-Stokes mass balance fused with
 *        the two-equation SST turbulence transport equations (k, omega).
 *
 * Extends Dumux::NavierStokesMassOnePLocalResidual's storage/flux with the k/omega storage
 * and diffusive flux terms. At boundary faces the "outside" (ghost) volume variables are
 * accessed via elemVars[scvf.outsideScvIdx()], because there is no sub control volume
 * in fvGeometry for them.
 *
 * \note The advective part of the k/omega flux is not added here, it is already contained in
 * the result of ParentType::computeFlux() through the Dumux::AdvectiveFlux specialization in
 * dumux/freeflow/rans/sst/advectiveflux.hh.
 *
 * The production/destruction/cross-diffusion source terms are provided by
 * Dumux::SSTMassProblem through problem.source(...).
 */
template<class TypeTag>
class SSTMassLocalResidual
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
        addTurbulenceFlux_(flux, problem, fvGeometry, elemVars, scvf);
        return flux;
    }

    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux = ParentType::fluxIntegral(fvGeometry, elemVars, scvf);
        addTurbulenceFlux_(flux, this->asImp().problem(), fvGeometry, elemVars, scvf);
        return flux;
    }

private:
    //! Adds the two-point diffusive fluxes of k and omega with the diffusion coefficients
    //! (mu + sigma_k*mu_t) and (mu + sigma_omega*mu_t), where sigma_k and sigma_omega are the
    //! blended coefficients of the selected model version (see volumevariables.hh).
    void addTurbulenceFlux_(NumEqVector& flux,
                            const Problem& problem,
                            const FVElementGeometry& fvGeometry,
                            const ElementVariables& elemVars,
                            const SubControlVolumeFace& scvf) const
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVars = elemVars[insideScv];
        const auto& outsideVars = elemVars[scvf.outsideScvIdx()];

        const auto area = Extrusion::area(fvGeometry, scvf)*insideVars.extrusionFactor();

        const auto diffCoeffKInside = insideVars.viscosity() + insideVars.sigmaK(problem)*insideVars.dynamicEddyViscosity();
        const auto diffCoeffKOutside = outsideVars.viscosity() + outsideVars.sigmaK(problem)*outsideVars.dynamicEddyViscosity();
        const auto diffCoeffOmegaInside = insideVars.viscosity() + insideVars.sigmaOmega(problem)*insideVars.dynamicEddyViscosity();
        const auto diffCoeffOmegaOutside = outsideVars.viscosity() + outsideVars.sigmaOmega(problem)*outsideVars.dynamicEddyViscosity();

        // The diffusion coefficient is interpolated to the face with distance-based weights,
        // which is more accurate than a plain average on graded meshes, in particular next to
        // the wall where the eddy viscosity varies strongly between neighboring cells.
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
