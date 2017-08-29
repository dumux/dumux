// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Element-wise calculation of the residual and its derivatives
 *        for a two-phase, incompressible test problem.
 */
#ifndef DUMUX_2P_INCOMPRESSIBLE_TEST_LOCAL_RESIDUAL_HH
#define DUMUX_2P_INCOMPRESSIBLE_TEST_LOCAL_RESIDUAL_HH

#include <dumux/porousmediumflow/immiscible/localresidual.hh>

namespace Dumux
{

template<class TypeTag>
class TwoPIncompressibleLocalResidual : public ImmiscibleLocalResidual<TypeTag>
{
    using ParentType = ImmiscibleLocalResidual<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementResidualVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    // first index for the mass balance
    enum
    {
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx
    };

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

public:
    using ParentType::ParentType;

    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars) const
    {
        // we know that these values are constant throughout the simulation
        static const auto phi = curVolVars.porosity();
        static const auto phi_rho_w = phi*curVolVars.density(FluidSystem::wPhaseIdx);
        static const auto phi_rho_n = phi*curVolVars.density(FluidSystem::nPhaseIdx);

        const auto volume = element.geometry().volume();

        // partial derivative of wetting phase storage term w.r.t. p_w
        partialDerivatives[contiWEqIdx][pressureIdx] += 0.0;
        // partial derivative of wetting phase storage term w.r.t. S_n
        partialDerivatives[contiWEqIdx][saturationIdx] -= volume*phi_rho_w/this->timeLoop().timeStepSize();
        // partial derivative of non-wetting phase storage term w.r.t. p_w
        partialDerivatives[contiNEqIdx][pressureIdx] += 0.0;
        // partial derivative of non-wetting phase storage term w.r.t. S_n
        partialDerivatives[contiNEqIdx][saturationIdx] += volume*phi_rho_n/this->timeLoop().timeStepSize();
    }

    template<class PartialDerivativeMatrix>
    void addSourceDerivatives(PartialDerivativeMatrix& partialDerivatives,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& curVolVars) const {}

    template<class PartialDerivativeMatrices>
    void addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                            const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& curElemVolVars,
                            const ElementFluxVariablesCache& elemFluxVarsCache,
                            const SubControlVolumeFace& scvf) const
    {
        using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
        static const Scalar upwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, UpwindWeight);

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        const auto wettingflux = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                     scvf, FluidSystem::wPhaseIdx, elemFluxVarsCache);
        const auto insideUpwindWeightWetting = std::signbit(wettingflux) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideUpwindWeightWetting = 1.0 - insideUpwindWeightWetting;

        const auto nonwettingflux = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                        scvf, FluidSystem::nPhaseIdx, elemFluxVarsCache);
        const auto insideUpwindWeightNonWetting = std::signbit(nonwettingflux) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideUpwindWeightNonWetting = 1.0 - insideUpwindWeightNonWetting;

        // get references to the two participating vol vars & partial derivative matrices
        const auto& insideVolVars = curElemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = curElemVolVars[scvf.outsideScvIdx()];
        auto& dI_dJ_inside = derivativeMatrices[scvf.insideScvIdx()];
        auto& dI_dJ_outside = derivativeMatrices[scvf.outsideScvIdx()];

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho_w = insideVolVars.density(FluidSystem::wPhaseIdx);
        static const auto rho_n = insideVolVars.density(FluidSystem::nPhaseIdx);
        static const auto rhow_muw = rho_w/insideVolVars.viscosity(FluidSystem::wPhaseIdx);
        static const auto rhon_mun = rho_n/insideVolVars.viscosity(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_inside = rho_w*insideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_inside = rho_n*insideVolVars.mobility(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_outside = rho_w*outsideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_outside = rho_n*outsideVolVars.mobility(FluidSystem::nPhaseIdx);

        // we know the material parameters are only position-dependent
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideMaterialParams = problem.spatialParams().materialLawParamsAtPos(insideScv.center());
        const auto& outsideMaterialParams = problem.spatialParams().materialLawParamsAtPos(outsideScv.center());

        // derivative w.r.t. to Sn is the negative of the one w.r.t. Sw
        const auto insideSw = insideVolVars.saturation(FluidSystem::wPhaseIdx);
        const auto outsideSw = outsideVolVars.saturation(FluidSystem::wPhaseIdx);

        const auto dKrw_dSn_inside = -1.0*MaterialLaw::dkrw_dsw(insideMaterialParams, insideSw);
        const auto dKrw_dSn_outside = -1.0*MaterialLaw::dkrw_dsw(outsideMaterialParams, outsideSw);
        const auto dKrn_dSn_inside = -1.0*MaterialLaw::dkrn_dsw(insideMaterialParams, insideSw);
        const auto dKrn_dSn_outside = -1.0*MaterialLaw::dkrn_dsw(outsideMaterialParams, outsideSw);

        const auto dpc_dSn_inside = -1.0*MaterialLaw::dpc_dsw(insideMaterialParams, insideSw);
        const auto dpc_dSn_outside = -1.0*MaterialLaw::dpc_dsw(outsideMaterialParams, outsideSw);

        const auto tij = elemFluxVarsCache[scvf].advectionTij();
        // partial derivative of the wetting phase flux w.r.t. p_w
        dI_dJ_inside[contiWEqIdx][pressureIdx] += tij*(rhowKrw_muw_inside*insideUpwindWeightWetting
                                                       + rhowKrw_muw_outside*outsideUpwindWeightWetting);
        dI_dJ_outside[contiWEqIdx][pressureIdx] -= tij*(rhowKrw_muw_inside*insideUpwindWeightWetting
                                                        + rhowKrw_muw_outside*outsideUpwindWeightWetting);

        // partial derivative of the wetting phase flux w.r.t. S_n
        dI_dJ_inside[contiWEqIdx][saturationIdx] += rhow_muw*wettingflux*dKrw_dSn_inside*insideUpwindWeightWetting;
        dI_dJ_outside[contiWEqIdx][saturationIdx] += rhow_muw*wettingflux*dKrw_dSn_outside*outsideUpwindWeightWetting;

        // partial derivative of the non-wetting phase flux w.r.t. p_w
        dI_dJ_inside[contiNEqIdx][pressureIdx] += tij*(rhonKrn_mun_inside*insideUpwindWeightNonWetting
                                                       + rhonKrn_mun_outside*outsideUpwindWeightNonWetting);
        dI_dJ_outside[contiNEqIdx][pressureIdx] -= tij*(rhonKrn_mun_inside*insideUpwindWeightNonWetting
                                                        + rhonKrn_mun_outside*outsideUpwindWeightNonWetting);

        // partial derivative of the non-wetting phase flux w.r.t. S_n (relative permeability derivative contribution)
        dI_dJ_inside[contiNEqIdx][saturationIdx] += rhon_mun*nonwettingflux*dKrn_dSn_inside*insideUpwindWeightNonWetting;
        dI_dJ_outside[contiNEqIdx][saturationIdx] += rhon_mun*nonwettingflux*dKrn_dSn_outside*outsideUpwindWeightNonWetting;

        // partial derivative of the non-wetting phase flux w.r.t. S_n (capillary pressure derivative contribution)
        dI_dJ_inside[contiNEqIdx][saturationIdx] += tij*dpc_dSn_inside*(rhonKrn_mun_inside*insideUpwindWeightNonWetting
                                                                        + rhonKrn_mun_outside*outsideUpwindWeightNonWetting);
        dI_dJ_outside[contiNEqIdx][saturationIdx] -= tij*dpc_dSn_outside*(rhonKrn_mun_inside*insideUpwindWeightNonWetting
                                                                          + rhonKrn_mun_outside*outsideUpwindWeightNonWetting);
    }

    template<class PartialDerivativeMatrices>
    void addDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& curElemVolVars,
                                     const ElementFluxVariablesCache& elemFluxVarsCache,
                                     const SubControlVolumeFace& scvf) const
    {
        using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
        static const Scalar upwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, UpwindWeight);

        // evaluate the current wetting phase Darcy flux and resulting upwind weights
        const auto wettingflux = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                     scvf, FluidSystem::wPhaseIdx, elemFluxVarsCache);
        const auto insideUpwindWeightWetting = std::signbit(wettingflux) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideUpwindWeightWetting = 1.0 - insideUpwindWeightWetting;

        const auto nonwettingflux = AdvectionType::flux(problem, element, fvGeometry, curElemVolVars,
                                                        scvf, FluidSystem::nPhaseIdx, elemFluxVarsCache);
        const auto insideUpwindWeightNonWetting = std::signbit(nonwettingflux) ? (1.0 - upwindWeight) : upwindWeight;
        const auto outsideUpwindWeightNonWetting = 1.0 - insideUpwindWeightNonWetting;

        // get references to the two participating vol vars & partial derivative matrices
        const auto& insideVolVars = curElemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = curElemVolVars[scvf.outsideScvIdx()];
        auto& dI_dJ_inside = derivativeMatrices[scvf.insideScvIdx()];

        // some quantities to be reused (rho & mu are constant and thus equal for all cells)
        static const auto rho_w = insideVolVars.density(FluidSystem::wPhaseIdx);
        static const auto rho_n = insideVolVars.density(FluidSystem::nPhaseIdx);
        static const auto rhow_muw = rho_w/insideVolVars.viscosity(FluidSystem::wPhaseIdx);
        static const auto rhon_mun = rho_n/insideVolVars.viscosity(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_inside = rho_w*insideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_inside = rho_n*insideVolVars.mobility(FluidSystem::nPhaseIdx);
        const auto rhowKrw_muw_outside = rho_w*outsideVolVars.mobility(FluidSystem::wPhaseIdx);
        const auto rhonKrn_mun_outside = rho_n*outsideVolVars.mobility(FluidSystem::nPhaseIdx);

        // we know the material parameters are only position-dependent
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideMaterialParams = problem.spatialParams().materialLawParamsAtPos(insideScv.center());

        // derivative w.r.t. to Sn is the negative of the one w.r.t. Sw
        const auto insideSw = insideVolVars.saturation(FluidSystem::wPhaseIdx);
        const auto dKrw_dSn_inside = -1.0*MaterialLaw::dkrw_dsw(insideMaterialParams, insideSw);
        const auto dKrn_dSn_inside = -1.0*MaterialLaw::dkrn_dsw(insideMaterialParams, insideSw);
        const auto dpc_dSn_inside = -1.0*MaterialLaw::dpc_dsw(insideMaterialParams, insideSw);

        const auto tij = elemFluxVarsCache[scvf].advectionTij();
        // partial derivative of the wetting phase flux w.r.t. p_w
        dI_dJ_inside[contiWEqIdx][pressureIdx] += tij*(rhowKrw_muw_inside*insideUpwindWeightWetting
                                                       + rhowKrw_muw_outside*outsideUpwindWeightWetting);

        // partial derivative of the wetting phase flux w.r.t. S_n
        dI_dJ_inside[contiWEqIdx][saturationIdx] += rhow_muw*wettingflux*dKrw_dSn_inside*insideUpwindWeightWetting;

        // partial derivative of the non-wetting phase flux w.r.t. p_w
        dI_dJ_inside[contiNEqIdx][pressureIdx] += tij*(rhonKrn_mun_inside*insideUpwindWeightNonWetting
                                                       + rhonKrn_mun_outside*outsideUpwindWeightNonWetting);

        // partial derivative of the non-wetting phase flux w.r.t. S_n (relative permeability derivative contribution)
        dI_dJ_inside[contiNEqIdx][saturationIdx] += rhon_mun*nonwettingflux*dKrn_dSn_inside*insideUpwindWeightNonWetting;

        // partial derivative of the non-wetting phase flux w.r.t. S_n (capillary pressure derivative contribution)
        dI_dJ_inside[contiNEqIdx][saturationIdx] += tij*dpc_dSn_inside*(rhonKrn_mun_inside*insideUpwindWeightNonWetting
                                                                        + rhonKrn_mun_outside*outsideUpwindWeightNonWetting);
    }

    template<class PartialDerivativeMatrices>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    {
        // we have no Robin-type boundary conditions here. Do nothing.
    }
};

} // end namespace Dumux

#endif
