// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassOnePNCFluxVariables
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1PNC_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MASS_1PNC_FLUXVARIABLES_HH

#include <dumux/common/typetraits/problem.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/upwindscheme.hh>
#include <dumux/freeflow/navierstokes/scalarfluxvariables.hh>

#include "advectiveflux.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the single-phase flow,
 *        multi-component Navier-Stokes model.
 */
template<class Problem,
         class ModelTraits,
         class FluxTs,
         class ElementVolumeVariables,
         class ElementFluxVariablesCache,
         class UpwindScheme = UpwindScheme<typename ProblemTraits<Problem>::GridGeometry>>
class NavierStokesMassOnePNCFluxVariables
: public NavierStokesScalarConservationModelFluxVariables<Problem,
                                                          ModelTraits,
                                                          FluxTs,
                                                          ElementVolumeVariables,
                                                          ElementFluxVariablesCache,
                                                          UpwindScheme>
{
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
    using NumEqVector = typename VolumeVariables::PrimaryVariables;
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
    using Indices = typename ModelTraits::Indices;

    static constexpr bool enableMolecularDiffusion = ModelTraits::enableMolecularDiffusion();
    static constexpr auto replaceCompEqIdx = ModelTraits::replaceCompEqIdx();
    static constexpr bool useTotalMassBalance = replaceCompEqIdx < ModelTraits::numFluidComponents();

    using FluidSystem = typename VolumeVariables::FluidSystem;

    using ParentType = NavierStokesScalarConservationModelFluxVariables<Problem,
                                                                        ModelTraits,
                                                                        FluxTs,
                                                                        ElementVolumeVariables,
                                                                        ElementFluxVariablesCache,
                                                                        UpwindScheme>;

public:

    static constexpr auto numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = ModelTraits::useMoles();
    using MolecularDiffusionType = typename FluxTs::MolecularDiffusionType;

    /*!
     * \brief Returns the diffusive fluxes computed by the respective law.
     */
    NumEqVector molecularDiffusionFlux(int phaseIdx = 0) const
    {
        NumEqVector result(0.0);
        if constexpr (enableMolecularDiffusion)
        {
            const auto diffusiveFluxes = MolecularDiffusionType::flux(this->problem(),
                                                                      this->element(),
                                                                      this->fvGeometry(),
                                                                      this->elemVolVars(),
                                                                      this->scvFace(),
                                                                      phaseIdx,
                                                                      this->elemFluxVarsCache());

            static constexpr auto referenceSystemFormulation = MolecularDiffusionType::referenceSystemFormulation();

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = Indices::conti0EqIdx + compIdx;
                if (eqIdx == replaceCompEqIdx)
                    continue;

                //check for the reference system and adapt units of the diffusive flux accordingly.
                if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                    result[eqIdx] += useMoles ? diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx)
                                        : diffusiveFluxes[compIdx];
                else if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                    result[eqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                        : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
                else
                    DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
            }

            // in case one balance is substituted by the total mass balance
            if constexpr(useTotalMassBalance)
            {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    //check for the reference system and adapt units of the diffusive flux accordingly.
                    if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                        result[replaceCompEqIdx] += diffusiveFluxes[compIdx];
                    else if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                        result[replaceCompEqIdx] += diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
                    else
                        DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
                }
            }
        }

        return result;
    }

    /*!
     * \brief Returns the flux of enthalpy in J/s carried by diffusing molecules.
     */
    Scalar diffusiveEnthalpyFlux(int phaseIdx = 0) const
    {
        const auto diffusiveFlux = MolecularDiffusionType::flux(this->problem(),
                                                                this->element(),
                                                                this->fvGeometry(),
                                                                this->elemVolVars(),
                                                                this->scvFace(),
                                                                phaseIdx,
                                                                this->elemFluxVarsCache());
        static constexpr auto referenceSystemFormulation = MolecularDiffusionType::referenceSystemFormulation();
        Scalar flux = 0.0;
        const auto& scvf = this->scvFace();
        const auto& elemVolVars = this->elemVolVars();

        const auto componentEnthalpy = [](const auto& volVars, int compIdx)
        { return FluidSystem::componentEnthalpy(volVars.fluidState(), 0, compIdx); };

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // do a simple upwinding of enthalpy based on the direction of the diffusive flux
            using std::signbit;
            const bool insideIsUpstream = !signbit(diffusiveFlux[compIdx]);
            const auto& upstreamVolVars = insideIsUpstream ? elemVolVars[scvf.insideScvIdx()] : elemVolVars[scvf.outsideScvIdx()];

            if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                flux += diffusiveFlux[compIdx] * componentEnthalpy(upstreamVolVars, compIdx);
            else
                flux += diffusiveFlux[compIdx] * componentEnthalpy(upstreamVolVars, compIdx)* elemVolVars[scvf.insideScvIdx()].molarMass(compIdx);
        }

        return flux;
    }

    /*!
     * \brief Returns the advective mass flux in kg/s
     *        or the advective mole flux in mole/s.
     */
    NumEqVector advectiveFlux(int phaseIdx = 0) const
    {
        NumEqVector result(0.0);
        // g++ requires to capture 'this' by value
        const auto upwinding = [this](const auto& term) { return this->getAdvectiveFlux(term); };
        AdvectiveFlux<ModelTraits>::addAdvectiveFlux(result, upwinding);
        return result;
    }

    /*!
     * \brief Returns all fluxes for the single-phase flow, multi-component
     *        Navier-Stokes model: the advective mass flux in kg/s
     *        or the advective mole flux in mole/s and the energy flux
     *        in J/s (for nonisothermal models).
     */
    NumEqVector flux(int phaseIdx = 0) const
    {
        const auto diffusiveFlux = molecularDiffusionFlux(phaseIdx);
        NumEqVector flux = diffusiveFlux;
        // g++ requires to capture 'this' by value
        const auto upwinding = [this](const auto& term) { return this->getAdvectiveFlux(term); };
        AdvectiveFlux<ModelTraits>::addAdvectiveFlux(flux, upwinding);

        if constexpr (ModelTraits::enableEnergyBalance())
        {
            ParentType::addHeatFlux(flux);
            flux[ModelTraits::Indices::energyEqIdx] += diffusiveEnthalpyFlux();
        }

        return flux;
    }

};

} // end namespace Dumux

#endif
