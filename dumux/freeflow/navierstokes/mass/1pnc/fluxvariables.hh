// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassOnePNCFluxVariables
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1PNC_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MASS_1PNC_FLUXVARIABLES_HH

#include <dumux/common/typetraits/problem.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/upwindscheme.hh>
#include <dumux/freeflow/navierstokes/scalarfluxvariables.hh>
#include <dumux/discretization/extrusion.hh>

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
    static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < ModelTraits::numFluidComponents();

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

            // // in case one balance is substituted by the total mole balance
            // if constexpr(useTotalMoleOrMassBalance)
            // {
            //     for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            //     {
            //         //check for the reference system and adapt units of the diffusive flux accordingly.
            //         if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
            //             result[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx) : diffusiveFluxes[compIdx];
            //         else if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
            //             result[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]
            //                                     : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
            //         else
            //             DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
            //     }
            // }
        }

        return result;
    }

    /*!
     * \brief Returns the flux of enthalpy in J/s carried by diffusing molecules.
     */
    Scalar diffusiveEnthalpyFlux(const NumEqVector& diffusiveFlux, int phaseIdx = 0) const
    {
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

    template<class UpwindTerm>
    Scalar advectiveFluxForCellCenter(UpwindTerm upwindTerm) const
    {
        const auto& scvf = this->scvFace();
        const auto velocity = this->problem().faceVelocity(this->element(), this->fvGeometry(), scvf);
        const bool insideIsUpstream = !signbit(scvf.unitOuterNormal()*velocity);
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(this->problem().paramGroup(), "Flux.UpwindWeight");

        const auto& insideVolVars = this->elemVolVars()[scvf.insideScvIdx()];
        const auto& outsideVolVars = this->elemVolVars()[scvf.outsideScvIdx()];

        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? outsideVolVars : insideVolVars;

        using Extrusion = Extrusion_t<typename ProblemTraits<Problem>::GridGeometry>;
        const Scalar flux = (upwindWeight * upwindTerm(upstreamVolVars) +
                            (1.0 - upwindWeight) * upwindTerm(downstreamVolVars))
                            * velocity*scvf.unitOuterNormal() * Extrusion::area(scvf);

        return flux;
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
        // flux *= 0.0;
        // g++ requires to capture 'this' by value
        // const auto upwinding = [this](const auto& term) { return this->getAdvectiveFlux(term); };
        // AdvectiveFlux<ModelTraits>::addAdvectiveFlux(flux, upwinding);
        static constexpr bool useMoles = ModelTraits::useMoles();

        for (int compIdx = 0; compIdx < ModelTraits::numFluidComponents(); ++compIdx)
        {

        auto upwindTerm = [compIdx](const auto& volVars)
            {
                const auto density = useMoles ? volVars.molarDensity() : volVars.density();
                const auto fraction =  useMoles ? volVars.moleFraction(compIdx) : volVars.massFraction(compIdx);

                return density * fraction;
            };

            flux[compIdx] += advectiveFluxForCellCenter(upwindTerm);

        }




        // in case one balance is substituted by the total mass balance
        if constexpr (ModelTraits::replaceCompEqIdx() < ModelTraits::numFluidComponents())
        {
            flux[ModelTraits::replaceCompEqIdx()] = std::accumulate(flux.begin(), flux.end(), 0.0);
        }

        if constexpr (ModelTraits::enableEnergyBalance())
        {
            ParentType::addHeatFlux(flux);
            flux[ModelTraits::Indices::energyEqIdx] += diffusiveEnthalpyFlux(diffusiveFlux);
        }

        return flux;
    }

};

} // end namespace Dumux

#endif
