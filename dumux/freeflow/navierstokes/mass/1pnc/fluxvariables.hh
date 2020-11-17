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
#include <dumux/flux/upwindscheme.hh>
#include <dumux/freeflow/navierstokes/scalarfluxvariables.hh>

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

            // in case one balance is substituted by the total mole balance
            if constexpr(useTotalMoleOrMassBalance)
            {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    //check for the reference system and adapt units of the diffusive flux accordingly.
                    if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                        result[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx) : diffusiveFluxes[compIdx];
                    else if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                        result[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                                : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
                    else
                        DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
                }
            }
        }

        return result;
    }

    /*!
     * \brief Returns the advective mass flux in kg/s
     *        or the advective mole flux in mole/s.
     */
    NumEqVector advectiveFlux(int phaseIdx = 0) const
    {
        NumEqVector result(0.0);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // get equation index
            const auto eqIdx = Indices::conti0EqIdx + compIdx;

            if (eqIdx != replaceCompEqIdx)
            {
                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [this, phaseIdx = phaseIdx, compIdx] (const auto& volVars)
                { return massOrMoleDensity_(volVars, phaseIdx, compIdx) * massOrMoleFraction_(volVars, phaseIdx, compIdx); };

                // advective fluxes
                result[eqIdx] += ParentType::advectiveFlux(upwindTerm);
            }
            else
            {
                // in case one balance is substituted by the total mole balance
                if constexpr(useTotalMoleOrMassBalance)
                {
                    // the physical quantities for which we perform upwinding
                    const auto upwindTerm = [this, phaseIdx = phaseIdx, compIdx] (const auto& volVars)
                    { return massOrMoleDensity_(volVars, phaseIdx, compIdx); };

                    result[replaceCompEqIdx] += ParentType::advectiveFlux(upwindTerm);
                }
            }
        }

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
        NumEqVector flux = molecularDiffusionFlux(phaseIdx) + advectiveFlux(phaseIdx);
        ParentType::addHeatFlux(flux);
        return flux;
    }

private:
    Scalar massOrMoleDensity_(const VolumeVariables& volVars, const int phaseIdx, const int compIdx) const
    { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

    Scalar massOrMoleFraction_(const VolumeVariables& volVars, const int phaseIdx, const int compIdx) const
    { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };
};

} // end namespace Dumux

#endif
