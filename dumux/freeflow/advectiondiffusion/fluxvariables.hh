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
 * \ingroup AdvectionDiffusionModel
 */
#ifndef DUMUX_ADVECTION_DIFFUSION_FLUXVARIABLES_HH
#define DUMUX_ADVECTION_DIFFUSION_FLUXVARIABLES_HH

#include <array>
#include <optional>
#include <type_traits>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/flux/upwindscheme.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup AdvectionDiffusionModel
 * \brief The flux variables class for the advection diffusion transport model
 */

template<class TypeTag, class UpScheme = UpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>> >
class AdvectionDiffusionFluxVariables
: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                  typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView,
                  typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                  typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
public:
    using UpwindScheme = UpScheme;
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
    using MolecularDiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static constexpr auto numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = ModelTraits::useMoles();
    static constexpr bool enableAdvection = ModelTraits::enableAdvection();
    static constexpr bool enableMolecularDiffusion = ModelTraits::enableMolecularDiffusion();

    //! The constructor
    AdvectionDiffusionFluxVariables()
    { }

    /*!
     * \brief Returns the advective mass flux in kg/s
     *        or the advective mole flux in mole/s.
     */
    NumEqVector advectiveFlux(int phaseIdx = 0) const
    {
        static_assert(enableAdvection);

        NumEqVector advFlux(0.0);
        // formulation with mole balances
        if (useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // the physical quantities for which we perform upwinding
                auto upwindTerm = [compIdx](const auto& volVars)
                { return volVars.molarDensity() * volVars.moleFraction(0, compIdx); };

                const Scalar advFluxBeforeUpwind = AdvectionType::flux(this->problem(),
                                                                       this->element(),
                                                                       this->fvGeometry(),
                                                                       this->elemVolVars(),
                                                                       this->scvFace(),
                                                                       phaseIdx,
                                                                       this->elemFluxVarsCache());
                // advective fluxes
                advFlux[compIdx] += UpwindScheme::apply(*this, upwindTerm, advFluxBeforeUpwind, phaseIdx);
            }
        }
        // formulation with mass balances
        else
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // the physical quantities for which we perform upwinding
                auto upwindTerm = [compIdx](const auto& volVars)
                { return volVars.density()*volVars.massFraction(0, compIdx); };


                const Scalar advFluxBeforeUpwind = AdvectionType::flux(this->problem(),
                                                                       this->element(),
                                                                       this->fvGeometry(),
                                                                       this->elemVolVars(),
                                                                       this->scvFace(),
                                                                       phaseIdx,
                                                                       this->elemFluxVarsCache());
                // advective fluxes
                advFlux[compIdx] += UpwindScheme::apply(*this, upwindTerm, advFluxBeforeUpwind, phaseIdx);
            }
        }

        return advFlux;
    }


    /*!
     * \brief Returns the diffusive fluxes computed by the respective law.
     */
    NumEqVector molecularDiffusionFlux([[maybe_unused]] const int phaseIdx = 0) const
    {
        static_assert(enableMolecularDiffusion);

        NumEqVector flux(0.0);

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
            //check for the reference system and adapt units of the diffusive flux accordingly.
            if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                flux[compIdx] += useMoles ? diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx) : diffusiveFluxes[compIdx];
            else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                flux[compIdx] += useMoles ? diffusiveFluxes[compIdx] : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
            else
                DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
        }

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
        NumEqVector flux(0.0);
        if constexpr (enableMolecularDiffusion)
            flux += molecularDiffusionFlux(phaseIdx);
        if constexpr (enableAdvection)
            flux += advectiveFlux(phaseIdx);
        return flux;
    }

private:

};

} // end namespace Dumux

#endif
