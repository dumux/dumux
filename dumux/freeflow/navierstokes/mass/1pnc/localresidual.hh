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
 * \copydoc Dumux::NavierStokesResidualImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1P_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MASS_1P_LOCAL_RESIDUAL_HH

#include <dumux/assembly/cclocalresidual.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/freeflow/navierstokes/energy/localresidual.hh>

namespace Dumux {


/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using the staggered discretization
 */
template<class TypeTag>
class NavierStokesMassOnePNCLocalResidual : public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using EnergyLocalResidual = NavierStokesEnergyLocalResidual<ModelTraits>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static_assert(GridGeometry::discMethod == DiscretizationMethod::cctpfa);
    static constexpr bool useMoles = ModelTraits::useMoles();
    static constexpr auto numComponents = ModelTraits::numFluidComponents();

    static constexpr int replaceCompEqIdx = ModelTraits::replaceCompEqIdx();
    static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < numComponents;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param scv The sub-control volume over which we integrate the storage term
     * \param volVars The volume variables associated with the scv
     * \note has to be implemented by the model specific residual class
     *
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        const Scalar density = useMoles ? volVars.molarDensity() : volVars.density();

        // compute storage term of all components
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int eqIdx = compIdx;

            const Scalar massOrMoleFraction = useMoles ? volVars.moleFraction(compIdx) : volVars.massFraction(compIdx);
            const Scalar s =  density * massOrMoleFraction;

            if (eqIdx != ModelTraits::replaceCompEqIdx())
                storage[eqIdx] += s;
        }

        // in case one balance is substituted by the total mass balance
        if(ModelTraits::replaceCompEqIdx() < numComponents)
            storage[ModelTraits::replaceCompEqIdx()] = density;

        EnergyLocalResidual::fluidPhaseStorage(storage, volVars);

        return storage;
    }

    /*!
     * \brief Evaluatex the mass flux over a face of a sub control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        static constexpr auto referenceSystemFormulation = FluxVariables::MolecularDiffusionType::referenceSystemFormulation();
        // get upwind weights into local scope
        NumEqVector flux(0.0);
        static constexpr auto phaseIdx = 0;

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // get the diffusive flux of each component
        const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // get equation index
            const auto eqIdx = Indices::conti0EqIdx + compIdx;

            // the physical quantities for which we perform upwinding
            const auto upwindTerm = [&massOrMoleDensity, &massOrMoleFraction, phaseIdx = phaseIdx, compIdx] (const auto& volVars)
            { return massOrMoleDensity(volVars, phaseIdx)*massOrMoleFraction(volVars, phaseIdx, compIdx); };

            if (eqIdx != replaceCompEqIdx)
            {
                // advective fluxes
                flux[eqIdx] += fluxVars.advectiveFlux(upwindTerm);

                // diffusive fluxes (only for the component balances)
                //check for the reference system and adapt units of the diffusive flux accordingly.
                if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                    flux[eqIdx] += useMoles ? diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx)
                                        : diffusiveFluxes[compIdx];
                else if constexpr (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                    flux[eqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                        : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
                else
                    DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
            }
        }

        // in case one balance is substituted by the total mole balance
        if constexpr(useTotalMoleOrMassBalance)
        {
            // the physical quantities for which we perform upwinding
            const auto upwindTerm = [&massOrMoleDensity, phaseIdx = phaseIdx] (const auto& volVars)
            { return massOrMoleDensity(volVars, phaseIdx); };

            flux[replaceCompEqIdx] += fluxVars.advectiveFlux(upwindTerm);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                //check for the reference system and adapt units of the diffusive flux accordingly.
                if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                    flux[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx) : diffusiveFluxes[compIdx];
                else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                    flux[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                            : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
                else
                    DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
            }
        }

        //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConvectionFlux(flux, fluxVars);


        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;
    }
};

} // end namespace Dumux

#endif
