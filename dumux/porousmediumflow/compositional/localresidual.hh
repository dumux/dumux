// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 */

#ifndef DUMUX_COMPOSITIONAL_LOCAL_RESIDUAL_HH
#define DUMUX_COMPOSITIONAL_LOCAL_RESIDUAL_HH

#include <vector>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 */
template<class TypeTag>
class CompositionalLocalResidual: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = ModelTraits::useMoles();

    enum { conti0EqIdx = Indices::conti0EqIdx };

    //! The index of the component balance equation that gets replaced with the total mass balance
    static constexpr int replaceCompEqIdx = ModelTraits::replaceCompEqIdx();
    static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < numComponents;

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                auto eqIdx = conti0EqIdx + compIdx;
                if (eqIdx != replaceCompEqIdx)
                    storage[eqIdx] += volVars.porosity()
                                      * volVars.saturation(phaseIdx)
                                      * massOrMoleDensity(volVars, phaseIdx)
                                      * massOrMoleFraction(volVars, phaseIdx, compIdx);
            }

            // in case one balance is substituted by the total mole balance
            if (useTotalMoleOrMassBalance)
                storage[replaceCompEqIdx] += massOrMoleDensity(volVars, phaseIdx)
                                             * volVars.porosity()
                                             * volVars.saturation(phaseIdx);

            //! The energy storage in the fluid phase with index phaseIdx
            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }

        //! The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
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

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = conti0EqIdx + compIdx;

                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&massOrMoleDensity, &massOrMoleFraction, phaseIdx, compIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*massOrMoleFraction(volVars, phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                if (eqIdx != replaceCompEqIdx)
                    flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                // diffusive fluxes (only for the component balances)
                if(eqIdx != replaceCompEqIdx)
                {
                    //check for the reference system and adapt units of the diffusive flux accordingly.
                    if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                        flux[eqIdx] += useMoles ? diffusiveFluxes[compIdx]/FluidSystem::molarMass(compIdx)
                                            : diffusiveFluxes[compIdx];
                    else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                        flux[eqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                            : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
                    else
                        DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
                }
            }

            // in case one balance is substituted by the total mole balance
            if (useTotalMoleOrMassBalance)
            {
                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&massOrMoleDensity, phaseIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*volVars.mobility(phaseIdx); };

                flux[replaceCompEqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                for(int compIdx = 0; compIdx < numComponents; ++compIdx)
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
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);

            if constexpr (ModelTraits::enableCompositionalDispersion())
            {
                if constexpr (FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::box && numPhases == 1)
                {
                    const auto dispersionFluxes = fluxVars.compositionalDispersionFlux(phaseIdx);
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        flux[compIdx] += dispersionFluxes[compIdx];
                    }
                }
                else
                    DUNE_THROW(Dune::NotImplemented, "Dispersion Fluxes are only implemented for single phase flows using the Box method.");
            }

        }

        //! Add diffusive and dispersive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);
        EnergyLocalResidual::heatDispersionFlux(flux, fluxVars);

        return flux;
    }

protected:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }
};

} // end namespace Dumux

#endif
