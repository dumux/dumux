// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 */
template<class TypeTag>
class CompositionalLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = typename GridGeometry::GridView;
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
            EnergyLocalResidual::fluidPhaseStorage(storage, problem, scv, volVars, phaseIdx);
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
        static constexpr auto referenceSystemFormulationDiffusion = FluxVariables::MolecularDiffusionType::referenceSystemFormulation();
        // get upwind weights into local scope
        NumEqVector flux(0.0);

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction = [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // Check for the reference system and adapt units of the flux accordingly.
        const auto adaptFluxUnits = [](const Scalar referenceFlux, const Scalar molarMass,
                                       const ReferenceSystemFormulation referenceSystemFormulation)
        {
            if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                return useMoles ? referenceFlux/molarMass
                                : referenceFlux;
            else if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged)
                return useMoles ? referenceFlux
                                : referenceFlux*molarMass;
            else
                DUNE_THROW(Dune::NotImplemented, "other reference systems than mass and molar averaged are not implemented");
        };

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
            auto dispersiveFluxes = decltype(diffusiveFluxes)(0.0);
            if constexpr (ModelTraits::enableCompositionalDispersion())
                dispersiveFluxes = fluxVars.compositionalDispersionFlux(phaseIdx);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = conti0EqIdx + compIdx;

                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&massOrMoleDensity, &massOrMoleFraction, phaseIdx, compIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*massOrMoleFraction(volVars, phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                // advective fluxes (only for the component balances)
                if (eqIdx != replaceCompEqIdx)
                    flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                // diffusive and dispersive fluxes
                auto diffusiveAndDispersiveFluxes = adaptFluxUnits(diffusiveFluxes[compIdx],
                                                                   FluidSystem::molarMass(compIdx),
                                                                   referenceSystemFormulationDiffusion);
                if constexpr (ModelTraits::enableCompositionalDispersion())
                {
                    static constexpr auto referenceSystemFormulationDispersion =
                        FluxVariables::DispersionFluxType::referenceSystemFormulation();
                    diffusiveAndDispersiveFluxes += adaptFluxUnits(dispersiveFluxes[compIdx],
                                                                   FluidSystem::molarMass(compIdx),
                                                                   referenceSystemFormulationDispersion);
                }
                if(eqIdx != replaceCompEqIdx)
                    flux[eqIdx] += diffusiveAndDispersiveFluxes;
                if (useTotalMoleOrMassBalance)
                    flux[replaceCompEqIdx] += diffusiveAndDispersiveFluxes;
            }

            // in case one balance is substituted by the total mole balance
            if (useTotalMoleOrMassBalance)
            {
                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&massOrMoleDensity, phaseIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*volVars.mobility(phaseIdx); };

                flux[replaceCompEqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
            }

            //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        //! Add diffusive and dispersive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);
        EnergyLocalResidual::heatDispersionFlux(flux, fluxVars);

        return flux;
    }
};

} // end namespace Dumux

#endif
