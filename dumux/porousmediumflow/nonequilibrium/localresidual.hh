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
 * \ingroup NonEquilibriumModel
 * \brief The local residual for the kinetic mass transfer module of
 *        the compositional multi-phase model.
 */

#ifndef DUMUX_NONEQUILIBRIUM_LOCAL_RESIDUAL_HH
#define DUMUX_NONEQUILIBRIUM_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/nonequilibrium/thermal/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, bool enableChemicalNonEquilibrium>
class NonEquilibriumLocalResidualImplementation;

template <class TypeTag>
using NonEquilibriumLocalResidual = NonEquilibriumLocalResidualImplementation<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableChemicalNonEquilibrium()>;

/*!
 * \ingroup NonEquilibriumModel
 * \brief The mass conservation part of the nonequilibrium model for a model without chemical non-equilibrium
 */
template<class TypeTag>
class NonEquilibriumLocalResidualImplementation<TypeTag, false>: public GetPropType<TypeTag, Properties::EquilibriumLocalResidual>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ParentType = GetPropType<TypeTag, Properties::EquilibriumLocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int numComponents = ModelTraits::numFluidComponents();
    enum { conti0EqIdx = Indices::conti0EqIdx };
public:
    using ParentType::ParentType;

    /*!
     * \brief Calculates the source term of the equation.
     *
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux compuation
     *
     * \note This is the default implementation for all models as sources are computed
     *       in the user interface of the problem
     *
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
        // get upwind weights into local scope
        NumEqVector flux(0.0);

        const auto moleDensity = [](const auto& volVars, const int phaseIdx)
        { return volVars.molarDensity(phaseIdx); };

        const auto moleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return  volVars.moleFraction(phaseIdx, compIdx); };

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = conti0EqIdx + compIdx;

                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&moleDensity, &moleFraction, phaseIdx, compIdx] (const auto& volVars)
                { return moleDensity(volVars, phaseIdx)*moleFraction(volVars, phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                    flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                // diffusive fluxes (only for the component balances)
                    flux[eqIdx] += diffusiveFluxes[compIdx];
            }

            //! Add advective and diffusive phase energy fluxes
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }
        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;


    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);
        // Call the (kinetic) Energy module, for the source term.
        // it has to be called from here, because the mass transfered has to be known.
        EnergyLocalResidual::computeSourceEnergy(source,
                                                 element,
                                                 fvGeometry,
                                                 elemVolVars,
                                                 scv);
        return source;
    }

};

/*!
 * \brief The mass conservation part of the nonequilibrium model for a model assuming chemical non-equilibrium and two phases but assuming thermal equilibrium
 */
template<class TypeTag>
class NonEquilibriumLocalResidualImplementation<TypeTag, true>: public GetPropType<TypeTag, Properties::EquilibriumLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::EquilibriumLocalResidual>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MolecularDiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

    static constexpr auto conti0EqIdx = Indices::conti0EqIdx;
    static constexpr auto comp1Idx = FluidSystem::comp1Idx;
    static constexpr auto comp0Idx = FluidSystem::comp0Idx;
    static constexpr auto phase0Idx = FluidSystem::phase0Idx;
    static constexpr auto phase1Idx = FluidSystem::phase1Idx;

    static_assert(numPhases == numComponents,
                  "currently chemical non-equilibrium is only available when numPhases equals numComponents");
    static_assert(useMoles == true,
                  "chemical nonequilibrium can only be calculated based on mole fractions not mass fractions");

public:
     using ParentType::ParentType;
    /*!
     * \brief Calculates the storage for all mass balance equations.
     *
     * \param problem The object specifying the problem which ought to be simulated
     * \param scv The sub-control volume
     * \param volVars The volume variables
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
       NumEqVector storage(0.0);

        // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                 auto eqIdx = conti0EqIdx + phaseIdx*numComponents + compIdx;
                 storage[eqIdx] += volVars.porosity()
                                   * volVars.saturation(phaseIdx)
                                   * volVars.molarDensity(phaseIdx)
                                   * volVars.moleFraction(phaseIdx, compIdx);
            }
            //! The energy storage in the fluid phase with index phaseIdx
            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }
         //! The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);
        return storage;
    }

    /*!
     * \brief Calculates the flux for all mass balance equations.
     *
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux compuation
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
        // get upwind weights into local scope
        NumEqVector flux(0.0);

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = conti0EqIdx + phaseIdx*numComponents + compIdx;

                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [phaseIdx, compIdx] (const auto& volVars)
                { return volVars.molarDensity(phaseIdx)*volVars.moleFraction(phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                // do not add diffusive flux of main component, as that is not done in master as well
                if (compIdx == phaseIdx)
                        continue;
                flux[eqIdx] += diffusiveFluxes[compIdx];
            }
            //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;
    }

    /*!
     * \brief Calculates the source term of the equation.
     *
     * \param problem The source term
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scv The sub-control volume over which we integrate the source term
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);
        // In the case of a kinetic consideration, mass transfer
        // between phases is realized via source terms there is a
        // balance equation for each component in each phase
        const auto& volVars = elemVolVars[scv];
        std::array<std::array<Scalar, numComponents>, numPhases> componentIntoPhaseMassTransfer = {{{0.0},{0.0}}};

        //get characteristic length
        const Scalar characteristicLength   = volVars.characteristicLength()  ;
        const Scalar factorMassTransfer     = volVars.factorMassTransfer()  ;

        const Scalar awn = volVars.interfacialArea(phase0Idx, phase1Idx);

        for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const Scalar sherwoodNumber = volVars.sherwoodNumber(phaseIdx);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                if (compIdx <= numPhases)
                {
                    //we only have to do the source calculation for one component going into the other phase as for each component: n_wn -> n_n = - n_n -> n_wn
                    if (phaseIdx == compIdx)
                        continue;

                    const Scalar xNonEquil = volVars.moleFraction(phaseIdx, compIdx);

                    //additionally get equilibrium values from volume variables
                    const Scalar xEquil = volVars.xEquil(phaseIdx, compIdx);
                    //get the diffusion coefficient
                    const Scalar diffCoeff = volVars.diffusionCoefficient(phaseIdx, compIdx);

                    //now compute the flux
                    const Scalar compFluxIntoOtherPhase = factorMassTransfer * (xEquil-xNonEquil)*characteristicLength * awn * volVars.molarDensity(phaseIdx) * diffCoeff * sherwoodNumber;

                    componentIntoPhaseMassTransfer[phaseIdx][compIdx] += compFluxIntoOtherPhase;
                    componentIntoPhaseMassTransfer[compIdx][compIdx] += -compFluxIntoOtherPhase;
                }
            }
        }

        // Actually add the mass transfer to the sources which might
        // exist externally
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                const unsigned int eqIdx = conti0EqIdx + compIdx + phaseIdx*numComponents;
                source[eqIdx] += componentIntoPhaseMassTransfer[phaseIdx][compIdx];

                using std::isfinite;
                if (!isfinite(source[eqIdx]))
                    DUNE_THROW(NumericalProblem, "Calculated non-finite source");
            }
        }

        if (ModelTraits::enableThermalNonEquilibrium())
        {
            // Call the (kinetic) Energy module, for the source term.
            // it has to be called from here, because the mass transfered has to be known.
            EnergyLocalResidual::computeSourceEnergy(source,
                                                     element,
                                                     fvGeometry,
                                                     elemVolVars,
                                                     scv);
        }

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);
        Valgrind::CheckDefined(source);

        return source;
    }
 };

} // end namespace Dumux

#endif
