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
 * \ingroup PorousmediumflowModels
 * \brief Sub-control entity-local evaluation of the operators
 *        within in the systems of equations of n-phase immiscible models.
 */
#ifndef DUMUX_FV_IMMISCIBLE_OPERATORS_HH
#define DUMUX_FV_IMMISCIBLE_OPERATORS_HH

#include <dumux/assembly/fv/operators.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/porousmediumflow/nonisothermal/operators.hh>

namespace Dumux::Experimental {
namespace Detail {

template<class ModelTraits, class ElemVars,
         bool enable = ModelTraits::enableEnergyBalance()>
struct DefaultEnergyOperatorsChooser;

template<class ModelTraits, class ElemVars>
struct DefaultEnergyOperatorsChooser<ModelTraits, ElemVars, true>
{ using type = FVNonIsothermalOperators<ModelTraits, FVOperatorsContext<ElemVars>>; };

template<class ModelTraits, class ElemVars>
struct DefaultEnergyOperatorsChooser<ModelTraits, ElemVars, false>
{ using type = void; };

} // end namespace Detail

/*!
 * \ingroup PorousmediumflowModels
 * \brief Sub-control entity-local evaluation of the operators
 *        within in the systems of equations of n-phase immiscible models.
 * \tparam ModelTraits defines model-related types and variables (e.g. number of phases)
 * \tparam FluxVariables the type that is responsible for computing the individual
 *                       flux contributions, i.e., advective, diffusive, convective...
 * \tparam ElemVars The element-(stencil)-local variables
 * \tparam EnergyOperators optional template argument, specifying the class that
 *                         handles the operators related to non-isothermal effects.
 */
template<class ModelTraits,
         class FluxVariables,
         class ElemVars,
         class EnergyOperators = typename Detail::DefaultEnergyOperatorsChooser<ModelTraits, ElemVars>::type>
class FVImmiscibleOperators
: public FVOperators<ElemVars>
{
    using ParentType = FVOperators<ElemVars>;
    using LC = typename ParentType::LocalContext;

    // The variables required for the evaluation of the equation
    using ElementVariables = typename LC::ElementVariables;
    using GridVariables = typename ElementVariables::GridVariables;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using ElementVolumeVariables = typename ElementVariables::ElementVolumeVariables;
    using ElementFluxVariablesCache = typename ElementVariables::ElementFluxVariablesCache;

    // The grid geometry on which the scheme operates
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using Problem = std::decay_t<decltype(std::declval<GridVariables>().gridVolVars().problem())>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int conti0EqIdx = ModelTraits::Indices::conti0EqIdx;
    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();;

public:
    //! export the type of context on which this class operates
    using LocalContext = LC;

    //! export the type used to store scalar values for all equations
    using typename ParentType::NumEqVector;

    // export the types of the flux/source/storage terms
    using typename ParentType::FluxTerm;
    using typename ParentType::SourceTerm;
    using typename ParentType::StorageTerm;

    /*!
     * \brief Compute the storage term of the equations for the given sub-control volume
     * \param problem The problem to be solved (could store additionally required quantities)
     * \param context The element-stencil-local required data
     * \param scv The sub-control volume
     */
     template<class Problem>
     static StorageTerm storage(const Problem& problem,
                                const LocalContext& context,
                                const SubControlVolume& scv)
    {
        const auto& volVars = context.elementVariables().elemVolVars()[scv];

        StorageTerm storage;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            auto eqIdx = conti0EqIdx + phaseIdx;
            storage[eqIdx] = volVars.porosity()
                             * volVars.density(phaseIdx)
                             * volVars.saturation(phaseIdx);

            // The energy storage in the fluid phase with index phaseIdx
            if constexpr (enableEnergyBalance)
                EnergyOperators::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }

        // The energy storage in the solid matrix
        if constexpr (enableEnergyBalance)
            EnergyOperators::solidPhaseStorage(storage, scv, volVars);

        // multiply with volume
        storage *= Extrusion::volume(scv)*volVars.extrusionFactor();

        return storage;
    }

    /*!
     * \brief Compute the flux term of the equations for the given sub-control volume face
     * \param problem The problem to be solved (could store additionally required quantities)
     * \param context The element-stencil-local required data
     * \param scvf The sub-control volume face for which the flux term is to be computed
     * \note This must be overloaded by the implementation
     */
    template<class Problem>
    static FluxTerm flux(const Problem& problem,
                         const LocalContext& context,
                         const SubControlVolumeFace& scvf)
    {
        FluxVariables fluxVars;
        fluxVars.init(context, scvf);

        NumEqVector flux;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // the physical quantities for which we perform upwinding
            auto upwindTerm = [phaseIdx](const auto& volVars)
                              { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx); };

            auto eqIdx = conti0EqIdx + phaseIdx;
            flux[eqIdx] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);

            // Add advective phase energy fluxes. For isothermal model the contribution is zero.
            if constexpr (enableEnergyBalance)
                EnergyOperators::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        // Add diffusive energy fluxes. For isothermal model the contribution is zero.
        if constexpr (enableEnergyBalance)
            EnergyOperators::heatConductionFlux(flux, fluxVars);

        return flux;
    }
};

} // end namespace Dumux::Experimental

#endif
