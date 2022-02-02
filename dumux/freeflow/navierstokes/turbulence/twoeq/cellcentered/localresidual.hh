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
 * \ingroup TwoEqModel
 * \copydoc Dumux::TwoEqResidualImpl
 */
#ifndef DUMUX_CC_TWOEQ_LOCAL_RESIDUAL_HH
#define DUMUX_CC_TWOEQ_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/mass/1p/localresidual.hh>
#include <dumux/freeflow/navierstokes/turbulence/twoeq/cellcentered/sources.hh>

namespace Dumux {

/*!
 * \ingroup TwoEqModel
 * \brief Element-wise calculation of the residual for twoeq models using the cctpfa discretization
 */

// forward declaration
template<class TypeTag, class BaseLocalResidual, class DiscretizationMethod>
class TwoEqResidualImpl;

template<class TypeTag, class BaseLocalResidual>
class TwoEqResidualImpl<TypeTag, BaseLocalResidual, DiscretizationMethods::CCTpfa>
: public BaseLocalResidual
{
    using ParentType = BaseLocalResidual;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr int turbulentKineticEnergyEqIdx = Indices::turbulentKineticEnergyEqIdx;
    static constexpr int dissipationEqIdx = Indices::dissipationEqIdx;

public:
    using ParentType::ParentType;

    //! Evaluate fluxes entering or leaving the cell center control volume.
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage = ParentType::computeStorage(problem, scv, volVars);

        storage[turbulentKineticEnergyEqIdx] = volVars.turbulentKineticEnergy() * volVars.density();
        storage[dissipationEqIdx] = volVars.dissipation() * volVars.density();

        return storage;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        if (problem.twoEqTurbulenceModelName() == "Wilcox")
        {
            source[turbulentKineticEnergyEqIdx] += TwoEqSources::wilcoxTKESource(problem, element, fvGeometry, elemVolVars, scv);
            source[dissipationEqIdx] += TwoEqSources::wilcoxDissipationSource(problem, element, fvGeometry, elemVolVars, scv);
        }
        else if (problem.twoEqTurbulenceModelName() == "SST" || problem.twoEqTurbulenceModelName() == "BSL")
        {
            source[turbulentKineticEnergyEqIdx] += TwoEqSources::shearStressTransportTKESource(problem, element, fvGeometry, elemVolVars, scv);
            source[dissipationEqIdx] += TwoEqSources::shearStressTransportDissipationSource(problem, element, fvGeometry, elemVolVars, scv);
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "This turbulence model " << problem.twoEqTurbulenceModelName() <<
                                             " is not available, try a different two-eq model.");
        }

        return source;
    }

    /*!
     * \brief Evaluatex the mass or mole flux over a face of a sub control volume.
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
        return fluxVars.flux(0);
    }

};
} // end namespace Dumux

#endif
