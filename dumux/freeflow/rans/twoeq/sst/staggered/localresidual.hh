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
 * \ingroup SSTModel
 * \copydoc Dumux::SSTResidualImpl
 */
#ifndef DUMUX_SST_STAGGERED_KOMEGA_LOCAL_RESIDUAL_HH
#define DUMUX_SST_STAGGERED_KOMEGA_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup SSTModel
 * \brief Element-wise calculation of the residual for SST models using the staggered discretization
 */

// forward declaration
template<class TypeTag, class BaseLocalResidual, DiscretizationMethod discMethod>
class SSTResidualImpl;

template<class TypeTag, class BaseLocalResidual>
class SSTResidualImpl<TypeTag, BaseLocalResidual, DiscretizationMethod::staggered>
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

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using CellCenterResidual = CellCenterPrimaryVariables;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr int turbulentKineticEnergyEqIdx = Indices::turbulentKineticEnergyEqIdx - ModelTraits::dim();
    static constexpr int dissipationEqIdx = Indices::dissipationEqIdx - ModelTraits::dim();

public:
    using ParentType::ParentType;

    //! Evaluate fluxes entering or leaving the cell center control volume.
    CellCenterPrimaryVariables computeStorageForCellCenter(const Problem& problem,
                                                           const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage = ParentType::computeStorageForCellCenter(problem, scv, volVars);

        storage[turbulentKineticEnergyEqIdx] = volVars.turbulentKineticEnergy()*volVars.density();
        storage[dissipationEqIdx] = volVars.dissipation()*volVars.density();

        return storage;
    }

    CellCenterPrimaryVariables computeSourceForCellCenter(const Problem& problem,
                                                          const Element &element,
                                                          const FVElementGeometry& fvGeometry,
                                                          const ElementVolumeVariables& elemVolVars,
                                                          const ElementFaceVariables& elemFaceVars,
                                                          const SubControlVolume &scv) const
    {
        CellCenterPrimaryVariables source = ParentType::computeSourceForCellCenter(problem, element, fvGeometry,
                                                                                   elemVolVars, elemFaceVars, scv);

        using std::min;
        const auto& volVars = elemVolVars[scv];

        const std::string turbulenceModel = problem.readInTurbulenceModel();
        if(turbulenceModel == "BSL"){
            // production
            Scalar productionTerm = 2.0 * volVars.dynamicEddyViscosity() * volVars.stressTensorScalarProduct();
            source[turbulentKineticEnergyEqIdx] += productionTerm;
            source[dissipationEqIdx] +=  volVars.gammaBSL() / volVars.kinematicEddyViscosity() * productionTerm;

            // destruction
            source[turbulentKineticEnergyEqIdx] -= volVars.betaStarBSL() * volVars.density() * volVars.turbulentKineticEnergy() * volVars.dissipation();
            source[dissipationEqIdx] -= volVars.betaBSL() * volVars.density() * volVars.dissipation() * volVars.dissipation();

            // cross-diffusion term
            Scalar gradientProduct = 0.0;
            for (unsigned int i = 0; i < ModelTraits::dim(); ++i)
                gradientProduct += volVars.storedTurbulentKineticEnergyGradient()[i]
                                * volVars.storedDissipationGradient()[i];

            source[dissipationEqIdx] += 2.0 * volVars.density() * (1.0-volVars.F1()) * volVars.sigmaOmega2() / volVars.dissipation() * gradientProduct;
        }
        else if(turbulenceModel == "SST"){
            // production
            Scalar productionTerm = 2.0 * volVars.dynamicEddyViscosity() * volVars.stressTensorScalarProduct();
            source[turbulentKineticEnergyEqIdx] += productionTerm;
            source[dissipationEqIdx] += volVars.gammaSST() / volVars.kinematicEddyViscosity() * productionTerm;

            // destruction
            source[turbulentKineticEnergyEqIdx] -= volVars.betaStarSST() * volVars.density() * volVars.turbulentKineticEnergy() * volVars.dissipation();
            source[dissipationEqIdx] -= volVars.betaSST() * volVars.density() * volVars.dissipation() * volVars.dissipation();

            // cross-diffusion term
            Scalar gradientProduct = 0.0;
            for (unsigned int i = 0; i < ModelTraits::dim(); ++i)
                gradientProduct += volVars.storedTurbulentKineticEnergyGradient()[i]
                                * volVars.storedDissipationGradient()[i];
            source[dissipationEqIdx] += 2.0 * volVars.density() * (1.0-volVars.F1()) * volVars.sigmaOmega2() / volVars.dissipation() * gradientProduct;
        }

        return source;
    }
};
} // end namespace Dumux

#endif
