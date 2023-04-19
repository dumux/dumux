// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SSTModel
 * \copydoc Dumux::SSTResidualImpl
 */
#ifndef DUMUX_SST_STAGGERED_KOMEGA_LOCAL_RESIDUAL_HH
#define DUMUX_SST_STAGGERED_KOMEGA_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup SSTModel
 * \brief Element-wise calculation of the residual for SST models using the staggered discretization
 */

// forward declaration
template<class TypeTag, class BaseLocalResidual, class DiscretizationMethod>
class SSTResidualImpl;

template<class TypeTag, class BaseLocalResidual>
class SSTResidualImpl<TypeTag, BaseLocalResidual, DiscretizationMethods::Staggered>
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

        if(problem.sstModelVersion() == SSTModel::BSL)
        {
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
        else if(problem.sstModelVersion() == SSTModel::SST)
        {
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
        else
            DUNE_THROW(Dune::NotImplemented, "\nThis SST Model is not implemented.\n");

        return source;
    }
};
} // end namespace Dumux

#endif
