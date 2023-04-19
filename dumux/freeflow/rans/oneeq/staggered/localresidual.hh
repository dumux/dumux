// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OneEqModel
 * \copydoc Dumux::OneEqResidualImpl
 */
#ifndef DUMUX_STAGGERED_ONEEQ_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_ONEEQ_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup OneEqModel
 * \brief Element-wise calculation of the residual for one-equation turbulence models
 *        using the staggered discretization
 */

// forward declaration
template<class TypeTag, class BaseLocalResidual, class DiscretizationMethod>
class OneEqResidualImpl;

template<class TypeTag, class BaseLocalResidual>
class OneEqResidualImpl<TypeTag, BaseLocalResidual, DiscretizationMethods::Staggered>
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
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr int viscosityTildeEqIdx = Indices::viscosityTildeEqIdx - ModelTraits::dim();

public:
    using ParentType::ParentType;

    //! Evaluate fluxes entering or leaving the cell center control volume.
    CellCenterPrimaryVariables computeStorageForCellCenter(const Problem& problem,
                                                           const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage = ParentType::computeStorageForCellCenter(problem, scv, volVars);
        storage[viscosityTildeEqIdx] = volVars.viscosityTilde() * volVars.density();
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

        const auto& volVars = elemVolVars[scv];

        source[viscosityTildeEqIdx] += volVars.cb1() * (1.0 - volVars.ft2())
                                       * volVars.stressTensorScalarProductTilde()
                                       * volVars.viscosityTilde() * volVars.density();

        source[viscosityTildeEqIdx] -= (volVars.cw1() * volVars.fW()
                                        - volVars.cb1() * volVars.ft2() / problem.karmanConstant() / problem.karmanConstant())
                                       * volVars.viscosityTilde() * volVars.viscosityTilde()
                                       / volVars.wallDistance() / volVars.wallDistance() * volVars.density();;

        for (unsigned int axisIdx = 0; axisIdx < ModelTraits::dim(); ++axisIdx)
        {
            source[viscosityTildeEqIdx] += volVars.cb2() / volVars.sigma()
                                           * volVars.storedViscosityTildeGradient()[axisIdx]
                                           * volVars.storedViscosityTildeGradient()[axisIdx]
                                           * volVars.density();
        }

        return source;
    }
};
} // end namespace Dumux

#endif
