// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 * \copydoc Dumux::LongWaveLocalResidual
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_LOCALRESIDUAL_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_LOCALRESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

#include "discharge.hh"

namespace Dumux {

/*!
 * \ingroup LongWaveModel
 * \brief Element-wise calculation of the residual for the long-wave approximations of the
 *        shallow water equations
 *
 * On grid geometries with interior boundaries (facet coupling), the flux over faces on an
 * interior boundary is provided by the coupling manager of the problem through
 * `computeBulkFlux`.
 */
template<class TypeTag>
class LongWaveLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
        "The long-wave local residual is implemented for control-volume finite element schemes");

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluate the storage term, the water depth \f$ h \f$
     *
     * Since the bed does not move, \f$ \partial_t h = \partial_t (h + z) \f$.
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage;
        storage[Indices::massBalanceIdx] = volVars.waterDepth();
        return storage;
    }

    /*!
     * \brief Evaluate the volumetric discharge over a sub-control-volume face
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        if constexpr (requires { scvf.interiorBoundary(); })
            if (scvf.interiorBoundary())
                return problem.couplingManager().computeBulkFlux(
                    problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache[scvf]
                );

        NumEqVector flux(0.0);
        flux[Indices::massBalanceIdx] = LongWave::discharge(
            problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache
        );
        return flux;
    }
};

} // end namespace Dumux

#endif
