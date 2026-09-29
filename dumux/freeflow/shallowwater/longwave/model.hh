// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 *
 * \brief Long-wave approximations of the shallow water equations
 *
 * For long waves in shallow water, the inertial terms of the momentum balance are small
 * compared to gravity and friction. Dropping them, the friction slope equals the slope of the
 * free surface \f$ H = h + z \f$, and the momentum balance reduces to a friction law that
 * determines the discharge per unit width \f$ \mathbf{q} \f$ from the water depth \f$ h \f$.
 * With Manning's law, the remaining mass balance reads
 *
 * \f[
 * \frac{\partial h}{\partial t} + \nabla \cdot \mathbf{q} = r, \qquad
 * \mathbf{q} = -\frac{h R^{2/3}}{n} \frac{\mathbf{s}}{\sqrt{|\mathbf{s}|}}, \qquad
 * \mathbf{s} = \nabla z + w \nabla h,
 * \f]
 *
 * with the bed surface elevation \f$ z \f$ (in \f$ m \f$), Manning's coefficient \f$ n \f$ (in \f$ s m^{-1/3} \f$),
 * the hydraulic radius \f$ R = h \f$ of unconfined sheet flow, and a source term \f$ r \f$
 * (in \f$ m s^{-1} \f$), e.g. rainfall. The weight \f$ w \f$ selects the approximation
 * (parameter `LongWave.Approximation`, see Dumux::LongWave::WaveApproximation):
 *
 * - `diffusive` (default), \f$ w = 1 \f$: the diffusive wave, driven by the free-surface gradient,
 * - `kinematic`, \f$ w = 0 \f$: the kinematic wave, driven by the bed slope alone,
 * - `inertiacorrected`, \f$ w = \max(0, 1 - \mathsf{V}^2) \f$: the diffusive wave with the
 *   diffusivity of the linearized shallow water equations about uniform flow, where
 *   \f$ \mathsf{V} \f$ is the Vedernikov number of normal flow.
 *
 * The flux coefficient \f$ D = h R^{2/3}/(n\sqrt{|\mathbf{s}|}) \f$ diverges as the free
 * surface levels and its derivative vanishes as the depth goes to zero. Both limits are
 * regularized, see Dumux::LongWave::regularInvSqrt, Dumux::LongWave::gradHThreshold and
 * Dumux::LongWave::regularizedDepth.
 *
 * The model is implemented for control-volume finite element schemes (e.g. the box method)
 * in one and two dimensions. In one dimension, the extrusion factor is the width of the channel.
 * The spatial parameters have to provide `bedSurface(element, scv)`, `manningN(element)` and
 * `extrusionFactor(element, scv, elemSol)`.
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_MODEL_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "indices.hh"
#include "iofields.hh"
#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

/*!
 * \ingroup LongWaveModel
 * \brief Specifies a number of properties of the long-wave models
 */
struct LongWaveModelTraits
{
    using Indices = LongWaveIndices;

    static constexpr int numEq() { return 1; }
};

/*!
 * \ingroup LongWaveModel
 * \brief Traits class for the volume variables of the long-wave models
 *
 * \tparam PV The type used for primary variables
 * \tparam MT The model traits
 */
template<class PV, class MT>
struct LongWaveVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

namespace Properties {

//! Type tag for the long-wave models
namespace TTag {
struct LongWave { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::LongWave>
{ using type = LongWaveModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::LongWave>
{ using type = LongWaveLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::LongWave>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
public:
    using type = LongWaveVolumeVariables<LongWaveVolumeVariablesTraits<PV, MT>>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::LongWave>
{ using type = LongWaveIOFields; };

} // end namespace Properties
} // end namespace Dumux

#endif
