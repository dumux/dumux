// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief Nonisothermal mass-domain model for the zero-equation RANS turbulence closure.
 *
 * The isothermal zero-equation closure (dumux/freeflow/rans/zeroeq/problem.hh) adds its eddy
 * viscosity purely on the momentum side and needs no dedicated mass-domain TypeTag. The
 * nonisothermal variant needs the mass domain to know the momentum-side eddy viscosity in order
 * to compute an eddy thermal conductivity. The TypeTag `NavierStokesMassOneZeroEqNI` therefore
 * inherits from `NavierStokesMassOnePNI` and only specializes `VolumeVariables` and
 * `ThermalConductivityModel`. As no additional transport equations are solved, everything else
 * is inherited unchanged.
 */
#ifndef DUMUX_RANS_ZEROEQ_MASS_MODEL_HH
#define DUMUX_RANS_ZEROEQ_MASS_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/rans/common/thermalconductivitymodel.hh>

#include "volumevariables.hh"

namespace Dumux {
namespace Properties {

namespace TTag {
//! The type tag for the single-phase, non-isothermal Navier-Stokes mass model with the
//! eddy viscosity of the zero-equation turbulence closure.
struct NavierStokesMassOneZeroEqNI { using InheritsFrom = std::tuple<NavierStokesMassOnePNI>; };
} // end namespace TTag

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOneZeroEqNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using BaseTraits = NavierStokesMassOnePVolumeVariablesTraits<PV, FSY, FST, MT>;
    using ETCM = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using HCT = GetPropType<TypeTag, Properties::HeatConductionType>;
    struct NITraits : public BaseTraits
    {
        using EffectiveThermalConductivityModel = ETCM;
        using HeatConductionType = HCT;
    };
public:
    using type = ZeroEqMassVolumeVariables<NITraits>;
};

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::NavierStokesMassOneZeroEqNI>
{ using type = RANSThermalConductivityModel; };

} // end namespace Properties
} // end namespace Dumux

#endif
