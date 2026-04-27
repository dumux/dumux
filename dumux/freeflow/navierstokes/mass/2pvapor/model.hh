// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_2PVAPOR_MODEL_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_2PVAPOR_MODEL_HH

#include <dumux/freeflow/navierstokes/mass/2p/model.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "localresidual.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \brief Model traits for the 2p Cahn-Hilliard model extended with vapor transport.
 *
 * Adds one equation (vapor convection-diffusion) to the base 2p model.
 * Primary variables: [p, φ, μ, c_v]
 */
struct NavierStokesMassTwoPVaporModelTraits : public NavierStokesMassTwoPModelTraits
{
    static constexpr int numEq() { return 4; }
    using Indices = NavierStokesMassTwoPVaporIndices;
};

namespace Properties {

namespace TTag {
struct NavierStokesMassTwoPVapor
{ using InheritsFrom = std::tuple<NavierStokesMassTwoP>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassTwoPVapor>
{ using type = NavierStokesMassTwoPVaporModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassTwoPVapor>
{ using type = NavierStokesMassTwoPVaporLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassTwoPVapor>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = NavierStokesMassTwoPVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = NavierStokesMassTwoPVaporVolumeVariables<Traits>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassTwoPVapor>
{ using type = NavierStokesMassTwoPVaporIOFields; };

} // end namespace Properties
} // end namespace Dumux

#endif
