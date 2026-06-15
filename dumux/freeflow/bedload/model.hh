// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransportModel
 * \brief A bedload transport model.
 */
#ifndef DUMUX_BEDLOAD_MODEL_HH
#define DUMUX_BEDLOAD_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/flux/fluxvariablescaching.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "indices.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup BedloadModel
 * \brief Traits for the bedload transport model.
 */
 template<int nGrainClasses>
struct BedloadModelTraits
{
    using Indices = BedloadIndices;

    static constexpr int numGrainClasses() { return nGrainClasses; }
    static constexpr int numEq() { return nGrainClasses; }
    static constexpr int numPhases() { return 1; }

    // Enable sediment
    static constexpr bool enableBedload() { return true; }
};

/*!
 * \ingroup BedloadModel
 * \brief Traits class for the volume variables of the bedload transport model.
 *
 * \tparam PV The type used for primary variables
 * \tparam MT The model traits
 */
template<class PV,
         class MT>
struct BedloadVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};


namespace Properties {

//! Type tag for the bedload transport model
// Create new type tags
namespace TTag {
struct Bedload{ using InheritsFrom = std::tuple<ModelProperties>; };
}// end namespace TTag

//////////////////////////////////////////////////////////////////
// Define properties
//////////////////////////////////////////////////////////////////

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Bedload>
{ using type = BedloadModelTraits<TypeTag::nGrainClasses>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Bedload>
{ using type = BedloadLocalResidual<TypeTag>; };

template<class TypeTag>
struct FluxVariables<TypeTag, TTag::Bedload>
{ using type = BedloadFluxVariables<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Bedload>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = BedloadVolumeVariablesTraits<PV, MT>;
public:
    using type = BedloadVolumeVariables<Traits>;
};

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::Bedload>
{ using type = FluxVariablesCaching::EmptyCache< GetPropType<TypeTag, Properties::Scalar> >; };

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::Bedload>
{ using type = FluxVariablesCaching::EmptyCacheFiller; };

template<class TypeTag>
struct IOFields<TypeTag, TTag::Bedload>
{ using type = BedloadIOFields; };

} // namespace Properties
} // namespace Dumux

#endif
