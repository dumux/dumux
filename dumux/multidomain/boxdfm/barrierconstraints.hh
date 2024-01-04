// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup BoxDFMModel
 * \brief Predefined constraint implementations for the box-dfm model with support
 *        for barriers. Each implementation is physics-specific.
 */

#ifndef DUMUX_MULTIDOMAIN_BOXDFM_BARRIER_CONSTRAINTS_HH
#define DUMUX_MULTIDOMAIN_BOXDFM_BARRIER_CONSTRAINTS_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

namespace Dumux {

namespace Properties {
namespace TTag { struct BoxDfmImmisibleConstraintModel { using InheritsFrom = std::tuple<ModelProperties>; }; }

template<class TypeTag>
struct IOFields<TypeTag, TTag::BoxDfmImmisibleConstraintModel> {
    using type = OnePIOFields;
};

} // namespace Properties

/*!
 * \ingroup MultiDomain
 * \ingroup BoxDFMModel
 * \brief Implementation of the box-dfm constraint on barriers for immiscible pm-flow models.
 */
template<...>
class BoxDfmImmiscibleBarrierConstraint
{

};

} // end namespace Dumux

#endif
