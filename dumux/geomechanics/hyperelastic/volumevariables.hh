// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Volume variables for the hyperelasticity model
 */

#ifndef DUMUX_GEOMECHANICS_HYPERELASTIC_VOLUME_VARIABLES_HH
#define DUMUX_GEOMECHANICS_HYPERELASTIC_VOLUME_VARIABLES_HH

#include <dumux/common/volumevariables.hh>

namespace Dumux {
/*!
 * \ingroup Hyperelastic
 * \brief Volume variables for the hyperelasticity model
 */
template <class Traits>
class HyperelasticVolumeVariables
: public BasicVolumeVariables<Traits>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;

    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

    Scalar displacement(int i) const
    { return this->priVar(i); }

    const PrimaryVariables& displacement() const
    { return this->priVars(); }
};

} // end namespace Dumux

#endif
