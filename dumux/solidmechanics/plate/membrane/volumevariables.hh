// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MembranePlateModel
 * \brief Volume variables for the membrane plate model
 */
#ifndef DUMUX_MEMBRANE_PLATE_VOLUME_VARIABLES_HH
#define DUMUX_MEMBRANE_PLATE_VOLUME_VARIABLES_HH

#include <dumux/common/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup MembranePlateModel
 * \brief Volume variables for the membrane plate model
 */
template<class Traits>
class MembranePlateVolumeVariables
: public BasicVolumeVariables<Traits>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    //! Returns the vertical deformation \f$ w \f$
    Scalar deformation() const
    { return this->priVar(Indices::deformationIdx); }
};

} // end namespace Dumux

#endif
