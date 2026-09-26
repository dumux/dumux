// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FoepplVonKarmanPlate
 * \brief Volume variables for the Föppl-von Kármán model
 */

#ifndef DUMUX_FOEPPL_VON_KARMAN_PLATE_VOLUME_VARIABLES_HH
#define DUMUX_FOEPPL_VON_KARMAN_PLATE_VOLUME_VARIABLES_HH

#include <dumux/common/volumevariables.hh>
#include <dumux/solidmechanics/plate/kirchhoff_love/volumevariables.hh>

namespace Dumux {

template<class Traits>
using FoepplVonKarmanPlateDeformationVolumeVariables
    = KirchhoffLovePlateDeformationVolumeVariables<Traits>;

template<class Traits>
using FoepplVonKarmanPlateRotationVolumeVariables
    = KirchhoffLovePlateRotationVolumeVariables<Traits>;

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief Volume variables for the in-plane displacements
 */
template<class Traits>
class FoepplVonKarmanPlateInPlaneVolumeVariables
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
};

} // end namespace Dumux

#endif
