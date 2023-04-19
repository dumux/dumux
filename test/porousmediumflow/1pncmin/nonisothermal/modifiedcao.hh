// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCMinTests
 * \brief Corrected material properties of pure Calcium-Oxide \f$CaO\f$ without considering a porosity
 * change in the reaction of Calciumoxyde and Calciumhydroxyde.
 */

#ifndef DUMUX_MODIFIED_CAO_HH
#define DUMUX_MODIFIED_CAO_HH


#include <dumux/material/components/cao.hh>

namespace Dumux {
namespace Components {
/*!
 * \ingroup OnePNCMinTests
 * \brief A class for the ModifiedCaO properties.
 *
 * This class uses a different CaO density. It is to be used for calculating the
 * chemical reaction of CaO to Ca(OH)2 without considering the porosity change
 * according to Shao et al. (2013) \cite shao2013.
 */
template <class Scalar>
class ModifiedCaO : public  Components::CaO<Scalar>
{
public:

    /*!
     * \brief The corrected mass density \f$\mathrm{[kg/m^3]}\f$ of CaO.
     *
     * This density is to be used for calculating the chemical reaction
     * of CaO to Ca(OH)2 without considering the solid volume change.
     * See Shao et al. (2013) \cite shao2013.
     */
    static Scalar solidDensity(Scalar temperature)
    {
        return 1656;
    }

};
} // end namespace Components
} // end namespace Dumux

#endif
