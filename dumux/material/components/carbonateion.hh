// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief A class for the CO3 ion properties.
 */
#ifndef DUMUX_CO3_ION_HH
#define DUMUX_CO3_ION_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the CO3 fluid properties.
 */
template <class Scalar>
class CarbonateIon
: public Components::Base<Scalar, CarbonateIon<Scalar> >
, public Components::Ion<Scalar, CarbonateIon<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the CO3 ion.
     */
    static std::string name()
    { return "CO3-"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the CO3 ion.
     */
    static constexpr Scalar molarMass()
    { return 60.0092e-3; } // kg/mol

    /*!
     * \brief The charge balance of the CO3 ion.
     */
    static constexpr int charge()
    { return -2; }

};

} // end namespace Components
} // end namespace Dumux

#endif
