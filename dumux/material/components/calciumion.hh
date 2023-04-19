// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief A class for the Ca2+ (Calcium ion) component properties
 */
#ifndef DUMUX_CA_ION_HH
#define DUMUX_CA_ION_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Ca2+ (Calcium ion) component properties.
 */
template <class Scalar>
class CalciumIon
: public Components::Base<Scalar, CalciumIon<Scalar> >
, public Components::Ion<Scalar, CalciumIon<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the Ca ion.
     */
    static std::string name()
    { return "Ca2+"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the Ca ion.
     */
    static constexpr Scalar molarMass()
    { return 40.078e-3; } // kg/mol

    /*!
     * \brief The charge balance of the Ca ion.
     */
    static constexpr int charge()
    { return +2; }

};

} // end namespace Components
} // end namespace Dumux

#endif
