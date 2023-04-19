// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup Components
 * \brief A class for the Glucose component properties
 */
#ifndef DUMUX_MATERIAL_COMPONENTS_GLUCOSE_HH
#define DUMUX_MATERIAL_COMPONENTS_GLUCOSE_HH

#include <dumux/material/components/base.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Glucose component properties
 */
template <class Scalar>
class Glucose
: public Components::Base<Scalar, Glucose<Scalar> >
{
public:

    /*!
     * \brief A human readable name for Glucose.
     */
    static std::string name()
    { return "Glucose"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Glucose.
     */
    static Scalar molarMass()
    { return 0.18016; } // kg/mol
};

} // end namespace Components
} // end namespace Dumux

#endif
