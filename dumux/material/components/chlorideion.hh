// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief A class for the Cl- (Chloride ion) component properties
 */
#ifndef DUMUX_MATERIAL_COMPONENTS_CL_ION_HH
#define DUMUX_MATERIAL_COMPONENTS_CL_ION_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Cl- (Chloride ion) component properties
 */
template <class Scalar>
class ChlorideIon
: public Components::Base<Scalar, ChlorideIon<Scalar> >
, public Components::Ion<Scalar, ChlorideIon<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the Cl- ion.
     */
    static std::string name()
    { return "Cl-"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the Cl- ion.
     */
    static Scalar molarMass()
    { return 35.453e-3; }

    /*!
     * \brief The charge of the Cl- ion.
     */
    static constexpr int charge()
    {
        return -1;
    }

};

} // end namespace Components
} // end namespace Dumux

#endif
