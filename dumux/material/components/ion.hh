// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Interface for components that are ions.
 */
#ifndef DUMUX_COMPONENT_ION_HH
#define DUMUX_COMPONENT_ION_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Interface for components that are ions.
 */
template<class Scalar, class Component>
class Ion
{
public:
    /*!
     * \brief Returns the charge of the ion.
     */
    template<class C = Component>
    static constexpr int charge()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: charge()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }
};

} // end namespace Components
} // end namespace Dumux

#endif
