// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Interface shared by all time integration methods
 */
#ifndef DUMUX_TIMESTEPPING_TIME_INTEGRATION_METHOD_HH
#define DUMUX_TIMESTEPPING_TIME_INTEGRATION_METHOD_HH

#include <string>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \brief Interface shared by all time integration methods
 *
 * Assemblers only need to know whether the spatial operator enters the
 * implicit solves of a method, because this determines the sparsity pattern
 * of the Jacobian and how its entries are computed.
 */
template<class Scalar>
class TimeIntegrationMethod
{
public:
    //! Whether the spatial operator contributes to the Jacobian of the implicit solves
    virtual bool implicit () const = 0;

    //! The name of the method
    virtual std::string name () const = 0;

    virtual ~TimeIntegrationMethod() = default;
};

} // end namespace Dumux::Experimental

#endif
