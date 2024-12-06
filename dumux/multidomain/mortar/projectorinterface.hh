// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Interface for projectors between a mortar and an adjacent subdomain trace.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_PROJECTOR_INTERFACE_HH
#define DUMUX_MULTIDOMAIN_MORTAR_PROJECTOR_INTERFACE_HH

namespace Dumux::Mortar {

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Interface for projectors between a subdomain trace and a mortar.
 */
template<typename SolutionVector>
struct Projector
{
    virtual ~Projector() = default;

    //! Project a mortar solution to a subdomain trace
    virtual SolutionVector toTrace(const SolutionVector& x) const { return toTrace_(x); }
    //! Project a subdomain trace to a mortar
    virtual SolutionVector fromTrace(const SolutionVector& x) const { return fromTrace_(x); }

 private:
    virtual SolutionVector toTrace_(const SolutionVector&) const = 0;
    virtual SolutionVector fromTrace_(const SolutionVector&) const = 0;
};

} // end namespace Dumux::Mortar

#endif
