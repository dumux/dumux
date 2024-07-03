// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Interface for projectors between subdomains and mortars.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_PROJECTOR_INTERFACE_HH
#define DUMUX_MULTIDOMAIN_MORTAR_PROJECTOR_INTERFACE_HH

namespace Dumux::Mortar {

/*!
 * \ingroup MultiDomain
 * \brief Interface for projectors between subdomains and mortars.
 */
template<typename S>
class ProjectorInterface
{
public:
    using SolutionVector = S;

    //! project the vector x into the target domain
    virtual SolutionVector project(const SolutionVector& x) const = 0;
};

} // end namespace Dumux::Mortar

#endif
