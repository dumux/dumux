// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concepts related to integration point data
 */
#ifndef DUMUX_CONCEPTS_IPDATA__HH
#define DUMUX_CONCEPTS_IPDATA__HH

#include <concepts>

namespace Dumux::Concept {

template<class T>
concept IpData = requires(T a)
{
    a.ipGlobal();
    a.ipLocal();
};

template<class T>
concept LocalDofIndexProvider = requires(T a)
{
    a.localDofIndex();
};

template<class T>
concept LocalDofIpData = IpData<T> && LocalDofIndexProvider<T>;

} // end namespace Dumux::Concept

#endif
