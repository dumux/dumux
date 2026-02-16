// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concepts related to interpolation point data
 */
#ifndef DUMUX_CONCEPTS_IPDATA__HH
#define DUMUX_CONCEPTS_IPDATA__HH

#include <concepts>

namespace Dumux::Concept::Detail{

    template<class T>
    concept LocalDofIndexProvider = requires(T a)
    {
        a.localDofIndex();
    };

    template<class T>
    concept ScvfIndexProvider = requires(T a)
    {
        a.scvfIndex();
    };

    template<class T>
    concept QpIndexProvider = requires(T a)
    {
        a.qpIndex();
    };

    template<class T>
    concept IntersectionIndexProvider = requires(T a)
    {
        a.intersectionIndex();
    };

} // end namespace Dumux::Concept::Detail

namespace Dumux::Concept {

template<class T>
concept IpData = requires(T a)
{
    a.global();
    a.local();
};

template<class T>
concept LocalDofIpData = IpData<T> && Detail::LocalDofIndexProvider<T>;

template<class T>
concept ScvfIpData = IpData<T> && Detail::ScvfIndexProvider<T>;

template<class T>
concept ScvfQpIpData = ScvfIpData<T> && Detail::QpIndexProvider<T>;

template<class T>
concept QIpData = IpData<T> && Detail::QpIndexProvider<T>;

template<class T>
concept IntersectionIpData = IpData<T> && Detail::IntersectionIndexProvider<T>;

template<class T>
concept IntersectionQpIpData = IntersectionIpData<T> && Detail::QpIndexProvider<T>;

} // end namespace Dumux::Concept

#endif
