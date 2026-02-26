// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concepts related to grid variables cache interfaces
 */
#ifndef DUMUX_CONCEPTS_VARIABLES__HH
#define DUMUX_CONCEPTS_VARIABLES__HH

namespace Dumux::Concept {

template<class T>
concept HasGridVariablesCache = requires
{
    typename T::GridVariablesCache;
};

template<class T>
concept HasGridVolumeVariables = requires
{
    typename T::GridVolumeVariables;
};

template<class T,
         bool hasGridVariablesCache = HasGridVariablesCache<T>,
         bool hasGridVolumeVariables = HasGridVolumeVariables<T>>
struct GridVariablesCacheType;

template<class T, bool hasGridVolumeVariables>
struct GridVariablesCacheType<T, true, hasGridVolumeVariables>
{
    using type = typename T::GridVariablesCache;
};

template<class T>
struct GridVariablesCacheType<T, false, true>
{
    using type = typename T::GridVolumeVariables;
};

template<class T>
using GridVariablesCache_t = typename GridVariablesCacheType<T>::type;

template<class T>
concept HasVariables = requires
{
    typename T::Variables;
};

template<class T>
concept HasVolumeVariables = requires
{
    typename T::VolumeVariables;
};

template<class T,
         bool hasVariables = HasVariables<T>,
         bool hasVolumeVariables = HasVolumeVariables<T>>
struct VariablesType;

template<class T, bool hasVolumeVariables>
struct VariablesType<T, true, hasVolumeVariables>
{
    using type = typename T::Variables;
};

template<class T>
struct VariablesType<T, false, true>
{
    using type = typename T::VolumeVariables;
};

template<class T>
using Variables_t = typename VariablesType<T>::type;

} // end namespace Dumux::Concept

#endif
