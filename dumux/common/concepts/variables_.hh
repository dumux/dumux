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

template<class>
inline constexpr bool AlwaysFalse = false;

template<class T>
concept FVGridVariables =
    requires(T t)
    {
        typename T::GridVolumeVariables;
        typename T::VolumeVariables;
        typename T::GridFluxVariablesCache;
        t.curGridVolVars();
        t.prevGridVolVars();
    };

template<class T>
concept GridVariables =
    !FVGridVariables<T>
    && requires(T t)
    {
        typename T::GridVariablesCache;
        typename T::Variables;
        t.curGridVars();
        t.prevGridVars();
    };

template<class T>
struct VariablesType
{
    static_assert(AlwaysFalse<T>, "Type does not satisfy Concept::FVGridVariables or Concept::GridVariables");
};

template<FVGridVariables T>
struct VariablesType<T>
{
    using type = typename T::VolumeVariables;
};

template<GridVariables T>
struct VariablesType<T>
{
    using type = typename T::Variables;
};

template<class T>
using Variables_t = typename VariablesType<T>::type;

template<class T>
struct GridVariablesCacheType
{
    static_assert(AlwaysFalse<T>, "Type does not satisfy Concept::FVGridVariables or Concept::GridVariables");
};

template<FVGridVariables T>
struct GridVariablesCacheType<T>
{
    using type = typename T::GridVolumeVariables;
};

template<GridVariables T>
struct GridVariablesCacheType<T>
{
    using type = typename T::GridVariablesCache;
};

template<class T>
using GridVariablesCache_t = typename GridVariablesCacheType<T>::type;

} // end namespace Dumux::Concept

#endif
