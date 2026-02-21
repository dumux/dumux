// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concepts related to grid data cache interfaces
 */
#ifndef DUMUX_CONCEPTS_DATACACHE__HH
#define DUMUX_CONCEPTS_DATACACHE__HH

namespace Dumux::Concept {

template<class T>
concept HasFluxVariablesCache = requires
{
    typename T::FluxVariablesCache;
};

template<class T>
concept HasDataCache = requires
{
    typename T::DataCache;
};

template<class T>
concept HasGridDataCacheInterface = HasFluxVariablesCache<T> || HasDataCache<T>;

template<class T>
concept HasGridFluxVariablesCache = requires
{
    typename T::GridFluxVariablesCache;
};

template<class T>
concept HasGridDataCache = requires
{
    typename T::GridDataCache;
};

template<class T>
concept HasGridVariablesCacheInterface = HasGridFluxVariablesCache<T> || HasGridDataCache<T>;

template<class T,
         bool hasDataCache = HasDataCache<T>,
         bool hasFluxVariablesCache = HasFluxVariablesCache<T>>
struct DataCacheType;

template<class T, bool hasFluxVariablesCache>
struct DataCacheType<T, true, hasFluxVariablesCache>
{
    using type = typename T::DataCache;
};

template<class T>
struct DataCacheType<T, false, true>
{
    using type = typename T::FluxVariablesCache;
};

template<class T>
using DataCache_t = typename DataCacheType<T>::type;

template<class T,
         bool hasGridDataCache = HasGridDataCache<T>,
         bool hasGridFluxVariablesCache = HasGridFluxVariablesCache<T>>
struct GridCacheType;

template<class T, bool hasGridFluxVariablesCache>
struct GridCacheType<T, true, hasGridFluxVariablesCache>
{
    using type = typename T::GridDataCache;
};

template<class T>
struct GridCacheType<T, false, true>
{
    using type = typename T::GridFluxVariablesCache;
};

template<class T>
using GridCache_t = typename GridCacheType<T>::type;

} // end namespace Dumux::Concept

#endif
