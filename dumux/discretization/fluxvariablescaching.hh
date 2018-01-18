// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Discretization
 * \brief Classes related to flux variables caching
 */
#ifndef DUMUX_DISCRETIZATION_FLUXVAR_CACHING_HH
#define DUMUX_DISCRETIZATION_FLUXVAR_CACHING_HH

#include <dumux/common/properties.hh>

namespace Dumux {
namespace FluxVariablesCaching {

#ifndef DOXYGEN // hide the empty caches from doxygen

//! The empty filler class corresponding to EmptyCache
template<class TypeTag>
class EmptyCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
public:
    //! Fill, signature for advection filler
    template<class FluxVariablesCacheFiller>
    static void fill(FluxVariablesCache& scvfFluxVarsCache,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolumeFace& scvf,
                     const FluxVariablesCacheFiller& fluxVarsCacheFiller)
    {}

    //! Fill, signature for diffusion filler
    template<class FluxVariablesCacheFiller>
    static void fill(FluxVariablesCache& scvfFluxVarsCache,
                     unsigned int phaseIdx, unsigned int compIdx,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolumeFace& scvf,
                     const FluxVariablesCacheFiller& fluxVarsCacheFiller)
    {}
};

// an empty cache filler
// \note Never use the _EmptyCache directly as it lead to ambiguous definitions
template<class TypeTag> struct _EmptyCache
{ using Filler = EmptyCacheFiller<TypeTag>; };

#endif // DOXYGEN

/*!
 * \ingroup Discretization
 * \brief Empty caches to use in a constitutive flux law/process, e.g. Darcy's law
 */
template<class TypeTag> class EmptyAdvectionCache : public _EmptyCache<TypeTag> {};
template<class TypeTag> class EmptyDiffusionCache : public _EmptyCache<TypeTag> {};
template<class TypeTag> class EmptyHeatConductionCache : public _EmptyCache<TypeTag> {};

} // end namespace FluxVariablesCaching
} // end namespace Dumux

#endif
