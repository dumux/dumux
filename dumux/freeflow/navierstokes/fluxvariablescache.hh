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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::FreeFlowFluxVariablesCache
 */
#ifndef DUMUX_FREEFLOW_IMPLICIT_FLUXVARIABLESCACHE_HH
#define DUMUX_FREEFLOW_IMPLICIT_FLUXVARIABLESCACHE_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethod Method>
class FreeFlowFluxVariablesCacheImplementation
{};

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables cache classes for the Navier-Stokes model.
 *        Store flux stencils and data required for flux calculation
 */
template<class TypeTag>
using FreeFlowFluxVariablesCache = FreeFlowFluxVariablesCacheImplementation<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables cache classes for the Navier-Stokes model.
 *        Store flux stencils and data required for flux calculation. <BR>
 *        Specialization for the staggered grid discretization.
 */
template<class TypeTag>
class FreeFlowFluxVariablesCacheImplementation<TypeTag, DiscretizationMethod::staggered>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Do nothing for the staggered grid specialization.
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace &scvf)
    {}
};

} // end namespace

#endif
