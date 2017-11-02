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
 * \brief Base classes for interaction volume of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_HYBRID_FPS_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_HYBRID_FPS_INTERACTIONVOLUME_HH

#include <dumux/common/math.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/discretization/cellcentered/mpfa/omethodfps/interactionvolume.hh>

namespace Dumux
{
//! Specialization of the interaction volume traits class
template<class TypeTag>
class CCMpfaOHybridFpsInteractionVolumeTraits : public CCMpfaOInteractionVolumeTraits<TypeTag>
{
public:
    // we will use the o method's interaction volume where necessary
    using BoundaryInteractionVolume = CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethodFps>;
};

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa-o method with full pressure support.
 */
template<class TypeTag>
class CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethodFpsHybrid> : public CCMpfaOInteractionVolume<TypeTag, CCMpfaOHybridFpsInteractionVolumeTraits<TypeTag>>
{
    using Traits = CCMpfaOHybridFpsInteractionVolumeTraits<TypeTag>;
    using ParentType = CCMpfaOInteractionVolume<TypeTag, Traits>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using IVSeed = typename Traits::Seed;

public:

    CCMpfaInteractionVolumeImplementation(const IVSeed& seed,
                                          const Problem& problem,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars)
    : ParentType(seed, problem, fvGeometry, elemVolVars)
    {}
};

} // end namespace

#endif
