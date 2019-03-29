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
 *
 * \brief Defines the primary variable and equation indices used by
 *        the Richards linear elasticity model.
 */
#ifndef DUMUX_ELASTICRICHARDS_INDICES_HH
#define DUMUX_ELASTICRICHARDS_INDICES_HH

#include <dumux/geomechanics/elastic/indices.hh>
#include <dumux/porousmediumflow/richards/implicit/indices.hh>
#include "properties.hh"

namespace Dumux
{
// \{

namespace Properties
{

/*!
 * \ingroup ElRichardsBoxModel
 * \ingroup ImplicitIndices
 * \brief The indices for the Richards linear elasticity model.
 *
 * This class inherits from the RichardsIndices and from the ElasticIndices
 */

// PVOffset is set to 0 for the RichardsIndices and to 2 for the ElasticIndices since
// the first two primary variables are the primary variables of the Richards
// model followed by the primary variables of the elastic model
// template <class TypeTag,
// int formulation = 0,
// int PVOffset =  2>
// //class ElRichardsIndices : public ElasticIndices<PVOffset>, public RichardsIndices<TypeTag,0>
// class ElRichardsIndices : public ElasticIndices<PVOffset>, public RichardsIndices<TypeTag>
// //{};


//template <class TypeTag>
//struct ElRichardsIndices
template <class TypeTag,
int formulation = 0,
int PVOffset =  1>
class ElRichardsIndices : public ElasticIndices<PVOffset>, public RichardsIndices<TypeTag>
{};

}

// //here on from richards
// {
//     typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
//
//     //////////
//     // primary variable indices
//     //////////
//
//     //! Primary variable index for the wetting phase pressure
//     static const int pwIdx = 0;
//     //! Primary variable index for the wetting phase pressure head (used for pressure head formulation)
//     static const int hIdx = 0;
//     //////////
//     // equation indices
//     //////////
//     //! Equation index for the mass conservation of the wetting phase
//     static const int contiEqIdx = 0;
//
//     //////////
//     // phase indices
//     //////////
//     static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the wetting phase;
//     static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the non-wetting phase;
// };

}


#endif
