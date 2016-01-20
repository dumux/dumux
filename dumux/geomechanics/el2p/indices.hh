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
 *        the two-phase linear elasticity model.
 */
#ifndef DUMUX_ELASTIC2P_INDICES_HH
#define DUMUX_ELASTIC2P_INDICES_HH

#include <dumux/geomechanics/elastic/indices.hh>
#include <dumux/porousmediumflow/2p/implicit/indices.hh>

namespace Dumux
{
// \{

namespace Properties
{

/*!
 * \ingroup ElTwoPBoxModel
 * \ingroup ImplicitIndices
 * \brief The indices for the two-phase linear elasticity model.
 *
 * This class inherits from the TwoPIndices and from the ElasticIndices
 */

// PVOffset is set to 0 for the TwoPIndices and to 2 for the ElasticIndices since
// the first two primary variables are the primary variables of the two-phase
// model followed by the primary variables of the elastic model
template <class TypeTag,
int formulation = 0,
int PVOffset =  2>
class ElTwoPIndices : public ElasticIndices<PVOffset>, public TwoPIndices<TypeTag,0>
{};

}
}

#endif
