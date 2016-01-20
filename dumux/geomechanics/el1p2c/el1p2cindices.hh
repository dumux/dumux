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
 *        the one-phase two-component linear elasticity model.
 */

#ifndef DUMUX_ELASTIC1P2C_INDICES_HH
#define DUMUX_ELASTIC1P2C_INDICES_HH

#include <dumux/geomechanics/elastic/elasticindices.hh>
#include <dumux/porousmediumflow/1p2c/implicit/indices.hh>

namespace Dumux
{
// \{

/*!
 * \ingroup ElOnePTwoCBoxModel
 * \ingroup ImplicitIndices
 * \brief The indices for the one-phase two-component linear elasticity model.
 *
 * This class inherits from the OnePTwoCIndices and from the ElasticIndices
 */
template <class TypeTag>
// PVOffset is set to 0 for the OnePTwoCIndices and to 2 for the ElasticIndices since
// the first two primary variables are the primary variables of the one-phase two-component
// model followed by the primary variables of the elastic model
class ElOnePTwoCIndices : public OnePTwoCIndices<TypeTag, 0>, public ElasticIndices<2>
{
};

} // namespace Dumux

#endif

