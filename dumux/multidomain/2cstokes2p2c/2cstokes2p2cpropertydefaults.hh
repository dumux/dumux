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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup MultidomainModel
 *
 * \brief Defines default values for the properties required by the
 *        coupled 2cstokes2p2c model.
 */
#ifndef DUMUX_TWOCSTOKESTWOPTWOC_PROPERTY_DEFAULTS_HH
#define DUMUX_TWOCSTOKESTWOPTWOC_PROPERTY_DEFAULTS_HH

#include "2cstokes2p2cproperties.hh"
#include <dumux/multidomain/2cstokes2p2c/2cstokes2p2cnewtoncontroller.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////

// Specify the multidomain gridview
SET_TYPE_PROP(TwoCStokesTwoPTwoC, GridView,
              typename GET_PROP_TYPE(TypeTag, MultiDomainGrid)::LeafGridView);

// Specify the type of the used solution vector
SET_TYPE_PROP(TwoCStokesTwoPTwoC, SolutionVector,
              typename GET_PROP_TYPE(TypeTag, MultiDomainGridOperator)::Traits::Domain);

// Specif the used Newton controller
SET_TYPE_PROP(TwoCStokesTwoPTwoC, NewtonController, Dumux::TwoCStokesTwoPTwoCNewtonController<TypeTag>);

// Set this to one here (must fit to the structure of the coupled matrix which has block length 1)
SET_INT_PROP(TwoCStokesTwoPTwoC, NumEq, 1);

// Specify the used boundary layer model
SET_INT_PROP(TwoCStokesTwoPTwoC, BoundaryLayerModel, 0);

// Specify the used mass transfer model
SET_INT_PROP(TwoCStokesTwoPTwoC, MassTransferModel, 0);

} // end namespace properties

} // end namespace Dumux

#endif
