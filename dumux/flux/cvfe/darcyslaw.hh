// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Cvfe
 * \brief Specialization of Darcy's Law for the control-volume finite-element method.
 *
 * This file contains the data which is required to calculate
 * volume and mass fluxes of fluid phases over a face of a finite volume by means
 * of the Darcy approximation.
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CVFE_DARCYS_LAW_HH

#include <dumux/flux/box/darcyslaw.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class DarcysLawImplementation;

/*!
 * \ingroup Cvfe
 * \brief Specialization of Darcy's Law for the control-volume finite-element method
 */
template<class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::Cvfe>
: public BoxDarcysLaw<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::GridGeometry>>
{};

} // end namespace Dumux

#endif
