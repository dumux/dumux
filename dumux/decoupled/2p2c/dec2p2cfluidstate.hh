// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Benjamin Faigle                                   *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Calculates the 2p2c phase state for compositional models.
 */
#ifndef DUMUX_DEC2P2C_FLUID_STATE_HH
#define DUMUX_DEC2P2C_FLUID_STATE_HH

#include "2p2cproperties.hh"
#include "2p2cfluidstate.hh"

namespace Dumux
{
/*!
 * \ingroup multiphysics multiphase
 * \brief Calculates the phase state from the primary variables in the
 *        sequential 2p2c model.
 *
 *        This boils down to so-called "flash calculation", in this case isothermal and isobaric.
 *
 *  \tparam TypeTag The property Type Tag
 */
template <class TypeTag>
class DecoupledTwoPTwoCFluidState : public TwoPTwoCFluidState<TypeTag>
{
public:
    DUNE_DEPRECATED_MSG("use TwoPTwoCFluidState instead")
    DecoupledTwoPTwoCFluidState() : TwoPTwoCFluidState<TypeTag>() {}
};

} // end namepace

#endif
