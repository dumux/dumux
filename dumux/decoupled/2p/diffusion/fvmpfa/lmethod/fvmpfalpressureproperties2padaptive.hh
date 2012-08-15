// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
 *   Institute of Applied Analysis and Numerical Simulation                  *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@mathematik.uni-stuttgart.de                   *
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
 * \ingroup MPFA
 * \ingroup Properties
 * \file
 *
 * \brief Properties for the adaptive MPFA-L method.
 */
#ifndef DUMUX_FVMPFALPROPERTIES2P_ADAPTIVE_HH
#define DUMUX_FVMPFALPROPERTIES2P_ADAPTIVE_HH

// dumux environment
#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>
#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>

namespace Dumux
{
namespace Properties
{
NEW_TYPE_TAG(FVMPFALPressureTwoPAdaptive, INHERITS_FROM(PressureTwoP, MPFAProperties));
}
}

#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfal2pfaboundvelocity2padaptive.hh>

namespace Dumux
{
namespace Properties
{
SET_TYPE_PROP(FVMPFALPressureTwoPAdaptive, PressureModel, Dumux::FVMPFAL2PFABoundVelocity2PAdaptive<TypeTag>);
}
}// end of Dune namespace
#endif
