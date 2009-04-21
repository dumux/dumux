/*****************************************************************************
 *   Copyright (C) 2008,2009 by Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: klaus.mosthaf _at_ iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file 
 *
 * \brief Contains the quantities which are are constant within a
 *        finite element in the two-phase, two-component model.
 *
 * By for the plain 2p2c model everything is given on the finite
 * volumes, so this class is empty.
 */
#ifndef DUMUX_2P2C_ELEMENT_DATA_HH
#define DUMUX_2P2C_ELEMENT_DATA_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>
#include <dumux/new_models/2p2c/2p2ctraits.hh>
#include <dumux/auxiliary/math.hh>

#include <dumux/auxiliary/apis.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{

/*!
 * \brief This template class contains the quantities which are are
 *        constant within a finite element in the two-phase,
 *        two-component model.
 *
 * By for the plain 2p2c model everything is given on the finite
 * volumes, so this class is empty.
 */
template <class TwoPTwoCTraits, 
          class ProblemT>
class TwoPTwoCElementData
{
};

} // end namepace

#endif
