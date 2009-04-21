/*****************************************************************************
 *   Copyright (C) 2008,2009 by Melanie Darcis                               *
 *                              Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: melanie.darcis _at_ iws.uni-stuttgart.de                         *
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
 * \brief This template class contains the quantities which are are
 *        constant within a finite element in the non-isothermal
 *        two-phase, two-component model.
 */
 */
#ifndef DUMUX_2P2CNI_ELEMENT_DATA_HH
#define DUMUX_2P2CNI_ELEMENT_DATA_HH

#include <dumux/new_models/2p2c/2p2celementdata.hh>

namespace Dune
{

/*!
 * \brief This template class contains the quantities which are are
 *        constant within a finite element in the non-isothermal
 *        two-phase, two-component model.
 */
template <class TwoPTwoCNITraits, 
          class Problem>
class TwoPTwoCNIElementData : public TwoPTwoCElementData<TwoPTwoCNITraits, Problem>
{
};

} // end namepace

#endif
