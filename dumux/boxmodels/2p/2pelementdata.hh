/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *        finite element in the two-phase model.
 *
 * For the plain twophase model everything is given on the finite
 * volumes, so this class is empty.
 */
#ifndef DUMUX_2P_ELEMENT_DATA_HH
#define DUMUX_2P_ELEMENT_DATA_HH

namespace Dune
{
/*!
 * \ingroup TwoPBoxModel
 * \brief This template class contains the quantities which are are
 *        constant within a finite element in the two-phase model.
 *
 * For the plain two-phase model everything is given on the finite
 * volumes, so this class is empty.
 */
template <class TypeTag>
class TwoPElementData
{
};

} // end namepace

#endif
