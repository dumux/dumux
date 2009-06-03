//$Id:$
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 *        finite element in the richards model.
 *
 * For the plain richards model everything is given on the finite
 * volumes, so this class is empty.
 */
#ifndef DUMUX_RICHARDS_ELEMENT_DATA_HH
#define DUMUX_RICHARDS_ELEMENT_DATA_HH

namespace Dune
{
/*!
 * \brief This template class contains the quantities which are are
 *        constant within a finite element in the richards model.
 *
 * For the plain richards model everything is given on the finite
 * volumes, so this class is empty.
 */
template <class TypeTag>
class RichardsElementData
{
};

} // end namepace

#endif
