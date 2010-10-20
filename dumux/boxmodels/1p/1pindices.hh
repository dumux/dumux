// $Id: 1pproperties.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief  Defines the indices for the one-phase box model.
 */
#ifndef DUMUX_1P_INDICES_HH
#define DUMUX_1P_INDICES_HH

namespace Dumux
{
/*!
 * \addtogroup OnePBoxModel
 */
// \{

/*!
 * \brief Indices for the one-phase model.
 */
struct OnePIndices
{
    static const int pressureIdx = 0;
};

// \}
} // end namepace

#endif
