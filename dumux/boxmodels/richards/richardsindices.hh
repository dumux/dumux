// $Id: richardsproperties.hh 3840 2010-07-15 10:14:15Z bernd $
/*****************************************************************************
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
 * \brief Indices for the Richards model.
 */
#ifndef DUMUX_RICHARDS_INDICES_HH
#define DUMUX_RICHARDS_INDICES_HH

namespace Dumux
{
/*!
 * \addtogroup RichardsModel
 */
// \{

/*!
 * \brief Indices for the Richards model.
 */
struct RichardsIndices
{
    //////////
    // primary variable indices
    //////////
    
    //! wetting phase pressure
    static const int pwIdx = 0;

    //////////
    // equation indices
    //////////
    //! continuity
    static const int contiEqIdx = 0;

    //////////
    // phase indices
    //////////
    static const int wPhaseIdx = 0;
    static const int nPhaseIdx = 1;
};
// \}

} // end namepace

#endif
