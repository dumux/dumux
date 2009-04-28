/*****************************************************************************
 *   Copyright (C) 2008 by Onur Dogan, Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_RICHARDSTRAITS_HH
#define DUMUX_RICHARDSTRAITS_HH

#include <dune/common/fvector.hh>

namespace Dune
{
///////////////////////////////////////////////////////////////////////////
// The traits for the twophase model
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief The traits for the richards model.
 */

class RichardsTraits
{
public:
    enum {
        numEq = 1,        //!< Number of primary variables
        numPhases = 1,    //!< Number of fluid phases
    };
    enum {
        pWIdx = 0  //!< Index for the fluid pressure in a field vector
    };
};

}

#endif

