// $Id: mpfaproperties.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2007-2010 by Yufei Cao                                    *
 *   Institute of Applied Analysis and Numerical Simulation                  *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@mathematik.uni-stuttgart.de                   *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_MPFAOPROPERTIES_HH
#define DUNE_MPFAOPROPERTIES_HH

// dumux environment
#include <dumux/decoupled/2p/2pproperties.hh>

namespace Dumux
{
/*!
 * \brief Indices denoting the different grid types.
 */
struct GridTypes
{
public:
    //SGrid
    static const int sGrid = 0;
    //YaspGrid
    static const int yaspGrid = 1;
    //UGGrid
    static const int ugGrid = 2;
    //ALUGrid
    static const int aluGrid = 3;
};

namespace Properties
{
NEW_TYPE_TAG(MPFAProperties);

NEW_PROP_TAG( GridTypeIndices );
NEW_PROP_TAG( GridImplementation ); //returns kind of grid implementation

SET_INT_PROP(MPFAProperties, GridImplementation, GridTypes::sGrid);

SET_TYPE_PROP(MPFAProperties, GridTypeIndices, GridTypes);

}
}// end of Dune namespace
#endif
