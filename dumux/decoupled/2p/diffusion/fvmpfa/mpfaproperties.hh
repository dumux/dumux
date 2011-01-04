// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2007-2010 by Yufei Cao                                    *
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
 * \ingroup MPFA2p
 * \file
 *
 * \brief Properties for the MPFA-O method.
 */
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
