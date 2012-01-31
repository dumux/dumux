// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Holger Class                                 *
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \file
 *
 * \brief Defines the indices used by the 3p3cni box model
 */
#ifndef DUMUX_3P3CNI_INDICES_HH
#define DUMUX_3P3CNI_INDICES_HH

#include <dumux/boxmodels/3p3c/3p3cindices.hh>

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCNIModel
 */
// \{

/*!
 * \brief Enumerations for the non-isothermal 3-phase 3-component model
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class ThreePThreeCNIIndices : public ThreePThreeCIndices<TypeTag, PVOffset>
{
public:
    static const int temperatureIdx = PVOffset + 3; //! The index for temperature in primary variable vectors.
    static const int energyEqIdx = PVOffset + 3; //! The index for energy in equation vectors.
};

}
#endif
