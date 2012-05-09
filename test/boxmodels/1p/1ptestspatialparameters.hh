// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 *        DEPRECATED use OnePTestSpatialParams
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMETERS_HH
#define DUMUX_1P_TEST_SPATIALPARAMETERS_HH

#include "1ptestspatialparams.hh"

#warning include 1ptestspatialparams.hh instead

namespace Dumux
{

template<class TypeTag>
class OnePTestSpatialParameters : public OnePTestSpatialParams<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    DUMUX_DEPRECATED_MSG("use OnePTestSpatialParams instead")
    OnePTestSpatialParameters(const GridView &gridView)
    : OnePTestSpatialParams<TypeTag>(gridView)
    { }
};

} // end namespace
#endif

