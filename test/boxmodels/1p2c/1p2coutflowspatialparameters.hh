// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Definition of the spatial parameters for the 1p2c
 *        outlfow problem.
 *        DEPRECATED use OnePTwoCOutflowSpatialParams
 */
#ifndef DUMUX_1P2C_OUTFLOW_SPATIAL_PARAMETERS_HH
#define DUMUX_1P2C_OUTFLOW_SPATIAL_PARAMETERS_HH

#include "1p2coutflowspatialparams.hh"
#warning include 1p2coutflowspatialparams.hh instead

namespace Dumux
{

template<class TypeTag>
class OnePTwoCOutflowSpatialParameters : public OnePTwoCOutflowSpatialParams<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    DUMUX_DEPRECATED_MSG("use OnePTwoCOutflowSpatialParams instead")
    OnePTwoCOutflowSpatialParameters(const GridView &gridView)
    : OnePTwoCOutflowSpatialParams<TypeTag>(gridView)
};

}
#endif
