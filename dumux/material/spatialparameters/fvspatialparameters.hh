// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
#ifndef DUMUX_FV_SPATIAL_PARAMETERS_HH
#define DUMUX_FV_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fvspatialparams.hh>

#warning include dumux/material/spatialparams/fvspatialparams.hh instead

namespace Dumux
{

/**
 * \brief The base class for spatial parameters of a multi-phase problem using the
 *        fv method.
 */
template<class TypeTag>
class FVSpatialParameters: public FVSpatialParams<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    DUNE_DEPRECATED_MSG("use FVSpatialParams instead")
    FVSpatialParameters(const GridView &gridView)
    :FVSpatialParams<TypeTag>(gridView)
    { }
};

} // namespace Dumux

#endif
