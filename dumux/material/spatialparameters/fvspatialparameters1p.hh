// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
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
#ifndef DUMUX_FV_SPATIAL_PARAMETERS_ONE_P_HH
#define DUMUX_FV_SPATIAL_PARAMETERS_ONE_P_HH

#include <dumux/material/spatialparams/fvspatialparams1p.hh>

#warning include dumux/material/spatialparams/fvspatialparams1p.hh instead

namespace Dumux
{

/**
 * \brief The base class for spatial parameters of problems using the
 *        fv method.
 */
template<class TypeTag>
class FVSpatialParametersOneP : public FVSpatialParamsOneP<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    DUNE_DEPRECATED_MSG("use FVSpatialParamsOneP instead")
    FVSpatialParametersOneP(const GridView &gridView)
    : FVSpatialParamsOneP<TypeTag>(gridView)
    { }
};

} // namespace Dumux

#endif
