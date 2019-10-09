// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Definition of the spatial parameters for the linear elasticity problem.
 */

#ifndef DUMUX_ELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_ELASTIC_SPATIAL_PARAMS_HH

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/material/spatialparams/fvelastic.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsTests
 * \brief Definition of the spatial parameters for the linear elasticity problem.
 */
template<class Scalar, class GridGeometry>
class ElasticSpatialParams : public FVSpatialParamsElastic< Scalar,
                                                            GridGeometry,
                                                            ElasticSpatialParams<Scalar, GridGeometry> >
{
    using ThisType = ElasticSpatialParams<Scalar, GridGeometry>;
    using ParentType = FVSpatialParamsElastic<Scalar, GridGeometry, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Export the type of the lame parameters
    using LameParams = Dumux::LameParams<Scalar>;

    ElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        lameParams_.setLambda(3e9);
        lameParams_.setMu(3e9);
    }

    //! Defines the Lame parameters.
    const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const
    { return lameParams_; }

private:
    LameParams lameParams_;
};
} // end namespace Dumux
#endif
