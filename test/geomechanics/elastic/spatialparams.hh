// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Definition of the spatial parameters for the linear elasticity problem.
 */

#ifndef DUMUX_ELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_ELASTIC_SPATIAL_PARAMS_HH

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/geomechanics/elastic/fvspatialparams.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsTests
 * \brief Definition of the spatial parameters for the linear elasticity problem.
 */
template<class Scalar, class GridGeometry>
class ElasticSpatialParams : public FVElasticSpatialParams< GridGeometry,
                                                            Scalar,
                                                            ElasticSpatialParams<Scalar, GridGeometry> >
{
    using ThisType = ElasticSpatialParams<Scalar, GridGeometry>;
    using ParentType = FVElasticSpatialParams<GridGeometry, Scalar, ThisType>;

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
