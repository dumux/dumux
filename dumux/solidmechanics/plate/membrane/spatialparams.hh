// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MembranePlate
 * \brief Default spatial parameters for the membrane plate model
 */
#ifndef DUMUX_MEMBRANE_PLATE_SPATIAL_PARAMS_HH
#define DUMUX_MEMBRANE_PLATE_SPATIAL_PARAMS_HH

#include <dumux/common/fvspatialparams.hh>
#include <dumux/common/functionfromstringexpression.hh>

namespace Dumux {

/*!
 * \ingroup MembranePlate
 * \brief Spatial parameters for the membrane plate model
 *
 * Reads the membrane tension \f$ T \f$ from the parameter `Problem.Tension`.
 */
template<class GridGeometry, class Scalar>
class MembranePlateSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, MembranePlateSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar, MembranePlateSpatialParams<GridGeometry, Scalar>>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PositionFunction = FunctionFromStringExpression<GlobalPosition::dimension, Scalar>;
public:
    MembranePlateSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , tension_(getParam<std::string>("Problem.Tension", "1.0"), "xy")
    {}

    //! Returns the membrane tension \f$ T(x) \f$ (force per unit length)
    Scalar tension(const GlobalPosition& globalPos) const
    { return tension_(globalPos); }

private:
    PositionFunction tension_;
};

template<class GridGeometry, class Scalar>
class MembranePlateDefaultSpatialParams
: public MembranePlateSpatialParams<GridGeometry, Scalar>
{
    using ParentType = MembranePlateSpatialParams<GridGeometry, Scalar>;
public:
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
