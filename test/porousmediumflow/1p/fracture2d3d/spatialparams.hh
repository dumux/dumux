// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model
 */

#ifndef DUMUX_ONEP_FRACTURE_TEST_SPATIALPARAMS_HH
#define DUMUX_ONEP_FRACTURE_TEST_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model.
 */
template<class GridGeometry, class Scalar>
class FractureSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         FractureSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;

    using ThisType = FractureSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    FractureSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$.
     *
     * \param globalPos The global position
     */
    Scalar permeabilityAtPos (const GlobalPosition& globalPos) const
    { return 1e-12; }

    /*!
     * \brief Returns the porosity \f$[-]\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }
};

} // end namespace Dumux

#endif
