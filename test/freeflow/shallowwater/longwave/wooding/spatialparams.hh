// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters of Wooding's V-catchment for the long-wave models
 */
#ifndef DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_SPATIALPARAMS_HH
#define DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_SPATIALPARAMS_HH

#include <memory>

#include <dumux/freeflow/spatialparams.hh>

#include "catchment.hh"

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters of Wooding's V-catchment for the long-wave models
 */
template<class GridGeometry, class Scalar>
class WoodingSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, WoodingSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = WoodingSpatialParams<GridGeometry, Scalar>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    WoodingSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    Scalar bedSurface(const Element& element, const SubControlVolume& scv) const
    { return catchment_.bedElevation(scv.dofPosition()); }

    Scalar manningN(const Element& element) const
    { return catchment_.manningN(element.geometry().center()); }

    const Wooding::Catchment<Scalar>& catchment() const
    { return catchment_; }

private:
    Wooding::Catchment<Scalar> catchment_;
};

} // end namespace Dumux

#endif
