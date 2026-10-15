// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef BENCHMARKS_ROCK_SPATIALPARAMS_HH
#define BENCHMARKS_ROCK_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class RockSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, RockSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RockSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    RockSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
    * \brief Defines the porosity \f$\mathrm{[-]}\f$.
    *
    * \param element The current finite element
    * \param scv The sub control volume
    * \param elemSol The current element solution vector
    */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return 0.0;
    }
};

} // end namespace Dumux

#endif
