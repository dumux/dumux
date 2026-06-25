// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef BENCHMARKS_WELLBORECOUPLINGMANAGER_HH
#define BENCHMARKS_WELLBORECOUPLINGMANAGER_HH

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

namespace Dumux {

/*!
 * \brief Extends Embedded1d3dCouplingManager with an innerRadius() accessor.
 *
 * The spatial params of the low-dim domain must provide innerRadius(eIdx).
 * radius() (inherited) returns the outer radius used for coupling surface placement.
 */
template<class MDTraits, class CouplingMode>
class WellboreCouplingManager
: public Embedded1d3dCouplingManager<MDTraits, CouplingMode>
{
    using ParentType = Embedded1d3dCouplingManager<MDTraits, CouplingMode>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

public:
    using ParentType::ParentType;

    Scalar innerRadius(std::size_t id) const
    {
        const auto& data = this->pointSourceData()[id];
        return this->problem(lowDimIdx).spatialParams().innerRadius(data.lowDimElementIdx());
    }
};

} // end namespace Dumux

#endif
