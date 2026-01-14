// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for low-dimensional domains embedded in the bulk domain.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_LINE_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_LINE_HH

#include <vector>

#include <dumux/common/tag.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3dbase.hh>

namespace Dumux {

namespace Embedded1d3dCouplingMode {
struct Line : public Utility::Tag<Line> {
    static std::string name() { return "line"; }
};

inline constexpr Line line{};
} // end namespace Embedded1d3dCouplingMode

// forward declaration
template<class MDTraits, class CouplingMode>
class Embedded1d3dCouplingManager;

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 * \note Specialization for coupling method using line sources with 3d quantities evaluated on the line
 * \note This is the simplest method but it mathematically not well defined as the 3d quantity is evaluated
 *       where the solution to the continuous problem has a singularity
 */
template<class MDTraits>
class Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Line>
: public Embedded1d3dCouplingManagerBase<MDTraits, Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Line>>
{
    using ThisType = Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Line>;
    using ParentType = Embedded1d3dCouplingManagerBase<MDTraits, ThisType>;

public:
    static constexpr Embedded1d3dCouplingMode::Line couplingMode{};

    using ParentType::ParentType;
};

//! we support multithreaded assembly
template<class MDTraits>
struct CouplingManagerSupportsMultithreadedAssembly<Embedded1d3dCouplingManager<MDTraits, Embedded1d3dCouplingMode::Line>>
: public std::true_type {};

} // end namespace Dumux

#endif
