// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief A thin axis-aligned box geometry
 */
#ifndef DUMUX_GEOMETRY_BOUNDINGBOX_GEOMETRY_HH
#define DUMUX_GEOMETRY_BOUNDINGBOX_GEOMETRY_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dumux/geometry/geometrytraits.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief A thin axis-aligned box geometry
 *
 * Provides the minimal geometry interface required by Dumux::BoundingBoxTree
 * (corner access and geometry type), backed only by the lower and upper
 * corner coordinates. This is useful to build a bounding box tree over a set
 * of plain bounding boxes (e.g. the per-process partition boxes of a
 * distributed tree) without the overhead of a full Dune geometry.
 */
template<class ct, int dimworld>
class BoundingBoxGeometry
{
public:
    static constexpr int mydimension = dimworld;
    static constexpr int coorddimension = dimworld;
    using ctype = ct;
    using GlobalCoordinate = Dune::FieldVector<ct, dimworld>;

    BoundingBoxGeometry() = default;

    //! Construct from the lower and upper corner of the box
    BoundingBoxGeometry(const GlobalCoordinate& lower, const GlobalCoordinate& upper)
    : lower_(lower), upper_(upper) {}

    //! the number of corners of the box
    int corners() const { return 1 << dimworld; }

    //! the i-th corner in Dune reference-cube ordering
    GlobalCoordinate corner(int i) const
    {
        GlobalCoordinate c;
        for (int d = 0; d < dimworld; ++d)
            c[d] = (i & (1 << d)) ? upper_[d] : lower_[d];
        return c;
    }

    //! the geometry type (an axis-aligned cube)
    Dune::GeometryType type() const { return Dune::GeometryTypes::cube(dimworld); }

    //! the lower (min) corner of the box
    const GlobalCoordinate& lower() const { return lower_; }
    //! the upper (max) corner of the box
    const GlobalCoordinate& upper() const { return upper_; }

private:
    GlobalCoordinate lower_, upper_;
};

//! \brief Specialization for the thin axis-aligned box geometry
template<class ct, int dimworld>
struct GeometryTraits<BoundingBoxGeometry<ct, dimworld>>
{ static constexpr bool axisAligned = true; };

} // end namespace Dumux

#endif
