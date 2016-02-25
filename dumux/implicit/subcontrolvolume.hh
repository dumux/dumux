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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Base class for a sub control volume
 */
#ifndef DUMUX_SUBCONTROLVOLUME_HH
#define DUMUX_SUBCONTROLVOLUME_HH

#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for a sub control volume, i.e a part of the control
 *        volume we are making the balance for.
 */
template<class G, typename I>
class SubControlVolume
{
public:
    // exported types
    using Geometry = G;
    using IndexType = I;

private:
    using Scalar = typename Geometry::ctype;
    enum { dimworld = Geometry::coorddimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    /*\brief The constructor
     * \param geometry The geometry of the sub control volume
     * \param elementIndex The index of the element the scv is part of
     * \param indexInElement The local index of the scv in the element
     */
    SubControlVolume(Geometry geometry,
                     IndexType elementIndex,
                     IndexType indexInElement)
    : geometry_(std::move(geometry)),
      elementIndex_(std::move(elementIndex)),
      indexInElement_(std::move(indexInElement)) {}

    // The move contructor
    SubControlVolume(SubControlVolume&& other)
    : geometry_(std::move(other.geometry_)),
      elementIndex_(std::move(other.elementIndex_)),
      indexInElement_(std::move(other.indexInElement_)) {}

    // The copy contructor
    SubControlVolume(const SubControlVolume& other)
    : geometry_(other.geometry_),
      elementIndex_(other.elementIndex_),
      indexInElement_(other.indexInElement_) {}

    //! The center of the sub control volume
    GlobalPosition center() const
    {
        return geometry_.center();
    }

    //! The volume of the sub control volume
    Scalar volume() const
    {
        return geometry_.volume();
    }

    //! The geometry of the sub control volume
    // e.g. for integration
    const Geometry& geometry() const
    {
        return geometry_;
    }

    IndexType elementIndex() const
    {
        return elementIndex_;
    }

    IndexType indexInElement() const
    {
        return indexInElement_;
    }

private:
    Geometry geometry_;
    IndexType elementIndex_;
    IndexType indexInElement_;
};

} // end namespace

#endif
