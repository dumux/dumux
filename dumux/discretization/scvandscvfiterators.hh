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
 * \ingroup Discretization
 * \brief Class providing iterators over sub control volumes and sub control volume faces of an element.
 */
#ifndef DUMUX_SCV_AND_SCVF_ITERATORS_HH
#define DUMUX_SCV_AND_SCVF_ITERATORS_HH

#include <dune/common/iteratorrange.hh>
#include <dune/common/iteratorfacades.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Iterators over sub control volumes
 * \note usage: for(const auto& scv : scvs(fvGeometry))
 */
template<class SubControlVolume, class Vector, class FVElementGeometry>
class ScvIterator : public Dune::ForwardIteratorFacade<ScvIterator<SubControlVolume,
                                                                   Vector,
                                                                   FVElementGeometry>,
                                                       const SubControlVolume>
{
    using ThisType = ScvIterator<SubControlVolume, Vector, FVElementGeometry>;
    using Iterator = typename Vector::const_iterator;
public:
    ScvIterator(const Iterator& it, const FVElementGeometry& fvGeometry)
    : it_(it), fvGeometryPtr_(&fvGeometry) {}

    ScvIterator() : it_(Iterator()), fvGeometryPtr_(nullptr) {}

    //! dereferencing yields a subcontrol volume
    const SubControlVolume& dereference() const
    {
        return fvGeometryPtr_->scv(*it_);
    }

    bool equals(const ThisType& other) const
    {
        return it_ == other.it_;
    }

    void increment()
    {
        ++it_;
    }

private:
    Iterator it_;
    const FVElementGeometry* fvGeometryPtr_;
};

/*!
 * \ingroup Discretization
 * \brief Iterators over sub control volume faces of an fv geometry
 * \note usage: for(const auto& scvf : scvfs(fvGeometry))
 */
template<class SubControlVolumeFace, class Vector, class FVElementGeometry>
class ScvfIterator : public Dune::ForwardIteratorFacade<ScvfIterator<SubControlVolumeFace,
                                                                     Vector,
                                                                     FVElementGeometry>,
                                                        const SubControlVolumeFace>
{
    using ThisType = ScvfIterator<SubControlVolumeFace, Vector, FVElementGeometry>;
    using Iterator = typename Vector::const_iterator;
public:
    ScvfIterator(const Iterator& it, const FVElementGeometry& fvGeometry)
    : it_(it), fvGeometryPtr_(&fvGeometry) {}

    ScvfIterator() : it_(Iterator()), fvGeometryPtr_(nullptr) {}

    //! dereferencing yields a subcontrol volume face
    const SubControlVolumeFace& dereference() const
    {
        return fvGeometryPtr_->scvf(*it_);
    }

    bool equals(const ThisType& other) const
    {
        return it_ == other.it_;
    }

    void increment()
    {
        it_++;
    }

private:
    Iterator it_;
    const FVElementGeometry* fvGeometryPtr_;
};

} // end namespace Dumux

#endif
