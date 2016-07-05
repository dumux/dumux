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
 * \brief Class providing iterators over sub control volumes and sub control
 *        volume faces of an element.
 */
#ifndef DUMUX_DISCRETIZATION_FV_ELEMENTGEOMETRY_HH
#define DUMUX_DISCRETIZATION_FV_ELEMENTGEOMETRY_HH

#include <dune/common/iteratorrange.hh>
#include <dune/common/iteratorfacades.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(SubControlVolume);
NEW_PROP_TAG(SubControlVolumeFace);
NEW_PROP_TAG(FVElementGeometry);
NEW_PROP_TAG(FVElementGeometryVector);
}

/*!
 * \ingroup Discretization
 * \brief An iterator over sub control volumes
 */
template<class SubControlVolume, class Vector, class FVElementGeometryVector>
class ScvIterator : public Dune::ForwardIteratorFacade<ScvIterator<SubControlVolume,
                                                                   Vector,
                                                                   FVElementGeometryVector>,
                                                       const SubControlVolume>
{
    using ThisType = ScvIterator<SubControlVolume, Vector, FVElementGeometryVector>;
    using Iterator = typename Vector::const_iterator;
public:
    ScvIterator(const Iterator& it, const FVElementGeometryVector& fvGeometryVector)
    : it_(it), fvGeometryVector_(&fvGeometryVector) {}

    //! default constructor
    ScvIterator() : it_(Iterator()), fvGeometryVector_(nullptr) {}

    const SubControlVolume& dereference() const
    {
        return fvGeometryVector_->subControlVolume(*it_);
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
    const FVElementGeometryVector* fvGeometryVector_;
};

/*!
 * \ingroup ImplcititModel
 * \brief An iterator over sub control volume faces
 */
template<class SubControlVolumeFace, class Vector, class FVElementGeometryVector>
class ScvfIterator : public Dune::ForwardIteratorFacade<ScvfIterator<SubControlVolumeFace,
                                                                     Vector,
                                                                     FVElementGeometryVector>,
                                                        const SubControlVolumeFace>
{
    using ThisType = ScvfIterator<SubControlVolumeFace, Vector, FVElementGeometryVector>;
    using Iterator = typename Vector::const_iterator;
public:
    ScvfIterator(const Iterator& it, const FVElementGeometryVector& fvGeometryVector)
    : it_(it), fvGeometryVector_(&fvGeometryVector) {}

    //! default constructor
    ScvfIterator() : it_(Iterator()), fvGeometryVector_(nullptr) {}

    const SubControlVolumeFace& dereference() const
    {
        return fvGeometryVector_->subControlVolumeFace(*it_);
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
    const FVElementGeometryVector* fvGeometryVector_;
};

/*!
 * \ingroup ImplcititModel
 * \brief Provide iterators over sub control volumes and sub control
 *        volume faces of an element..
 */
template<class TypeTag>
class FVElementGeometry
{
    using ThisType = Dumux::FVElementGeometry<TypeTag>;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometryVector = typename GET_PROP_TYPE(TypeTag, FVElementGeometryVector);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using ScvIterator = Dumux::ScvIterator<SubControlVolume, std::vector<IndexType>, FVElementGeometryVector>;
    using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, std::vector<IndexType>, FVElementGeometryVector>;

public:
    // Constructor with vectors
    FVElementGeometry(const FVElementGeometryVector& fvGeometryVector,
                      const std::vector<IndexType>& scvIndices,
                      const std::vector<IndexType>& scvfIndices)
    : fvGeometryVector_(fvGeometryVector), scvIndices_(scvIndices), scvfIndices_(scvfIndices)
    {}

    //! iterator range for sub control volumes
    inline Dune::IteratorRange<ScvIterator> scvs() const
    {
        return Dune::IteratorRange<ScvIterator>(ScvIterator(scvIndices_.begin(), fvGeometryVector_),
                                                ScvIterator(scvIndices_.end(), fvGeometryVector_));
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    {
        return scvIndices_.size();
    }

    //! iterator range for sub control volume faces
    inline Dune::IteratorRange<ScvfIterator> scvfs() const
    {
        return Dune::IteratorRange<ScvfIterator>(ScvfIterator(scvfIndices_.begin(), fvGeometryVector_),
                                                 ScvfIterator(scvfIndices_.end(), fvGeometryVector_));
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    {
        return scvfIndices_.size();
    }

private:
    const FVElementGeometryVector& fvGeometryVector_;
    std::vector<IndexType> scvIndices_;
    std::vector<IndexType> scvfIndices_;
};

}

#endif
