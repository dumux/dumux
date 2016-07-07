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
    // This class in not default-constructible
    FVElementGeometry() = delete;

    // Constructor with vectors
    FVElementGeometry(const FVElementGeometryVector& localFvGeometry,
                      const std::vector<IndexType>& scvIndices,
                      const std::vector<IndexType>& scvfIndices)
    : localFvGeometry_(localFvGeometry),
      scvIndices_(scvIndices),
      scvfIndices_(scvfIndices)
    {}

    //! iterator range for sub control volumes
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (const auto& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<ScvIterator>
    scvs(const FVElementGeometry& g)
    {
        return Dune::IteratorRange<ScvIterator>(ScvIterator(g.scvIndices().begin(), g.localFvGeometry()),
                                                ScvIterator(g.scvIndices().end(), g.localFvGeometry()));
    }

    //! iterator range for sub control volume faces
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (const auto& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<ScvfIterator>
    scvfs(const FVElementGeometry& g)
    {
        return Dune::IteratorRange<ScvfIterator>(ScvfIterator(g.scvfIndices().begin(), g.localFvGeometry()),
                                                 ScvfIterator(g.scvfIndices().end(), g.localFvGeometry()));
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    {
        return scvIndices_.size();
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    {
        return scvfIndices_.size();
    }

private:

    //! The sub control volume indices
    //! Depending on the type of discretization and caching
    //! these may be local or global indices
    const std::vector<IndexType>& scvIndices() const
    { return scvIndices_; }

    //! The sub control volume face indices
    //! Depending on the type of discretization and caching
    //! these may be local or global indices
    const std::vector<IndexType>& scvfIndices() const
    { return scvfIndices_; }

    //! The LocalFvGeometry object this FvElementGeometry belongs to
    const FVElementGeometryVector& localFvGeometry() const
    { return localFvGeometry_; }

    const FVElementGeometryVector& localFvGeometry_;
    std::vector<IndexType> scvIndices_;
    std::vector<IndexType> scvfIndices_;
};

}

#endif
