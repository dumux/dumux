// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

/*!
 * \ingroup Discretization
 * \brief Iterators over sub control volume faces of an fv geometry and a given sub control volume
 * \note usage: for(const auto& scvf : scvfs(fvGeometry, scv))
 */
template<class SubControlVolumeFace, class Vector, class FVElementGeometry>
class SkippingScvfIterator : public Dune::ForwardIteratorFacade<SkippingScvfIterator<SubControlVolumeFace,
                                                                                     Vector,
                                                                                     FVElementGeometry>,
                                                                 const SubControlVolumeFace>
{
    using ThisType = SkippingScvfIterator<SubControlVolumeFace, Vector, FVElementGeometry>;
    using Iterator = typename Vector::const_iterator;
public:

    SkippingScvfIterator() : it_(Iterator()), fvGeometryPtr_(nullptr) {}

    static ThisType makeBegin(const Vector& vector, const FVElementGeometry& fvGeometry, const std::size_t scvIdx)
    {
        auto begin = vector.begin();
        const auto end = vector.end();

        while (begin != end && fvGeometry.scvf(*begin).insideScvIdx() != scvIdx)
            ++begin;

        return SkippingScvfIterator(begin, end, fvGeometry, scvIdx);
    }

    static ThisType makeEnd(const Vector& vector, const FVElementGeometry& fvGeometry, const std::size_t scvIdx)
    {
        return SkippingScvfIterator(vector.end(), vector.end(), fvGeometry, scvIdx);
    }

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
        ++it_;
        while (it_ != itEnd_ && dereference().insideScvIdx() != scvIdx_)
            ++it_;
    }

private:

    SkippingScvfIterator(const Iterator& itBegin, const Iterator& itEnd, const FVElementGeometry& fvGeometry, const std::size_t scvIdx)
    : it_(itBegin), fvGeometryPtr_(&fvGeometry), itEnd_(itEnd), scvIdx_(scvIdx) {}

    Iterator it_;
    const FVElementGeometry* fvGeometryPtr_;
    const Iterator itEnd_;
    std::size_t scvIdx_;
};

} // end namespace Dumux

#endif
