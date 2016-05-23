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
 * \brief Base class for the finite volume geometry vector for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
#ifndef DUMUX_IMPLICIT_TPFA_FV_GEOMETRY_VECTOR_HH
#define DUMUX_IMPLICIT_TPFA_FV_GEOMETRY_VECTOR_HH

#include <dumux/implicit/subcontrolvolume.hh>
#include <dumux/implicit/subcontrolvolumeface.hh>
#include <dumux/implicit/fvelementgeometry.hh>
#include <dumux/common/elementmap.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool enableFVElementGeometryCache>
class CCTpfaFVElementGeometryVector
{};

// specialization in case the FVElementGeometries are stored
template<class TypeTag>
class CCTpfaFVElementGeometryVector<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Constructor
    CCTpfaFVElementGeometryVector(const GridView gridView)
    : gridView_(gridView), elementMap_(gridView) {}

    /* \brief Get the finite volume geometry of an element
     * \note The finite volume geometry offers iterators over the sub control volumes
     *       and the sub control volume faces of an element.
     */
    const FVElementGeometry& fvGeometry(IndexType eIdx) const
    {
        return fvGeometries_[eIdx];
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& subControlVolume(IndexType scvIdx) const
    {
        return *scvs_[scvIdx];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& subControlVolumeFace(IndexType scvfIdx) const
    {
        return scvfs_[scvfIdx];
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return scvs_.size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return scvfs_.size();
    }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    {
        return numBoundaryScvf_;
    }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        scvs_.clear();
        scvfs_.clear();
        fvGeometries_.clear();
        elementMap_.clear();

        // Build the SCV and SCV faces
        IndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        elementMap_.resize(gridView_.size(0));
        scvs_.resize(gridView_.size(0));
        for (const auto& element : elements(gridView_))
        {
            auto eIdx = problem.elementMapper().index(element);
            scvs_[eIdx] = std::make_shared<SubControlVolume>(element.geometry(), eIdx);

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element-wise index sets for finite volume geometry
            std::vector<IndexType> scvfsIndexSet;
            for (const auto& intersection : intersections(gridView_, element))
            {
                if (!intersection.boundary())
                {
                    auto nIdx = problem.elementMapper().index(intersection.outside());
                    scvfs_.push_back(SubControlVolumeFace(intersection.geometry(),
                                                          intersection.geometry().center(),
                                                          intersection.centerUnitOuterNormal(),
                                                          scvfIdx,
                                                          std::vector<IndexType>({eIdx, nIdx}),
                                                          false));
                    scvfsIndexSet.push_back(scvfIdx++);
                }
                else
                {
                    scvfs_.push_back(SubControlVolumeFace(intersection.geometry(),
                                                          intersection.geometry().center(),
                                                          intersection.centerUnitOuterNormal(),
                                                          scvfIdx,
                                                          std::vector<IndexType>({eIdx, gridView_.size(0) + numBoundaryScvf_++}),
                                                          true));
                    scvfsIndexSet.push_back(scvfIdx++);
                }
            }

            // Compute the finite volume element geometries
            fvGeometries_.push_back(FVElementGeometry(*this, {eIdx}, scvfsIndexSet));
        }
    }

    // Binding of an element, called by the local jacobian to prepare element assembly
    // For compatibility reasons with the FVGeometry cache disabled
    void bind(const Element& element) {}
    void bindElement(const Element& element) {}

private:
    GridView gridView_;
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::shared_ptr<SubControlVolume>> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<FVElementGeometry> fvGeometries_;
    IndexType numBoundaryScvf_;
};

// specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class CCTpfaFVElementGeometryVector<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! Constructor
    CCTpfaFVElementGeometryVector(const GridView gridView)
    : gridView_(gridView), elementMap_(gridView), eIdxBound_(-1) {}

    /* \brief Get the finite volume geometry of an element
     * \note The finite volume geometry offers iterators over the sub control volumes
     *       and the sub control volume faces of an element.
     */
    const FVElementGeometry& fvGeometry(IndexType eIdx) const
    {
        return fvGeometries_[getStencilScvIdx_(eIdx)];
    }

    const FVElementGeometry& fvGeometry(const Element& element) const
    {
        return fvGeometries_[getStencilScvIdx_(problem_().elementMapper().index(element))];
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& subControlVolume(IndexType scvIdx) const
    {
        return stencilScvs_[getStencilScvIdx_(scvIdx)];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& subControlVolumeFace(IndexType scvfIdx) const
    {
        return stencilScvfs_[getStencilScvFaceIdx_(scvfIdx)];
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return numScvs_;
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return numScvf_;
    }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    {
        return numBoundaryScvf_;
    }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        problemPtr_ = &problem;
        elementMap_.clear();
        scvFaceToElementMap_.clear();
        eIdxBound_ = -1;

        // Build the SCV and SCV faces
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        numScvs_ = gridView_.size(0);
        elementMap_.resize(numScvs_);
        scvFaceIndices_.resize(numScvs_);
        neighborVolVarIndices_.resize(numScvs_);
        for (const auto& element : elements(gridView_))
        {
            auto eIdx = problem.elementMapper().index(element);

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element-wise index sets for finite volume geometry
            std::vector<IndexType> scvfsIndexSet;
            std::vector<IndexType> neighborVolVarIndexSet;
            for (const auto& intersection : intersections(gridView_, element))
            {
                scvFaceToElementMap_.push_back(eIdx);
                scvfsIndexSet.push_back(numScvf_++);
                if (!intersection.boundary())
                {
                    auto nIdx = problem.elementMapper().index(intersection.outside());
                    neighborVolVarIndexSet.push_back(nIdx);
                }
                else
                    neighborVolVarIndexSet.push_back(numScvs_ + numBoundaryScvf_++);
            }

            // store the sets of indices in the data container
            scvFaceIndices_[eIdx] = scvfsIndexSet;
            neighborVolVarIndices_[eIdx] = neighborVolVarIndexSet;
        }
    }

    // Binding of an element preparing the geometries of the whole stencil
    // called by the local jacobian to prepare element assembly
    void bind(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);
        if (eIdxBound_ == eIdx && fvGeometries_.size() > 1)
            return;

        release_();
        eIdxBound_ = eIdx;

        makeElementGeometries_(element);
        for (const auto& intersection : intersections(gridView_, element))
        {
            if (!intersection.boundary())
                makeElementGeometries_(intersection.outside());
        }
    }

    // Binding of an element preparing the geometries only inside the element
    void bindElement(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);
        if (eIdx == eIdxBound_ || std::find(stencilScvIndices_.begin(), stencilScvIndices_.end(), eIdx) != stencilScvIndices_.end() )
            return;

        release_();
        eIdxBound_ = eIdx;

        makeElementGeometries_(element);
    }


private:

    void release_()
    {
        stencilScvs_.clear();
        stencilScvIndices_.clear();
        stencilScvfs_.clear();
        stencilScvFaceIndices_.clear();
        fvGeometries_.clear();
    }

    void makeElementGeometries_(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);

        std::vector<SubControlVolume> scv({SubControlVolume(element.geometry(), eIdx)});
        stencilScvs_.push_back(scv[0]);
        stencilScvIndices_.push_back(eIdx);

        const auto& scvFaceIndices = scvFaceIndices_[eIdx];
        const auto& neighborVolVarIndices = neighborVolVarIndices_[eIdx];

        int scvfCounter = 0;
        for (const auto& intersection : intersections(gridView_, element))
        {
            stencilScvfs_.push_back(SubControlVolumeFace(intersection.geometry(),
                                                         intersection.geometry().center(),
                                                         intersection.centerUnitOuterNormal(),
                                                         scvFaceIndices[scvfCounter],
                                                         std::vector<IndexType>({eIdx, neighborVolVarIndices[scvfCounter]}),
                                                         intersection.boundary()));

            stencilScvFaceIndices_.push_back(scvFaceIndices[scvfCounter]);
            scvfCounter++;
        }

        // Compute the finite volume element geometries
        fvGeometries_.push_back(FVElementGeometry(*this, {eIdx}, scvFaceIndices));
    }

    const int getStencilScvIdx_(const int scvIdx) const
    {
        auto it = std::find(stencilScvIndices_.begin(), stencilScvIndices_.end(), scvIdx);

        if (it != stencilScvIndices_.end())
            return std::distance(stencilScvIndices_.begin(), it);
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Could not find the sub control volume for scvIdx = " << scvIdx <<
                       ", make sure to properly bind the FVGeometries to the element before using them.");
    }

    const int getStencilScvFaceIdx_(const int scvfIdx) const
    {
        auto it = std::find(stencilScvFaceIndices_.begin(), stencilScvFaceIndices_.end(), scvfIdx);

        if (it != stencilScvFaceIndices_.end())
            return std::distance(stencilScvFaceIndices_.begin(), it);
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Could not find the sub control volume face for scvfIdx = " << scvfIdx <<
                       ", make sure to properly bind the FVGeometries to the to the element before using them.");
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
    GridView gridView_;

    // Information on the global number of geometries
    IndexType numScvs_;
    IndexType numScvf_;
    IndexType numBoundaryScvf_;

    // vectors that stores global data
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::vector<IndexType>> scvFaceIndices_;
    std::vector<std::vector<IndexType>> neighborVolVarIndices_;
    std::vector<IndexType> scvFaceToElementMap_;

    // vectors to store the geometries temporarily after binding an element
    IndexType eIdxBound_;
    std::vector<SubControlVolume> stencilScvs_;
    std::vector<IndexType> stencilScvIndices_;
    std::vector<SubControlVolumeFace> stencilScvfs_;
    std::vector<IndexType> stencilScvFaceIndices_;
    std::vector<FVElementGeometry> fvGeometries_;
};

} // end namespace

#endif
