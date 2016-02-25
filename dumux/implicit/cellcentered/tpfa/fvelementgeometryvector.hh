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

namespace Dumux
{

//! An index to element map
template <class GridView>
class ElementMap
  : public std::vector<typename GridView::Traits::Grid::template Codim<0>::EntitySeed>
{
    using Grid = typename GridView::Traits::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
public:
    ElementMap(const GridView& gridView_) : grid_(gridView_.grid()) {}

    Element element(IndexType eIdx)
    { return grid_.entity((*this)[eIdx]); }

private:
    const Grid& grid_;
};

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag>
class CCTpfaFVElementGeometryVector
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
    FVElementGeometry fvGeometry(IndexType eIdx) const
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
        return *scvfs_[scvfIdx];
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return scvs_.size();
    }

    //! The total number of sun control volume faces
    std::size_t numScvf() const
    {
        return scvfs_.size();
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
                    scvfs_.push_back(std::make_shared<SubControlVolumeFace>(intersection.geometry(),
                                                                            intersection.centerUnitOuterNormal(),
                                                                            scvfIdx,
                                                                            std::vector<IndexType>({eIdx, nIdx}),
                                                                            std::vector<IndexType>({eIdx, nIdx}),
                                                                            false));
                    scvfsIndexSet.push_back(scvfIdx++);
                }
                else
                {
                    scvfs_.push_back(std::make_shared<SubControlVolumeFace>(intersection.geometry(),
                                                                            intersection.centerUnitOuterNormal(),
                                                                            scvfIdx,
                                                                            std::vector<IndexType>({eIdx}),
                                                                            std::vector<IndexType>({eIdx}),
                                                                            true));
                    scvfsIndexSet.push_back(scvfIdx++);
                }
            }

            // Compute the finite volume element geometries
            fvGeometries_.push_back(FVElementGeometry(*this, {eIdx}, scvfsIndexSet));
        }
    }

private:
    GridView gridView_;
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::shared_ptr<SubControlVolume>> scvs_;
    std::vector<std::shared_ptr<SubControlVolumeFace>> scvfs_;
    std::vector<FVElementGeometry> fvGeometries_;
};

} // end namespace

#endif
