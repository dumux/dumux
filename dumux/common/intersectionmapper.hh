// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief defines intersection mappers.
 */
#ifndef DUMUX_INTERSECTIONITERATOR_HH
#define DUMUX_INTERSECTIONITERATOR_HH

#include <vector>
#include <unordered_map>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief defines a standard intersection mapper for mapping of global DOFs assigned
 *        to faces. It only works for conforming grids, without hanging nodes.
 */
template<class GridView>
class ConformingGridIntersectionMapper
{
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename GridView::IndexSet::IndexType;

    static constexpr int codimIntersection =  1;
public:

    ConformingGridIntersectionMapper(const GridView& gridView)
    : gridView_(gridView)
    , indexSet_(&gridView.indexSet())
    {}

    void update (const GridView& gridView)
    {
        gridView_ = gridView;
        indexSet_ = &gridView_.indexSet();
    }

    void update (GridView&& gridView)
    {
        gridView_ = std::move(gridView);
        indexSet_ = &gridView_.indexSet();
    }

    //! The total number of intersections
    std::size_t numIntersections() const
    {
        return gridView_.size(1);
    }

    /*!
     * \brief The number of faces for a given element
     *
     * \param element The element
     */
    std::size_t numFaces(const Element& element) const
    {
        return element.subEntities(1);
    }

    GridIndexType globalIntersectionIndex(const Element& element, const std::size_t localFaceIdx) const
    {
        return indexSet_->subIndex(element, localFaceIdx, codimIntersection);
    }

private:
    GridView gridView_;
    const typename GridView::IndexSet* indexSet_;
};

/*!
 * \ingroup Core
 * \brief defines an intersection mapper for mapping of global DOFs assigned
 *        to faces which also works for non-conforming grids and corner-point grids.
 */
template<class GridView>
class NonConformingGridIntersectionMapper
{
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GridIndexType = typename GridView::IndexSet::IndexType;

public:
    NonConformingGridIntersectionMapper(const GridView& gridview)
    : gridView_(gridview),
      indexSet_(&gridView_.indexSet()),
      numIntersections_(gridView_.size(1)),
      intersectionMapGlobal_(gridView_.size(0))
    {

    }

    //! The total number of intersections
    std::size_t numIntersections() const
    {
        return numIntersections_;
    }

    GridIndexType globalIntersectionIndex(const Element& element, const std::size_t localFaceIdx) const
    {
        return (intersectionMapGlobal_[index(element)].find(localFaceIdx))->second; //use find() for const function!
    }

    std::size_t numFaces(const Element& element)
    {
        return intersectionMapGlobal_[index(element)].size();
    }

    void update (const GridView& gridView)
    {
        gridView_ = gridView;
        indexSet_ = &gridView_.indexSet();
        update_();
    }

    void update (GridView&& gridView)
    {
        gridView_ = std::move(gridView);
        indexSet_ = &gridView_.indexSet();
        update_();
    }

private:
    GridIndexType index(const Element& element) const
    {
        return indexSet_->index(element);
    }

    void update_()
    {
        intersectionMapGlobal_.clear();
        intersectionMapGlobal_.resize(gridView_.size(0));

        int globalIntersectionIdx = 0;
        for (const auto& element : elements(gridView_))
        {
            int eIdx = index(element);
            int fIdx = 0;

            // run through all intersections with neighbors
            for (const auto& intersection : intersections(gridView_, element))
            {
                if (intersection.neighbor())
                {
                    auto neighbor = intersection.outside();
                    int eIdxN = index(neighbor);

                    if (element.level() > neighbor.level() || (element.level() == neighbor.level() && eIdx < eIdxN))
                    {
                        int fIdxN = 0;
                        for (const auto& intersectionNeighbor : intersections(gridView_, neighbor))
                        {
                            if (intersectionNeighbor.neighbor())
                            {
                                if (intersectionNeighbor.outside() == element)
                                {
                                    break;
                                }
                            }
                            fIdxN++;
                        }
                        intersectionMapGlobal_[eIdx][fIdx] = globalIntersectionIdx;
                        intersectionMapGlobal_[eIdxN][fIdxN] = globalIntersectionIdx;
                        globalIntersectionIdx++;
                    }
                }
                else
                {
                    intersectionMapGlobal_[eIdx][fIdx] = globalIntersectionIdx;
                    globalIntersectionIdx++;
                }

                fIdx++;
            }
        }
        numIntersections_ = globalIntersectionIdx;
    }

    GridView gridView_;
    const typename GridView::IndexSet* indexSet_;
    unsigned int numIntersections_;
    std::vector<std::unordered_map<int, int> > intersectionMapGlobal_;
};

/*!
 * \ingroup Core
 * \brief defines an intersection mapper for mapping of global DOFs assigned
 *        to faces which also works for adaptive grids.
 */
template<class GridView>
class IntersectionMapper
{
    using Grid = typename GridView::Grid;
    enum {dim=Grid::dimension};
    using Element = typename Grid::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

public:
    IntersectionMapper(const GridView& gridview)
    : gridView_(gridview)
    , elementMapper_(gridView_, Dune::mcmgElementLayout())
    , size_(gridView_.size(1))
    , intersectionMapGlobal_(gridView_.size(0))
    , intersectionMapLocal_(gridView_.size(0))
    {
        const auto element = *gridView_.template begin<0>();

        int fIdx = 0;
        for (const auto& intersection : intersections(gridView_, element))
        {
            int idxInInside = intersection.indexInInside();

            standardLocalIdxMap_[idxInInside] = fIdx;

            fIdx++;
        }
    }

    const ElementMapper& elementMapper() const
    {
        return elementMapper_;
    }

    /*!
     * \brief Map element to array index.
     *
     * \param element Reference to element which should be mapped.
     * \return An index in the range 0 ... Max number of elements in set - 1.
     */
    int index(const Element& element) const
    {
        return elementMapper_.index(element);
    }

    /*!
     * \brief Map interface fIdx'th interface of element index to array index.
     *
     * \param elemIdx Index of element.
     * \param fIdx Index of interface to map.
     * \return An index in the range 0 ... Max number of intersections in set - 1.
     */
    int subIndex(int elemIdx, int fIdx)
    {
        return intersectionMapGlobal_[elemIdx][fIdx];
    }

    /*!
     * \brief Map interface fIdx'th interface of element index to array index.
     *
     * \param elemIdx Index of element.
     * \param fIdx Index of interface to map.
     * \return An index in the range 0 ... Max number of intersections in set - 1.
     */
    int subIndex(int elemIdx, int fIdx) const
    {
        return (intersectionMapGlobal_[elemIdx].find(fIdx))->second;//use find() for const function!
    }

    /*!
     * \brief Map interface fIdx'th interface of element to array index.
     *
     * \param element Reference to element.
     * \param fIdx Index of interface to map.
     * \return An index in the range 0 ... Max number of intersections in set - 1.
     */
    int subIndex(const Element& element, int fIdx)
    {
        return intersectionMapGlobal_[index(element)][fIdx];
    }

    int subIndex(const Element& element, int fIdx) const
    {
        return intersectionMapGlobal_[index(element)].find(fIdx)->second;//use find() for const function!
    }

    int maplocal(int elemIdx, int fIdx)
    {
        return intersectionMapLocal_[elemIdx][fIdx];
    }

    int maplocal(int elemIdx, int fIdx) const
    {
        return (intersectionMapLocal_[elemIdx].find(fIdx))->second;//use find() for const function!
    }


    int maplocal(const Element& element, int fIdx)
    {
        return intersectionMapLocal_[index(element)][fIdx];
    }

    int maplocal(const Element& element, int fIdx) const
    {
        return (intersectionMapLocal_[index(element)].find(fIdx))->second;//use find() for const function!
    }

    // return number intersections
    unsigned int size () const
    {
        return size_;
    }

    unsigned int size (int elemIdx) const
    {
        return intersectionMapLocal_[elemIdx].size();
    }

    unsigned int size (const Element& element) const
    {
        return intersectionMapLocal_[index(element)].size();
    }

    void update (const GridView &gridView)
    {
        gridView_ = gridView;
        update_();
    }

    void update (GridView&& gridView)
    {
        gridView_ = std::move(gridView);
        update_();
    }

protected:
    void update_()
    {
        elementMapper_.update(gridView_);

        intersectionMapGlobal_.clear();
        intersectionMapGlobal_.resize(elementMapper_.size());
        intersectionMapLocal_.clear();
        intersectionMapLocal_.resize(elementMapper_.size());

        for (const auto& element : elements(gridView_))
        {
            int eIdxGlobal = index(element);

            int fIdx = 0;
            // run through all intersections with neighbors
            for (const auto& intersection : intersections(gridView_, element))
            {
                int indexInInside = intersection.indexInInside();
                intersectionMapLocal_[eIdxGlobal][fIdx] = indexInInside;

                fIdx++;
            }
        }
        int globalIntersectionIdx = 0;
        for (const auto& element : elements(gridView_))
        {
            int eIdxGlobal = index(element);

            int fIdx = 0;
            // run through all intersections with neighbors
            for (const auto& intersection : intersections(gridView_, element))
            {
                if (intersection.neighbor())
                {
                    auto neighbor = intersection.outside();
                    int globalIdxNeighbor = index(neighbor);

                    if (element.level() > neighbor.level() || (element.level() == neighbor.level() && eIdxGlobal < globalIdxNeighbor))
                    {

                        int faceIdxNeighbor = 0;
                        if (size(globalIdxNeighbor) == 2 * dim)
                        {
                            faceIdxNeighbor = standardLocalIdxMap_[intersection.indexInOutside()];
                        }
                        else
                        {
                            for (const auto& intersectionNeighbor : intersections(gridView_, neighbor))
                            {
                                if (intersectionNeighbor.neighbor())
                                {
                                    if (intersectionNeighbor.outside() == element)
                                    {
                                        break;
                                    }
                                }
                                faceIdxNeighbor++;
                            }
                        }

                        intersectionMapGlobal_[eIdxGlobal][fIdx] = globalIntersectionIdx;
                        intersectionMapGlobal_[globalIdxNeighbor][faceIdxNeighbor] = globalIntersectionIdx;
                        globalIntersectionIdx ++;
                    }
                }
                else
                {
                    intersectionMapGlobal_[eIdxGlobal][fIdx] = globalIntersectionIdx;
                    globalIntersectionIdx ++;
                }
                fIdx++;
            }
        }
        size_ = globalIntersectionIdx;
    }
    GridView gridView_;
    ElementMapper elementMapper_;
    unsigned int size_;
    std::vector<std::unordered_map<int, int> > intersectionMapGlobal_;
    std::vector<std::unordered_map<int, int> > intersectionMapLocal_;
    std::unordered_map<int, int> standardLocalIdxMap_;
};

} // end namespace Dumux

#endif
