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
 * \ingroup InputOutput
 * \brief A data handle for commucating grid data for gmsh grids
 */
#ifndef DUMUX_GMSH_GRID_DATA_HANDLE_HH
#define DUMUX_GMSH_GRID_DATA_HANDLE_HH

#include <memory>
#include <algorithm>
#include <map>

#include <dune/common/parallel/communication.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/datahandleif.hh>

// UGGrid specific includes
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A data handle for commucating grid data for gmsh grids
 */
template<class Grid, class GridFactory, class Data>
struct GmshGridDataHandle : public Dune::CommDataHandleIF<GmshGridDataHandle<Grid, GridFactory, Data>, typename Data::value_type>
{
    using GridView = typename Grid::LevelGridView;

    GmshGridDataHandle(const Grid& grid, const GridFactory& gridFactory, Data& elementMarkers, Data& boundaryMarkers, Data& faceMarkers)
    : gridView_(grid.levelGridView(0))
    , idSet_(grid.localIdSet())
    , elementMarkers_(elementMarkers)
    , boundaryMarkers_(boundaryMarkers)
    , faceMarkers_(faceMarkers)
    {
        const auto& indexSet = gridView_.indexSet();

        for (const auto& element : elements(gridView_, Dune::Partitions::interior))
           std::swap(elementMarkers_[gridFactory.insertionIndex(element)], data_[idSet_.id(element)]);

        for (const auto& face : entities(gridView_, Dune::Codim<1>()))
           std::swap(faceMarkers_[indexSet.index(face)], data_[idSet_.id(face)]);
    }

    ~GmshGridDataHandle()
    {
        const auto& indexSet = gridView_.indexSet();

        elementMarkers_.resize(indexSet.size(0));
        for (const auto& element : elements(gridView_))
           std::swap(elementMarkers_[indexSet.index(element)], data_[idSet_.id(element)]);

        faceMarkers_.resize(indexSet.size(1));
        for (const auto& face : entities(gridView_, Dune::Codim<1>()))
           std::swap(faceMarkers_[indexSet.index(face)], data_[idSet_.id(face)]);

        boundaryMarkers_.resize(gridView_.grid().numBoundarySegments(), 0);
        for (const auto& element : elements(gridView_.grid().leafGridView()))
        {
            for (const auto& intersection : intersections(gridView_.grid().leafGridView(), element))
            {
                if (intersection.boundary())
                {
                    const auto marker = faceMarkers_[indexSet.index(element.template subEntity<1>(intersection.indexInInside()))];
                    boundaryMarkers_[intersection.boundarySegmentIndex()] = marker;
                }
            }
       }
    }

    Dune::CommDataHandleIF<GmshGridDataHandle<Grid, GridFactory, Data>, typename Data::value_type>& interface()
    { return *this; }

    bool contains (int dim, int codim) const
    { return codim == 0 || codim == 1; }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedSize(int dim, int codim) const
    { return true; }

    template<class EntityType>
    std::size_t size (const EntityType& e) const
    { return 1; }

    template<class MessageBufferImp, class EntityType>
    void gather (MessageBufferImp& buff, const EntityType& e) const
    { buff.write(data_[idSet_.id(e)]); }

    template<class MessageBufferImp, class EntityType>
    void scatter (MessageBufferImp& buff, const EntityType& e, std::size_t n)
    { buff.read(data_[idSet_.id(e)]); }

private:
    using IdSet = typename Grid::LocalIdSet;

    const GridView gridView_;
    const IdSet &idSet_;
    Data& elementMarkers_;
    Data& boundaryMarkers_;
    Data& faceMarkers_;
    mutable std::map< typename IdSet::IdType, typename Data::value_type> data_;
};

#if HAVE_UG

/*!
 * \ingroup InputOutput
 * \brief A data handle for commucating grid data for gmsh grids (specialization for UGGrid)
 */
template<class GridFactory, class Data, int dimgrid>
struct GmshGridDataHandle<Dune::UGGrid<dimgrid>, GridFactory, Data>
: public Dune::CommDataHandleIF<GmshGridDataHandle<Dune::UGGrid<dimgrid>, GridFactory, Data>, typename Data::value_type>
{
    using Grid = Dune::UGGrid<dimgrid>;
    using GridView = typename Grid::LevelGridView;

    GmshGridDataHandle(const Grid& grid, const GridFactory& gridFactory, Data& elementMarkers, Data& boundaryMarkers)
    : gridView_(grid.levelGridView(0))
    , idSet_(grid.localIdSet())
    , elementMarkers_(elementMarkers)
    , boundaryMarkers_(boundaryMarkers)
    {
        for (const auto& element : elements(gridView_, Dune::Partitions::interior))
           std::swap(elementMarkers_[gridFactory.insertionIndex(element)], data_[idSet_.id(element)]);

        // Depending on the Dune version, the boundary markers are present on
        // all processes (<= 2.6) or on the root process only (>= 2.7). Try to
        // handle this in a flexible way: determine if the minimum size over
        // all processes of the boundary markers vector is zero. If yes, assume
        // that the root process contains all markers and broadcast them.
        auto bmSizeMin = boundaryMarkers_.size();
        Dune::MPIHelper::getCollectiveCommunication().min(&bmSizeMin, 1);
        if (bmSizeMin == 0)
        {
            auto bmSize = boundaryMarkers_.size();
            Dune::MPIHelper::getCollectiveCommunication().broadcast(&bmSize, 1, 0);
            boundaryMarkers_.resize(bmSize);
            Dune::MPIHelper::getCollectiveCommunication().broadcast(&boundaryMarkers_.front(), bmSize, 0);
        }
    }

    ~GmshGridDataHandle()
    {
        const auto& indexSet = gridView_.indexSet();
        elementMarkers_.resize(indexSet.size(0));
        for (const auto& element : elements(gridView_))
           std::swap(elementMarkers_[indexSet.index(element)], data_[idSet_.id(element)]);
    }

    Dune::CommDataHandleIF<GmshGridDataHandle<Grid, GridFactory, Data>, typename Data::value_type>& interface()
    { return *this; }

    bool contains (int dim, int codim) const
    { return codim == 0 || codim == 1; }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedSize(int dim, int codim) const
    { return true; }

    template<class EntityType>
    std::size_t size (const EntityType& e) const
    { return 1; }

    template<class MessageBufferImp, class EntityType>
    void gather (MessageBufferImp& buff, const EntityType& e) const
    { buff.write(data_[idSet_.id(e)]); }

    template<class MessageBufferImp, class EntityType>
    void scatter (MessageBufferImp& buff, const EntityType& e, std::size_t n)
    { buff.read(data_[idSet_.id(e)]); }

private:
    using IdSet = typename Grid::LocalIdSet;

    const GridView gridView_;
    const IdSet &idSet_;
    Data& elementMarkers_;
    Data& boundaryMarkers_;
    mutable std::map< typename IdSet::IdType, typename Data::value_type> data_;
};

#endif // HAVE_UG

} // namespace Dumux

#endif
