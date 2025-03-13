// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief A data handle for communicating grid data for VTK grids
 */
#ifndef DUMUX_VTK_GRID_DATA_HANDLE_HH
#define DUMUX_VTK_GRID_DATA_HANDLE_HH

#include <memory>
#include <algorithm>
#include <map>
#include <numeric>

#include <dune/common/parallel/communication.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dumux/io/vtk/vtkreader.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A data handle for communicating grid data for VTK grids
 */
template<class Grid, class GridFactory, class Data>
struct VtkGridDataHandle
: public Dune::CommDataHandleIF<VtkGridDataHandle<Grid, GridFactory, Data>, typename Data::value_type>
{
    using GridView = typename Grid::LevelGridView;

    VtkGridDataHandle(const Grid& grid, const GridFactory& gridFactory, VTKReader::Data& cellData, VTKReader::Data& pointData)
    : gridView_(grid.levelGridView(0))
    , idSet_(grid.localIdSet())
    , userCellData_(cellData)
    , userPointData_(pointData)
    {
        // For the following to work we assume a sorted map of keys to values in the user data.
        // This is not guaranteed by the VTKReader, so we need to sort the data first.
        for (const auto& [key, data] : userCellData_)
            cellData_[key] = std::move(userCellData_[key]);
        for (const auto& [key, data] : userPointData_)
            pointData_[key] = std::move(userPointData_[key]);

        // assume all data is on rank 0 (see grid manager)
        // First broadcast how many keys we have
        std::array<std::size_t, 2> numKeys{{ cellData_.size(), pointData_.size() }};
        Dune::MPIHelper::getCommunication().broadcast(numKeys.data(), 2, 0);

        // then broadcast the length of the individual key strings
        // and the number of component associated with each key (e.g. vector/tensor fields)
        std::vector<std::size_t> keyLengthAndComponents(2*(numKeys[0] + numKeys[1]), 0);
        {
            int n = 0;

            // key length
            for (const auto& [key, data] : cellData_)
                keyLengthAndComponents[n++] = key.size();
            for (const auto& [key, data] : pointData_)
                keyLengthAndComponents[n++] = key.size();

            // number of components
            for (const auto& [key, data] : cellData_)
                keyLengthAndComponents[n++] = gridView_.size(0) > 0 ? data.size()/gridView_.size(0) : 0;
            for (const auto& [key, data] : pointData_)
                keyLengthAndComponents[n++] = gridView_.size(Grid::dimension) > 0 ? data.size()/gridView_.size(Grid::dimension) : 0;

            // entries only exist on rank 0 and the data containers are empty on other ranks
            assert((Dune::MPIHelper::instance().rank() == 0) == (n == keyLengthAndComponents.size()));

            Dune::MPIHelper::getCommunication().broadcast(keyLengthAndComponents.data(), keyLengthAndComponents.size(), 0);
        }

        // save the number of components for each cell and point data array
        const auto begin = keyLengthAndComponents.begin() + numKeys[0] + numKeys[1];
        cellDataComponents_.assign(begin, begin + numKeys[0]);
        pointDataComponents_.assign(begin + numKeys[0], keyLengthAndComponents.end());
        numCellDataPerElement_ = std::accumulate(cellDataComponents_.begin(), cellDataComponents_.end(), 0UL);
        numPointDataPerVertex_ = std::accumulate(pointDataComponents_.begin(), pointDataComponents_.end(), 0UL);

        // then broadcast the actual keys
        std::string keys; keys.resize(std::accumulate(keyLengthAndComponents.begin(), begin, 0UL));
        {
            int n = 0;
            for (const auto& [key, data] : cellData_)
                for (const auto& c : key)
                    keys[n++] = c;
            for (const auto& [key, data] : pointData_)
                for (const auto& c : key)
                    keys[n++] = c;

            // entries only exist on rank 0 and the data containers are empty on other ranks
            assert((Dune::MPIHelper::instance().rank() == 0) == (n == keys.size()));

            Dune::MPIHelper::getCommunication().broadcast(keys.data(), keys.size(), 0);
        }

        // create the entries in the cellData and pointData maps on all processes
        std::size_t offset = 0;
        for (int keyIdx = 0; keyIdx < numKeys[0]; ++keyIdx)
        {
            if (std::string key{ keys, offset, keyLengthAndComponents[keyIdx] }; cellData_.count(key) == 0)
                cellData_[key] = Data{};

            offset += keyLengthAndComponents[keyIdx];
        }
        for (int keyIdx = numKeys[0]; keyIdx < numKeys[0] + numKeys[1]; ++keyIdx)
        {
            if (std::string key{ keys, offset, keyLengthAndComponents[keyIdx] }; pointData_.count(key) == 0)
                pointData_[key] = Data{};

            offset += keyLengthAndComponents[keyIdx];
        }

        // write data into an id map
        for (const auto& element : elements(gridView_, Dune::Partitions::interior))
        {
            data_[idSet_.id(element)].resize(numCellDataPerElement_);

            int n = 0, l = 0;
            for (const auto& [key, data] : cellData_)
            {
                const auto nComp = cellDataComponents_[l++];
                for (int k = 0; k < nComp; ++k)
                    std::swap(cellData_[key][k + nComp*gridFactory.insertionIndex(element)], data_[idSet_.id(element)][n++]);
            }

            assert(n == numCellDataPerElement_);
        }

        for (const auto& vertex : vertices(gridView_))
        {
            data_[idSet_.id(vertex)].resize(numPointDataPerVertex_);

            int n = 0, l = 0;
            for (const auto& [key, data] : pointData_)
            {
                const auto nComp = pointDataComponents_[l++];
                for (int k = 0; k < nComp; ++k)
                    std::swap(pointData_[key][k + nComp*gridFactory.insertionIndex(vertex)], data_[idSet_.id(vertex)][n++]);
            }

            assert(n == numPointDataPerVertex_);
        }
    }

    ~VtkGridDataHandle()
    {
        // resize arrays and unpack communicated data
        const auto& indexSet = gridView_.indexSet();

        {
            int n = 0;
            for (const auto& [key, data] : cellData_)
                cellData_[key].resize(indexSet.size(0)*cellDataComponents_[n++]);
        }
        {
            int n = 0;
            for (const auto& [key, data] : pointData_)
                pointData_[key].resize(indexSet.size(GridView::dimension)*pointDataComponents_[n++]);
        }

        for (const auto& element : elements(gridView_))
        {
            int n = 0, l = 0;
            for (const auto& [key, data] : cellData_)
            {
                const auto nComp = cellDataComponents_[l++];
                for (int k = 0; k < nComp; ++k)
                    std::swap(cellData_[key][k + nComp*indexSet.index(element)], data_[idSet_.id(element)][n++]);
            }
        }

        for (const auto& vertex : vertices(gridView_))
        {
            int n = 0, l = 0;
            for (const auto& [key, data] : pointData_)
            {
                const auto nComp = pointDataComponents_[l++];
                for (int k = 0; k < nComp; ++k)
                    std::swap(pointData_[key][k + nComp*indexSet.index(vertex)], data_[idSet_.id(vertex)][n++]);
            }
        }

        // move data back from sorted internal storage to user data
        for (const auto& [key, data] : cellData_)
            userCellData_[key] = std::move(cellData_[key]);
        for (const auto& [key, data] : pointData_)
            userPointData_[key] = std::move(pointData_[key]);
    }

    Dune::CommDataHandleIF<VtkGridDataHandle<Grid, GridFactory, Data>, typename Data::value_type>& interface()
    { return *this; }

    bool contains (int dim, int codim) const
    { return codim == 0 || codim == dim; }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedSize(int dim, int codim) const
    { return true; }

    template<class Entity>
    std::size_t size (const Entity&) const
    {
        if constexpr (Entity::codimension == 0)
            return numCellDataPerElement_;
        else if constexpr (Entity::codimension == GridView::dimension)
            return numPointDataPerVertex_;
        else
            return 0;
    }

    template<class MessageBufferImp, class Entity>
    void gather (MessageBufferImp& buff, const Entity& e) const
    {
        if constexpr (Entity::codimension == 0)
        {
            const auto& data = data_[idSet_.id(e)];
            for (int n = 0; n < numCellDataPerElement_; ++n)
                buff.write(data[n]);
        }

        if constexpr (Entity::codimension == GridView::dimension)
        {
            const auto& data = data_[idSet_.id(e)];
            for (int n = 0; n < numPointDataPerVertex_; ++n)
                buff.write(data[n]);
        }
    }

    template<class MessageBufferImp, class Entity>
    void scatter (MessageBufferImp& buff, const Entity& e, std::size_t n)
    {
        auto& data = data_[idSet_.id(e)];
        data.resize(n);

        if constexpr (Entity::codimension == 0)
            for (int k = 0; k < numCellDataPerElement_; ++k)
                buff.read(data[k]);

        if constexpr (Entity::codimension == GridView::dimension)
            for (int k = 0; k < numPointDataPerVertex_; ++k)
                buff.read(data[k]);
    }

private:
    using IdSet = typename Grid::LocalIdSet;

    const GridView gridView_;
    const IdSet &idSet_;

    VTKReader::Data& userCellData_;
    VTKReader::Data& userPointData_;

    std::map<std::string, VTKReader::Data::mapped_type> cellData_;
    std::map<std::string, VTKReader::Data::mapped_type> pointData_;

    std::vector<std::size_t> cellDataComponents_;
    std::vector<std::size_t> pointDataComponents_;

    std::size_t numCellDataPerElement_;
    std::size_t numPointDataPerVertex_;

    mutable std::map< typename IdSet::IdType, std::vector<typename Data::value_type> > data_;
};

} // namespace Dumux

#endif
