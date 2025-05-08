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
#include <array>
#include <vector>
#include <string>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <ranges>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/future.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dumux/io/vtk/vtkreader.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A data handle for communicating grid data for VTK grids
 */
template<class Grid, class GridInput, class Data>
struct VtkGridDataHandle
: public Dune::CommDataHandleIF<VtkGridDataHandle<Grid, GridInput, Data>, typename Data::value_type>
{
    using GridView = typename Grid::LevelGridView;

    VtkGridDataHandle(const Grid& grid, const GridInput& gridInput, VTKReader::Data& cellData, VTKReader::Data& pointData)
    : gridView_(grid.levelGridView(0))
    , idSet_(grid.localIdSet())
    , userCellData_(cellData)
    , userPointData_(pointData)
    {
        const auto rank = grid.comm().rank();

        // For the following to work we assume a sorted map of keys to values in the user data.
        // This is not guaranteed by the VTKReader, so we need to sort the data first.
        for (const auto& [key, data] : userCellData_)
            cellData_[key] = std::move(userCellData_[key]);
        for (const auto& [key, data] : userPointData_)
            pointData_[key] = std::move(userPointData_[key]);

#if DUMUX_HAVE_GRIDFORMAT
        // compute the number of components for each cell and point data array
        // and the serialization size per entity
        std::array<std::size_t, 2> numKeys{{ cellData_.size(), pointData_.size() }};
        std::vector<std::size_t> keyComponents(numKeys[0] + numKeys[1], 0);
        {
            int n = 0;
            for (const auto& [key, data] : cellData_)
                keyComponents[n++] = rank == 0 ? data.size()/gridInput.numElements() : 0;
            for (const auto& [key, data] : pointData_)
                keyComponents[n++] = rank == 0 ? data.size()/gridInput.numVertices() : 0;
            grid.comm().broadcast(keyComponents.data(), keyComponents.size(), 0);
        }

        const auto begin = keyComponents.begin();
        cellDataComponents_.assign(begin, begin + numKeys[0]);
        pointDataComponents_.assign(begin + numKeys[0], keyComponents.end());
        numCellDataPerElement_ = std::accumulate(cellDataComponents_.begin(), cellDataComponents_.end(), 0UL);
        numPointDataPerVertex_ = std::accumulate(pointDataComponents_.begin(), pointDataComponents_.end(), 0UL);
#else
        // assume all data is on rank 0 (see grid manager)
        // First broadcast how many keys we have
        std::array<std::size_t, 2> numKeys{{ cellData_.size(), pointData_.size() }};
        grid.comm().broadcast(numKeys.data(), 2, 0);

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
                keyLengthAndComponents[n++] = rank == 0 ? data.size()/gridInput.numElements() : 0;
            for (const auto& [key, data] : pointData_)
                keyLengthAndComponents[n++] = rank == 0 ? data.size()/gridInput.numVertices() : 0;

            grid.comm().broadcast(keyLengthAndComponents.data(), keyLengthAndComponents.size(), 0);
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

            grid.comm().broadcast(keys.data(), keys.size(), 0);
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
#endif

        // write data into an id map
        for (const auto& element : elements(gridView_, Dune::Partitions::interior))
        {
            data_[idSet_.id(element)].resize(numCellDataPerElement_);

            int n = 0, l = 0;
            for (const auto& [key, data] : cellData_)
            {
                const auto nComp = cellDataComponents_[l++];
                for (int k = 0; k < nComp; ++k)
                    std::swap(cellData_[key][k + nComp*gridInput.insertionIndex(element)], data_[idSet_.id(element)][n++]);
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
                    std::swap(pointData_[key][k + nComp*gridInput.insertionIndex(vertex)], data_[idSet_.id(vertex)][n++]);
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

    Dune::CommDataHandleIF<VtkGridDataHandle<Grid, GridInput, Data>, typename Data::value_type>& interface()
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

namespace Dumux::Detail::VtkData {

// for structured vtk data, we manually distribute the data to the other ranks
template<class Grid, class GridInput>
void communicateStructuredVtkData(const Grid& grid, const GridInput& gridInput, ::Dumux::VTKReader::Data& cellData, ::Dumux::VTKReader::Data& pointData)
{
#if HAVE_MPI // needed due to oversight in Dune::Communication interface
    const auto commSize = grid.comm().size();
    if (commSize <= 1)
        return;

    const auto rank = grid.comm().rank();

#ifndef NDEBUG
    if (rank == 0)
    {
        std::cout << "Communicating structured VTK data...\n\n" << std::endl;
        std::cout << "Grid has " << gridInput.numElements() << " elements and " << gridInput.numVertices() << " vertices." << std::endl;
    }
#endif

    // first some preliminary steps
    // we need to sort the data
    std::map<std::string, ::Dumux::VTKReader::Data::mapped_type> sortedCellData, sortedPointData;
    for (const auto& [key, data] : cellData)
        sortedCellData[key] = std::move(cellData[key]);
    for (const auto& [key, data] : pointData)
        sortedPointData[key] = std::move(pointData[key]);

    // and we need to compute the number of components for each cell and point data array
    // to know the message sizes
    std::array<std::size_t, 2> numKeys{{ sortedCellData.size(), sortedPointData.size() }};
    std::vector<std::size_t> keyComponents(numKeys[0] + numKeys[1], 0);
    {
        int n = 0;
        for (const auto& [key, data] : sortedCellData)
            keyComponents[n++] = rank == 0 ? data.size()/gridInput.numElements() : 0;
        for (const auto& [key, data] : sortedPointData)
            keyComponents[n++] = rank == 0 ? data.size()/gridInput.numVertices() : 0;
        grid.comm().broadcast(keyComponents.data(), keyComponents.size(), 0);
    }

    const auto begin = keyComponents.begin();
    const std::size_t numCellDataPerElement = std::accumulate(begin, begin + numKeys[0], 0UL);
    const std::size_t numPointDataPerVertex = std::accumulate(begin + numKeys[0], keyComponents.end(), 0UL);

#ifndef NDEBUG
    std::cout << "Rank " << rank << ": numbers of components for each cell and point data array: "
              << numCellDataPerElement << " " << numPointDataPerVertex << std::endl;
#endif

    // each process knows:
    // - total number of cells and vertices
    // - its data indices for each cell and vertex (global numbering)
    // process 0 has all the data

    // each rank decides which data they need
    const auto& gridView = grid.levelGridView(0);
    std::vector<std::size_t> requestedData(gridView.size(0) + gridView.size(Grid::dimension));
    const std::size_t numElements = gridView.size(0);
    const std::size_t numVertices = gridView.size(Grid::dimension);
    for (const auto& element : elements(gridView))
    {
        requestedData[gridView.indexSet().index(element)] = gridInput.insertionIndex(element);
        for (const auto& vertex : subEntities(element, Dune::Codim<Grid::dimension>{}))
            requestedData[numElements + gridView.indexSet().index(vertex)] = gridInput.insertionIndex(vertex);
    }

    // gather the sizes of the data on rank 0
    std::vector<std::size_t> numData_;
    if (rank == 0)
        numData_.resize(2*commSize);
    std::array<std::size_t, 2> localNumData{{ numElements, numVertices }};
    grid.comm().gather(localNumData.data(), numData_.data(), 2, 0);

#ifndef NDEBUG
    if (rank == 0)
    {
        std::cout << "Number of elements and vertices on each rank: ";
        for (std::size_t i = 0; i < commSize; ++i)
            std::cout << numData_[2*i] << " " << numData_[2*i + 1] << " ";
        std::cout << std::endl;
    }
#endif

    // send the data request to rank 0
    using RequestData = std::vector<std::size_t>;
    using FutureIndices = Dune::Future<RequestData>;
    std::unique_ptr<FutureIndices> sendRequest;
    if (rank != 0)
        sendRequest = std::make_unique<FutureIndices>(
            std::move(grid.comm().isend(std::move(requestedData), 0, /*tag*/0))
        );

    // receive the data on rank 0
    std::vector<RequestData> requestedDataAll(commSize);
    if (rank == 0)
    {
        std::vector<std::unique_ptr<FutureIndices>> receiveRequests(commSize-1);
        for (std::size_t i = 0; i < commSize; ++i)
        {
            requestedDataAll[i].resize(numData_[2*i] + numData_[2*i + 1]);

            if (i == 0)
                requestedDataAll[i] = std::move(requestedData);
            else
            {
                receiveRequests[i-1] = std::make_unique<FutureIndices>(
                    std::move(grid.comm().irecv(requestedDataAll[i], i, /*tag*/0))
                );
            }
        }

        /// TODO: actually we want to call MPI_Waitall here: how to do this with the futures?
        std::ranges::for_each(receiveRequests, [](auto& request) { request->wait(); });

#ifndef NDEBUG
        std::cout << "Received data from all ranks." << std::endl;
        for (std::size_t i = 0; i < commSize; ++i)
            std::cout << "From rank " << i << ": " << requestedDataAll[i].size() << " index size." << std::endl;
#endif
    }

    // Now we need to pack the data and send it to the other ranks
    using PackedData = std::vector<double>;
    PackedData receivedFlatData(localNumData[0]*numCellDataPerElement + localNumData[1]*numPointDataPerVertex);
    if (rank == 0)
    {
        using FutureData = Dune::Future<PackedData>;
        std::vector<std::unique_ptr<FutureData>> sendRequests(commSize-1);
        for (std::size_t i = 0; i < commSize; ++i)
        {
            const auto& requestedIndices = requestedDataAll[i];
            PackedData packedData(numData_[2*i]*numCellDataPerElement + numData_[2*i + 1]*numPointDataPerVertex);

            // pack the data
            std::size_t n = 0, l = 0;
            for (const auto& [key, data] : sortedCellData)
            {
                const auto nComp = keyComponents[l++];
                for (std::size_t k = 0; k < numData_[2*i]; ++k)
                {
                    const auto idx = requestedIndices[k];
                    for (std::size_t j = 0; j < nComp; ++j)
                        packedData[n++] = data[j + nComp*idx];
                }
            }

            const auto pointDataOffsett = numData_[2*i];
            for (const auto& [key, data] : sortedPointData)
            {
                const auto nComp = keyComponents[l++];
                for (std::size_t k = 0; k < numData_[2*i + 1]; ++k)
                {
                    const auto idx = requestedIndices[k + pointDataOffsett];
                    for (std::size_t j = 0; j < nComp; ++j)
                        packedData[n++] = data[j + nComp*idx];
                }
            }

            if (n != packedData.size())
                DUNE_THROW(Dune::InvalidStateException, "Packed data size does not match expected size.");

            if (i == 0)
                receivedFlatData = std::move(packedData);
            else
            {
                // send the data to rank i
                sendRequests[i-1] = std::make_unique<FutureData>(
                    std::move(grid.comm().isend(std::move(packedData), i, /*tag*/1))
                );
            }
        }

        /// TODO: actually we want to call MPI_Waitall here: how to do this with the futures?
        std::ranges::for_each(sendRequests, [](auto& request) { request->wait(); });
    }
    else
    {
        // receive the data from rank 0
        auto receiveRequest = grid.comm().irecv(receivedFlatData, 0, /*tag*/1);

        // finalize all communication
        sendRequest->wait();
        receiveRequest.wait();
    }

#ifndef NDEBUG
    std::cout << "On rank " << rank << ", the received data size is " << receivedFlatData.size() << std::endl;
#endif

    // unpack the data on each process into cellData and pointData
    std::size_t n = 0, l = 0;
    for (const auto& [key, data] : sortedCellData)
    {
        auto& out = cellData[key];
        const auto nComp = keyComponents[l++];
        out.resize(localNumData[0]*nComp);
        for (std::size_t k = 0; k < localNumData[0]; ++k)
            for (std::size_t j = 0; j < nComp; ++j)
                out[j + nComp*k] = receivedFlatData[n++];
    }

    for (const auto& [key, data] : sortedPointData)
    {
        auto& out = pointData[key];
        const auto nComp = keyComponents[l++];
        out.resize(localNumData[1]*nComp);
        for (std::size_t k = 0; k < localNumData[1]; ++k)
            for (std::size_t j = 0; j < nComp; ++j)
                out[j + nComp*k] = receivedFlatData[n++];
    }

    if (n != receivedFlatData.size())
        DUNE_THROW(Dune::InvalidStateException, "Unpacked data size does not match expected size.");

#ifndef NDEBUG
    if (rank == 0)
        std::cout << "\n\n+++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
#endif
#endif
}

} // end namespace Dumux::Detail::VtkData

#endif
