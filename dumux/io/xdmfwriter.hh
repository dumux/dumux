// -*- tab-width: 4; indent-tabs-mode: nil -*-
// TODO check license
// This file is from the BAW dune-swf module

#ifndef DUNE_SWF_XDMFWriter
#define DUNE_SWF_XDMFWriter

extern "C"
{
#include <hdf5.h>
}
#include <vector>
#include <tuple>
#include <numeric>
#include <strstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <numeric>
#include <dune/common/indent.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/function.hh>
namespace Dune {
namespace Swf {
namespace Detail
{
/// \brief A data handle to partition the vertices.
///
/// A vertex ends up in the partition of the process that has
/// the lowest rank of all processes that know this vertex.
/// \tparam GridView The type of the grid view used for communication.
template<class GridView>
class DetermineOwnerVerticesHandle
    : public Dune::CommDataHandleIF<DetermineOwnerVerticesHandle<GridView>,int>
{
public:
    typedef int DataType;

    DetermineOwnerVerticesHandle(int rank, std::vector<int>& isOwner,
                         MultipleCodimMultipleGeomTypeMapper<typename Grid::GridView> mapper (GridView, mcmgVertexLayout())
                         )
        : rank_(rank), isOwner_(isOwner), mapper_(mapper)
    {}
    bool fixedsize(int /* dim */, int /* codim */) const
    {
        return true;
    }
    bool contains(int /* dim */, int  codim) const
    {
        return codim==GridView::dimension;
    }
    template<class T>
    std::size_t size(const T&) const
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T& ) const
    {
        // We send our rank
        buffer.write(rank_);
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t)
    {
        int orank = -1;
        buffer.read(orank);
        // We do not own if we received a lower rank.
        if ( orank < rank_)
            isOwner_[mapper_.index(t)] = false;
    }
private:
    int rank_;
    std::vector<int>& isOwner_;
    MultipleCodimMultipleGeomTypeMapper< GridView, MCMGVertexLayout >& mapper_;
};

template<class GridView, class Id>
class CopyVertexDataHandle
    : public Dune::CommDataHandleIF<CopyVertexDataHandle<GridView,Id>, Id>
{
public:
    typedef std::size_t DataType;

    CopyVertexDataHandle(std::vector<Id>& data,
                         MultipleCodimMultipleGeomTypeMapper< GridView, MCMGVertexLayout >& mapper)
        : data_(data), mapper_(mapper)
    {}
    bool fixedsize(int /* dim */, int /* codim */) const
    {
        return true;
    }
    bool contains(int /* dim */, int  codim) const
    {
        return codim==GridView::dimension;
    }
    template<class T>
    std::size_t size(const T&) const
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t) const
    {
        // We send our rank
        buffer.write(data_[mapper_.index(t)]);
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t)
    {
        DataType d;
        buffer.read(d);
        if ( d < std::numeric_limits<DataType>::max() )
            data_[mapper_.index(t)]=d;
    }
private:
    std::vector<Id>& data_;
    MultipleCodimMultipleGeomTypeMapper< GridView, MCMGVertexLayout >& mapper_;
};

} // end namespace Detail
template<class GridView>
class XDMFWriter
{

    typedef typename GridView::Grid Grid;
    typedef typename GridView::ctype CoordinateType;
    static const int dimension = GridView::dimension;
    static const int dimensionworld = GridView::dimensionworld;

    typedef MultipleCodimMultipleGeomTypeMapper< GridView, MCMGVertexLayout > VertexMapper;
    static const PartitionIteratorType HDF5PartitionType = InteriorBorder_Partition;
public:
    typedef Dune::VTKFunction< GridView > VTKFunction;

    /**
     * @brief Creates an XDMFWriter.
     *
     * This is a collective operation that has to be called on all ranks
     * of the communicator used by the grid.
     * @param filename The filename for the XDMF file. The HDF5 file will be
     * named accordingly. If not present the extension xdmf will be added.
     * @param gridView The grid view that the data is attached on.
     * @param description Arbitrary string describing the data file.
     * @param root_rank The rank that writes the XDMF file stub (default 0).
     */
    XDMFWriter(const std::string& filename, const GridView& gridView,
               const std::string& description = std::string(),
               int root_rank=0);

    ~XDMFWriter();

    /**
     * @brief Begin writing a new time step with a different grid.
     *
     * Use this function only if the last grid specified with either the
     * Constructor or this method is not valid any more.
     * Otherwise the grid will be written to file another time.
     * All subsequent calls
     * to beginTimeStep(double) will assume that the same grid will be used.
     * @param time the current time
     */
    void beginTimeStep(double time, const GridView& gridView);

    /**
     * @brief Begin writing a new time step using the last grid.
     *
     * The last recorded grid will be ued to attache the data to.
     */
    void beginTimeStep(double time);

    void endTimeStep();

    /**
     * @brief Write cell data to the HDF5 file and record it in the xdmf file.
     * @param p The VTK function describing the data to be written.
     * @param unit A string describing the units (optional).
     */
    void writeCellData(const std::shared_ptr< const VTKFunction > & p,
                       const std::string& unit = std::string());


    /**
     * @brief Write cell data to the HDF5 file and record it in the xdmf file.
     * @param v Container describing containing the cell data.
     * @param name The name of the data.
     * @param ncomps The number of components per cell.
     * @param unit A string describing the units (optional).
     */
    template<class Container>
    void writeCellData(const Container& v, const std::string& name,
                       const std::string& unit=std::string(),
                       // Use decltype to prevent compiler of using this
                       // function with std::shared<const VTKFunktion>
                       decltype(v.size()) ncomps=1);

protected:
    std::tuple<std::size_t,std::size_t> precomputeVerticesAndCells(const GridView& gridView);
    void writeTopologyAndGeometry(const GridView& gridView);
    std::string getCoordinateType(const std::size_t& dim) const;
    void finishAndFlushXDMFFile();
    void writeStringAttribute(const auto& dataSet, const std::string& name,
                         const std::string& value);
    typedef typename GridView::template Codim< 0 >
    ::template Partition< HDF5PartitionType >::Iterator
    GridCellIterator;
    typedef typename GridView::template Codim< dimension >
    ::template Partition< HDF5PartitionType >::Iterator
    GridVertexIterator;

    typedef typename GridCellIterator::Reference EntityReference;

    struct CellIterator : public GridCellIterator
    {
        CellIterator(const GridCellIterator & it) : GridCellIterator(it) {}
        //! get the position of the center of the element, in element-local
        //! coordinates
        const FieldVector<CoordinateType,dimension> position() const
        {
            return ReferenceElements<CoordinateType,dimension>::general((*this)->type()).position(0,0);
        }
    };

    CellIterator cellBegin() const
    {
      return gridView_.template begin< 0, HDF5PartitionType >();
    }

    CellIterator cellEnd() const
    {
      return gridView_.template end< 0, HDF5PartitionType >();
    }

private:
    GridView gridView_;
    int gridIndex_;
    int timeStepIndex_;
    int dataIndex_;
    std::vector<int> uniqueVertexIndices_;
    std::ofstream xdmfFile_;
    std::string hdf5Name_;
    std::streampos xdmfFilePos_;
    Indent indent_;
    hid_t hdf5PropList_;
    hid_t hdf5FileId_;
    hid_t gridsId_;
    hid_t timeStepsId_;
    hid_t currentTimeStepId_;
    std::size_t noCells_;
    std::size_t noVertices_;
    std::vector<std::size_t> cellOffsets_;

};

template<class GridView>
XDMFWriter<GridView>::XDMFWriter(const std::string& filename,
                                 const GridView& gridView,
                                 const std::string& description,
                                 int root_rank)
    : gridView_(gridView), gridIndex_(-1), timeStepIndex_(0), dataIndex_(0)
{
    size_t dotIndex = filename.find_last_of(".");
    if ( dotIndex == std::string::npos )
        dotIndex = filename.length();

    std::ostringstream hdf5NameStr;
    std::ostrstream xdmfNameStr;
    hdf5NameStr << filename.substr(0, dotIndex);
    xdmfNameStr << filename.substr(0, dotIndex);
    xdmfNameStr << ".xdmf"<<std::ends;
    hdf5NameStr << ".hdf5";
    hdf5Name_ = hdf5NameStr.str();

    MPI_Info info =  MPI_INFO_NULL;

//----- Lustre MPI-IO problems https://wickie.hlrs.de/platforms/index.php/MPI-IO
/*  MPI_Info info;
    MPI_Info_create(&info);

    // Disables ROMIO's data-sieving
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");

    // Enable ROMIO's collective buffering
    MPI_Info_set(info, "romio_cb_read", "enable");
    MPI_Info_set(info, "romio_cb_write", "enable");
*/

//----end Lustre MPI-IO problem

    hdf5PropList_ = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(hdf5PropList_, gridView.comm(), info);
    hdf5FileId_ = H5Fcreate(hdf5Name_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, hdf5PropList_);
    // Write the description if it is not empty
    if( ! description.empty() )
        writeStringAttribute(hdf5FileId_, std::string("description"), description);

    // Create grids and timesteps directories
    gridsId_ = H5Gcreate(hdf5FileId_, "/grids", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    timeStepsId_ = H5Gcreate(hdf5FileId_ , "/timesteps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if( gridView.grid().comm().rank() == root_rank)
    {
        xdmfFile_.open(xdmfNameStr.str());
    }

    writeTopologyAndGeometry(gridView);
    // Write closing tags to XDMF file
    finishAndFlushXDMFFile();

}

template<class GridView>
void XDMFWriter<GridView>::finishAndFlushXDMFFile()
{
    if( xdmfFile_.is_open() )
    {
        xdmfFilePos_ = xdmfFile_.tellp();
        auto oldIndent = indent_;
        xdmfFile_<<--indent_<<"</Domain>"<<std::endl;
        xdmfFile_<<--indent_<<"</Xdmf>"<<std::endl<<std::flush;
        indent_ = oldIndent;
    }
}

template<class GridView>
void XDMFWriter<GridView>::writeStringAttribute(const auto& dataSet, const std::string& name,
                         const std::string& value)
{
    auto stringDataType = H5Tcopy(H5T_C_S1 );
    H5Tset_size(stringDataType, H5T_VARIABLE);
    hsize_t dims[1] = { 1 };
    const char *strings[1] = { value.c_str() };
    auto nameId = H5Screate_simple (1, dims, NULL);
    auto attributeId = H5Acreate (dataSet, name.c_str(), stringDataType, nameId, H5P_DEFAULT, H5P_DEFAULT);
    auto status = H5Awrite(attributeId, stringDataType, strings);
    // Cleanup
    H5Tclose (stringDataType);
    H5Aclose(attributeId);
    H5Sclose(nameId);
}

template<class GridView>
void XDMFWriter<GridView>::beginTimeStep(double time, const GridView& gridView)
{
    gridView_ = gridView;
    writeTopologyAndGeometry(gridView);
    beginTimeStep(time);
}

template<class GridView>
void XDMFWriter<GridView>::beginTimeStep(double time)
{
    std::ostrstream timestr;
    timestr << "timestep_" <<timeStepIndex_<<std::ends;

    if( xdmfFile_.is_open() )
    {
        xdmfFile_.seekp(xdmfFilePos_);
        xdmfFile_<<indent_<<"<Grid Name=\""<<timestr.str()<<"\">"<<std::endl;
        xdmfFile_<<++indent_<<"<Time Value=\""<<time<<"\" />"<<std::endl;
        xdmfFile_<<indent_<<"<Topology Reference=\"/Xdmf/Domain/Topology["<<gridIndex_+1<<"]\" />"<<std::endl;
        xdmfFile_<<indent_<<"<Geometry Reference=\"/Xdmf/Domain/Geometry["<<gridIndex_+1<<"]\" />"<<std::endl;
        xdmfFilePos_ = xdmfFile_.tellp();
    }
    dataIndex_ = 0;
    currentTimeStepId_ = H5Gcreate(timeStepsId_, timestr.str(),
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Write grid name as attribute
    std::ostringstream gridName;
    gridName << "grid_" << gridIndex_;
    writeStringAttribute(currentTimeStepId_, std::string("grid"), gridName.str().c_str());
    hsize_t dims[1] = { 1 };
    auto timeId = H5Screate_simple(1, dims, NULL);
    auto attributeId = H5Acreate(currentTimeStepId_, "time", H5T_NATIVE_DOUBLE, timeId, H5P_DEFAULT, H5P_DEFAULT);
    auto status = H5Awrite(attributeId, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(attributeId);
    H5Sclose(timeId);

}

template<class GridView>
void XDMFWriter<GridView>::endTimeStep()
{
    H5Gclose(currentTimeStepId_);
    if( xdmfFile_.is_open() )
    {
        xdmfFile_<<--indent_<<"</Grid>"<<std::endl;
        xdmfFilePos_ = xdmfFile_.tellp();
    }
    ++timeStepIndex_;
    finishAndFlushXDMFFile();
}

template<class GridView>
std::string XDMFWriter<GridView>::getCoordinateType(const std::size_t& dim) const
{
    assert(dim<4);
    static std::array<std::string, 3> types =
        { "X", "XY", "XYZ"};
    return types[dim-1];
}

template<class GridView>
void XDMFWriter<GridView>::writeTopologyAndGeometry(const GridView& gridView)
{
    ++gridIndex_;
    std::tie(noCells_, noVertices_) =
        precomputeVerticesAndCells(gridView);
    if( xdmfFile_.is_open() )
    {
        xdmfFile_<<"<?xml version=\"1.0\" ?>"<<std::endl;
        xdmfFile_<<"<Xdmf Version=\"2.0\">"<<std::endl;
        xdmfFile_<<++indent_<<"<Domain>"<<std::endl;
        xdmfFile_<<++indent_<<"<Geometry GeometryType=\""
                 <<getCoordinateType(dimensionworld)<<"\">"<<std::endl;
        xdmfFile_<<++indent_<<"<DataItem Dimensions=\""<<noVertices_<<" "
                 <<dimensionworld<<"\" Format=\"HDF\" "
                 <<"NumberType=\"Float\" Precision=\"8\">"
                 <<hdf5Name_<<":/grids/grid_"<<gridIndex_<<"/vertices</DataItem>"
                 <<std::endl;
        xdmfFile_<<--indent_<<"</Geometry>"<<std::endl;
        xdmfFile_<<indent_<<"<Topology TopologyType=\"Triangle\" NumberOfElements=\""
                 <<noCells_<<"\">"<<std::endl;
        xdmfFile_<<++indent_<<"<DataItem Dimensions=\""<< (3*noCells_)
                 <<"\" Format=\"HDF\" NumberType=\"Int\" Precision=\"8\">"
                 <<hdf5Name_<<":/grids/grid_"<<gridIndex_<<"/topology</DataItem>"
                 <<std::endl;
        xdmfFile_<<--indent_<<"</Topology>"<<std::endl;
    }
}

template<class GridView>
std::tuple<std::size_t,std::size_t> XDMFWriter<GridView>::precomputeVerticesAndCells(const GridView& gridView)
{
    /*    if ( gridView.comm().size() == 1 )
    {
        return std::make_tuple(gridView_.size(0), gridView_.size(1));
    }
    else*/
    {
        std::size_t cells=0;
        using GlobalId = std::size_t; //GridView::Traits::GlobalId;
        std::vector<GlobalId> vertexGlobalIndex(gridView.size(dimension),
                                                std::numeric_limits<GlobalId>::max());
        std::vector<int> vertexOwned(gridView.size(dimension), false);
        VertexMapper vertexMapper(gridView);

        for(auto vertex : vertices(gridView) )
        {
            if( vertex.partitionType() == Dune :: InteriorEntity ||
                vertex.partitionType() == Dune :: BorderEntity)
                vertexOwned[vertexMapper.index(vertex)] = true;
        }

        Detail::DetermineOwnerVerticesHandle<GridView>
            handle(gridView.comm().rank(), vertexOwned, vertexMapper);
        gridView.communicate(handle, Dune::InteriorBorder_All_Interface,
                             Dune::ForwardCommunication);
        std::size_t ownedVertices = std::count_if(vertexOwned.begin(),
                                                  vertexOwned.end(),
                                                  [](bool a){ return a; });
        // Compute the offsets for the parallel write operation.
        std::vector<GlobalId> offsets(gridView.comm().size()+1, 0);
        GlobalId noVertices = ownedVertices;
        gridView.comm().allgather(&noVertices, 1, offsets.data()+1);
        std::partial_sum(offsets.begin()+1, offsets.end(), offsets.begin()+1);

        // The index of our first vertex.
        GlobalId vertexIndex = offsets[gridView.comm().rank()];

        // Create the coordinate array to write and renumber the vertices.
        std::vector<CoordinateType> coordinates(dimension*ownedVertices);
        std::vector<int> notVisited(gridView.size(dimension), true);
        auto currentCoordinate = coordinates.begin();

        for(auto cell = gridView_.template begin< 0, HDF5PartitionType >(),
                endCell = gridView_.template end< 0, HDF5PartitionType >();
            cell != endCell; ++cell)
        {
            ++cells;
            for(int corner = 0, corners = cell->subEntities(GridView::dimension);
                corner < corners; ++corner)
            {
                auto index = vertexMapper.subIndex(*cell, corner, GridView::dimension);
                if ( vertexOwned[index] && notVisited[index] )
                {
                    notVisited[index] = false;
                    vertexGlobalIndex[index] = vertexIndex;
                    const auto& coord = cell->geometry().corner(corner);
                    assert(currentCoordinate+dimensionworld<=coordinates.end());
                    std::copy(coord.begin(), coord.end(), currentCoordinate);
                    currentCoordinate += dimensionworld;
                    ++vertexIndex;
                }
            }
        }
        assert(vertexIndex == offsets[gridView.comm().rank()+1]);
        Detail::CopyVertexDataHandle<GridView, GlobalId>
            copyHandle(vertexGlobalIndex, vertexMapper);
        gridView.communicate(copyHandle, Dune::InteriorBorder_InteriorBorder_Interface,
                             Dune::ForwardCommunication);

        // Write the coordinates to file.
        std::ostrstream tmpString;
        tmpString << "grid_" << gridIndex_<<std::ends;
        auto gridName = tmpString.str(); // for reuse later.
        auto gridId = H5Gcreate(gridsId_, gridName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // the global dimensions
        auto mpi_size = gridView.comm().size();
        auto mpi_rank = gridView.comm().rank();
        hsize_t globalDims[2] = { offsets[mpi_size],
                                  dimensionworld };

        //using DataType = HDF5Datatype<CoordinateType>::type;
        //using DataType = H5T_NATIVE_DOUBLE;
        auto fileSpace = H5Screate_simple(2, globalDims, NULL);
        auto dataSetId = H5Dcreate(gridId, "vertices",
                                   H5T_NATIVE_DOUBLE, fileSpace,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(fileSpace);

        hsize_t localDims[2] = { offsets[mpi_rank+1] - offsets[mpi_rank], dimensionworld };
        hsize_t localOffsets[2] = { offsets[mpi_rank], 0};

        auto memorySpace = H5Screate_simple(2, localDims, NULL);
        // Select hyperslab in the file
        fileSpace = H5Dget_space(dataSetId);
        H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, localOffsets, NULL,
                            localDims, NULL);

        // Create property list for collective dataset write
        auto propertyListId = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(propertyListId, H5FD_MPIO_COLLECTIVE);

        auto status = H5Dwrite(dataSetId, H5T_NATIVE_DOUBLE, memorySpace, fileSpace,
                               propertyListId, coordinates.data());
        // Cleanup
        H5Dclose(dataSetId);
        H5Sclose(fileSpace);
        H5Sclose(memorySpace);
        //H5Pclose(propertyListId);
        // Create the topology array to write.
        // Currently only triangles are supported.
        std::vector<std::size_t> topology(3*cells);
        vertexIndex = 0;

        for(auto cell = gridView_.template begin< 0, HDF5PartitionType >(),
                endCell = gridView_.template end< 0, HDF5PartitionType >();
            cell != endCell; ++cell)
        {
           for(int corner = 0, corners = cell->subEntities(GridView::dimension);
                corner < corners; ++corner)
            {
                auto index = vertexMapper.subIndex(*cell, corner, GridView::dimension);
                topology[vertexIndex++] = vertexGlobalIndex[index];
            }
        }

        cellOffsets_.resize(gridView.comm().size()+1, 0);
        GlobalId ownedCells = cells;
        gridView.comm().allgather(&ownedCells, 1, cellOffsets_.data()+1);
        std::partial_sum(cellOffsets_.begin()+1, cellOffsets_.end(), cellOffsets_.begin()+1);

        hsize_t globalTopologyDims[1] = { 3 * cellOffsets_[mpi_size] };
        fileSpace = H5Screate_simple(1, globalTopologyDims, NULL);
        dataSetId = H5Dcreate(gridId, "topology",
                              H5T_NATIVE_ULONG, fileSpace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(fileSpace);

        hsize_t localTopologyDims[1] = { 3 * (cellOffsets_[mpi_rank+1] - cellOffsets_[mpi_rank]) };
        hsize_t localTopologyOffsets[1] = { 3 * cellOffsets_[mpi_rank] };

        memorySpace = H5Screate_simple(1, localTopologyDims, NULL);
        // Select hyperslab in the file
        fileSpace = H5Dget_space(dataSetId);
        H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, localTopologyOffsets, NULL,
                            localTopologyDims, NULL);

        // Create property list for collective dataset write
        //auto propertyListId = H5Pcreate(H5P_DATASET_XFER);
        //H5Pset_dxpl_mpio(propertyListId, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dataSetId, H5T_NATIVE_ULONG, memorySpace, fileSpace,
                               propertyListId, topology.data());
        hsize_t numberDims[1] = { 1 };
        auto numberId = H5Screate_simple(1, numberDims, NULL);
        auto attributeId = H5Acreate(dataSetId, "elements", H5T_NATIVE_ULONG, numberId, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite(attributeId, H5T_NATIVE_ULONG, &cellOffsets_[mpi_size]);
        writeStringAttribute(dataSetId, std::string("type"), std::string("Triangle"));
        // Cleanup
        H5Aclose(attributeId);
        H5Sclose(numberId);
        H5Dclose(dataSetId);
        H5Sclose(fileSpace);
        H5Sclose(memorySpace);
        H5Pclose(propertyListId);
        H5Gclose(gridId);
        // communicate vertex indices at border
        return std::tuple<std::size_t,std::size_t>(gridView.comm().sum(cells),
                                                   offsets[mpi_size]);
    }
}

template<class GridView>
void XDMFWriter<GridView>::writeCellData(const std::shared_ptr< const VTKFunction > & func,
                                         const std::string& unit)
{
    std::string center="Cell";
    if( xdmfFile_.is_open() )
    {
        xdmfFile_<<indent_<<"<Attribute AttributeType=\"";
        if ( func->ncomps()==1 )
            xdmfFile_<<"Scalar";
        else
            xdmfFile_<<"Vector";
        xdmfFile_<<"\" "<<"Center=\""<<center<<"\" Name=\""<<func->name()<<"\">"<<std::endl;
        xdmfFile_<<++indent_<<"<DataItem Dimensions=\""<<noCells_;

        if(func->ncomps()>1) xdmfFile_<<" "<<func->ncomps();

        xdmfFile_<<"\" Format=\"HDF\" NumberType=\"Float\" "
                 <<"Precision=\"8\">"<<hdf5Name_<<":/timesteps/timestep_"
                 <<timeStepIndex_<<"/"<<func->name()<<"</DataItem>"<<std::endl;
        xdmfFile_<<--indent_<<"</Attribute>"<<std::endl;
    }
    auto mpi_size = gridView_.comm().size();
    auto mpi_rank = gridView_.comm().rank();
    hsize_t globalDims[2] = { cellOffsets_[mpi_size], static_cast<int>(func->ncomps())  };
    hsize_t localDims[2] = { cellOffsets_[mpi_rank+1] - cellOffsets_[mpi_rank],
                             static_cast<int>(func->ncomps()) };
    hsize_t offsets[2] = { cellOffsets_[mpi_rank], 0 };
    int dim = func->ncomps() <= 1 ? 1 : 2;
    auto fileSpace = H5Screate_simple(dim, globalDims, NULL);
    auto dataSetId = H5Dcreate(currentTimeStepId_, func->name().c_str(),
                              H5T_NATIVE_DOUBLE, fileSpace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(fileSpace);
    auto memorySpace = H5Screate_simple(dim, localDims, NULL);
    // If there is a unit string write it to the HDF5 file
    if( ! unit.empty() )
        writeStringAttribute(dataSetId, std::string("unit"), unit);
    // Select hyperslab in the file
    fileSpace = H5Dget_space(dataSetId);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offsets, NULL,
                        localDims, NULL);
    // Create property list for collective dataset write
    auto propertyListId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(propertyListId, H5FD_MPIO_COLLECTIVE);
    // retrieve data
    std::vector<double> data(localDims[0] * func->ncomps(), -100+mpi_rank);
    auto entry = data.begin();

    for(auto cell = cellBegin(), endCell = cellEnd(); cell != endCell; ++cell)
    {
        // Retrieve element data
        for ( int i = 0; i < func->ncomps(); ++i)
            *(entry++) = func->evaluate(i, *cell, cell.position());
    }
    auto status = H5Dwrite(dataSetId, H5T_NATIVE_DOUBLE, memorySpace, fileSpace,
                               propertyListId, data.data());
    // Cleanup
    H5Dclose(dataSetId);
    H5Sclose(fileSpace);
    H5Sclose(memorySpace);
    H5Pclose(propertyListId);
}

template<class GridView>
template<class Container>
void XDMFWriter<GridView>::writeCellData(const Container& v, const std::string& name,
                                         const std::string& unit,
                                         decltype(v.size()) ncomps)
{
      typedef P0VTKFunction<GridView, Container> Function;
      for (decltype(v.size()) c=0; c<ncomps; ++c) {
        std::stringstream compName;
        compName << name;
        if (ncomps>1)
          compName << "[" << c << "]";
        VTKFunction* p = new Function(gridView_, v, compName.str(), ncomps, c);
        writeCellData(std::shared_ptr< const VTKFunction >(p), unit);
      }
}

template<class GridView>
XDMFWriter<GridView>::~XDMFWriter()
{
    H5Pclose(hdf5PropList_);
    H5Gclose(gridsId_);
    H5Gclose(timeStepsId_);
    H5Fclose(hdf5FileId_);
    xdmfFile_.close();
}

} // end namespace Swf
} // namespace Dune
#endif //DUNE_SWF_XDMFWriter
