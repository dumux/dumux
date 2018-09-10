#ifndef DUMUX_HDF5GRIDREADER_HH
#define DUMUX_HDF5GRIDREADER_HH

extern "C"
{
#include <hdf5.h>
}

#include <string>
#include <strstream>
#include <vector>
#include <numeric>
#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/vtk.hh>

namespace Dune
{
namespace Dumux
{
namespace Detail
{
/// \brief A datahandler to distribute the cell data during loadBalance.
template<class Grid>
class HDF5ReaderHandle
    : public Dune::CommDataHandleIF<HDF5ReaderHandle<Grid>,double>
{
public:
    typedef double DataType;
    HDF5ReaderHandle(std::vector<std::string>&& names,
                     const std::vector<std::vector<double>>&& dataArrays,
                     std::map<std::string, std::vector<double>>& distributedData,
                     const GridFactory<Grid>& factory,
                     const Grid& grid)
        : grid_(grid), names_(names), distributedData_(distributedData),
          factory_(factory), isRoot_(grid_.comm().rank() == 0)
    {
        // Broadcast the names of the data arrays to all processes
        std::size_t size=names_.size();
        std::size_t globalCells = grid_.size(0);
        grid_.comm().broadcast(&size, 1, 0);
        grid_.comm().broadcast(&globalCells, 1, 0);
        std::size_t bufferLength=0;
        multipliers_.resize(size);
        if ( isRoot_ )
        {
            auto dataItem = dataArrays.begin();
            auto multiplier = multipliers_.begin();

            for (const auto& name : names_)
            {
                bufferLength += name.length() + 1; //Beware of trailing '\0'
                *multiplier = dataItem->size()/globalCells;
                ++multiplier; ++dataItem;
            }
        }

        grid_.comm().broadcast(&bufferLength, 1, 0);
        grid_.comm().broadcast(multipliers_.data(), size, 0);
        elementDataSize_ = std::accumulate(multipliers_.begin(),
                                           multipliers_.end(), 0);
        char* namesBuffer = new char[bufferLength];
        if ( isRoot_ )
        {
            char* currentChar = namesBuffer;
            for (const auto& name : names_)
            {
                currentChar = std::copy(name.begin(), name.end(), currentChar);
                *(currentChar++) = '\0';
            }
            assert(currentChar - namesBuffer == bufferLength);
        }
        grid_.comm().broadcast(namesBuffer, bufferLength, 0);
        if ( ! isRoot_ )
        {
            // initialize the names array.
            names_.reserve(size);
            assert(names_.size() == 0);
            const char* currentChar = namesBuffer;
            for(std::size_t i = 0; i < size; ++i)
            {
                names_.emplace_back(currentChar);
                currentChar += names_.back().size()+1;
            }
        }
        delete[] namesBuffer;

        if( grid_.size(0) == 0 || names_.size() == 0) return;

        // Copy data from data arrays to map from local id to element data
        const auto& idSet = grid_.localIdSet();
        auto gridView = grid_.leafGridView();

        for( auto elem = gridView.template begin<0, Interior_Partition>(),
                 endElem = gridView.template end<0, Interior_Partition>();
             elem != endElem; ++elem)
        {
            const auto& id = idSet.id(*elem);
            const auto& index = factory_.insertionIndex(*elem);
            std::size_t dataIndex = 0;
            auto& dataOnId = elementData_[id];
            dataOnId.resize(elementDataSize_);

            for(const auto& array : dataArrays)
            {
                 dataOnId[dataIndex++] = array[index];
            }
        }
    }
    ~HDF5ReaderHandle()
    {
        distributedData_.clear();
        const auto& idSet = grid_.localIdSet();
        auto gridView = grid_.leafGridView();
        auto multiplier = multipliers_.begin();
        for(const auto name: names_)
        {
            distributedData_.emplace(name, std::vector<double>(gridView.size(0) *
                                                               *(multiplier++)));
        }
        auto mapper = MultipleCodimMultipleGeomTypeMapper< typename Grid::LeafGridView,
                                                       MCMGElementLayout >(gridView);

        // Remap the data from localId based to index based.
        for (const auto& element: elements(gridView))
        {
            const auto& toIndex = mapper.index(element);
            const auto& fromIndex = idSet.id(element);
            const auto& data = elementData_[fromIndex];
            std::size_t index = 0;

            for(const auto name: names_)
            {
                distributedData_[name][toIndex] = data[index++];
            }
        }
    }
    bool fixedsize(int /* dim */, int /* codim */) const
    {
        return true;
    }
    bool contains(int /* dim */, int  codim) const
    {
        return codim==0;
    }
    template<class T>
    std::size_t size(const T&) const
    {
        return elementDataSize_;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t ) const
    {
        auto iter = elementData_.find(grid_.localIdSet().id(t));
        assert(iter != elementData_.end());
        for(const auto& data : iter->second)
            buffer.write(data);
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t s)
    {
        assert(s==names_.size());
        auto& array = elementData_[grid_.localIdSet().id(t)];
        array.resize(s);
        for(auto& data: array)
        {
            buffer.read(data);
        }
    }
    private:
    const Grid& grid_;
    /// \brief The names associated with the vectors read from HDF5
    std::vector<std::string>& names_;
    typedef typename Grid::LocalIdSet LocalIdSet;
    /// \brief The data per element as a vector mapped from the local id.
    std::map<typename LocalIdSet::IdType, std::vector<double> > elementData_;
    /// \brief Upon exit the mapping of name to data array
    std::map<std::string, std::vector<double>>& distributedData_;
    const GridFactory<Grid>& factory_;
    /// \brief Whether  this is process 0.
    bool isRoot_;
    /// \brief The number of doubles per element
    std::size_t elementDataSize_;
    /// \brief The number of data items per element of a vector from HDF5
    std::vector<std::size_t> multipliers_;
};
} // end namespace Detail

/// \brief A rudimentary Reader reading a grid from a hdf5 file
///
/// This read will open the hdf5 file on rank 0 and upon request read the grid
/// information from it, create the global grid, read all the data of a given
/// time step, and loadbalance the grid with the data.
/// \code
///  auto reader = HDF5Reader("test.hdf5");
///  /* Read the first grid */
///  reader.readGrid();
///  /* Read all the cell data for the first timestep. */
///  reader.readAllCellData();
///  /* load balance the grid and data */
///  reader.loadBalance();
///  /* Retrieve the cell data */
///  const auto& reader.getCellData();
///  const auto& reader.grid();
/// \endcode
template<class Grid>
class HDF5Reader
{
public:
    using MPICommunicatorType = typename MPIHelper::MPICommunicator;
    /// \brief Create the HDF5 reader
    /// \param filename The filename including the extension.
    HDF5Reader(const std::string& filename, MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
        int rank = 0;
#if HAVE_MPI
      MPI_Comm_rank( comm, &rank );
#endif
      isRoot_ = ( rank == 0);

      if ( isRoot_ ){
            fileId_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        }
    }
    /// \brief Read the grid from the hdf file.
    void readGrid(const std::string& dataSetId = std::string("/grids/grid_0"));

    /// \brief Load balance the grid read on the master process and distribute date.
    void loadBalance();

    /// \brief Get the distributed data attached to the cells
    /// \warning A previous call to loadBalance is assumed.
    const std::map<std::string, std::vector<double> >& getCellData() const
    {
        return distributedCellData_;
    }

    /// \brief Get the distributed data attached to the cells
    /// \warning A previous call to loadBalance is assumed.
    /// \return A map that uses the associated name of the data as
    ///         the key to the data.
    std::map<std::string, std::vector<double> >& getCellData()
    {
        return distributedCellData_;
    }

    ~HDF5Reader()
    {
        if( isRoot_ )
            H5Fclose(fileId_);
    }
    /// \brief Get the created grid
    const Grid& grid() const
    {
        return *grid_;
    }
    /// \brief Get the created grid
    Grid& grid()
    {
        return *grid_;
    }

    /// \brief Read all data attached to the cells
    /// \param dataSetName The path to the data set of the time step in the hdf5 file.
    void readAllCellData(const std::string& dataSetName = std::string("/timesteps/timestep_0"));
private:
    bool isRoot_;
    GridFactory<Grid> factory_;
    std::unique_ptr<Grid> grid_;
    hid_t fileId_;
    std::vector<std::string> dataNames_;
    std::vector<std::vector<double> > cellData_;
    std::map<std::string, std::vector<double> > distributedCellData_;
    std::string readStringAttribute(hsize_t objectId, const char* name);
};

template<class Grid>
void HDF5Reader<Grid>::loadBalance()
{
    if ( grid_->comm().size()>1)
    {
        // Use loadbalance to distribute grid and data.
        // data handle for distributing cell data.
        Detail::HDF5ReaderHandle<Grid> handle(std::move(dataNames_),
                                              std::move(cellData_),
                                              distributedCellData_,
                                              factory_,
                                              *grid_);
        grid_->loadBalance(handle);
        // Cater for overlap/ghost cells
        auto gridView  = grid_->leafGridView();
        gridView.communicate(handle, InteriorBorder_All_Interface,ForwardCommunication);
    }
    else
    {
        // At this point the data arrays in cellData_ are ordered by the insertion index
        // but this must not by the same as the index of the element mapper
        // Therefore we resort.
        auto gridView  = grid_->leafGridView();
        auto dataArray = cellData_.begin();

        assert(distributedCellData_.size()==0);
        MultipleCodimMultipleGeomTypeMapper< typename Grid::LeafGridView,
                                             MCMGElementLayout > mapper(gridView);

        for( const auto& name : dataNames_ )
        {
            auto entry = distributedCellData_.emplace(name, std::vector<double>(dataArray->size()));
            assert(entry.second);
            auto& reorderedArray = entry.first->second;
            auto entriesPerIndex = dataArray->size()/grid_->size(0);

            for( const auto& elem : elements(gridView) )
            {
                const auto& insertionIndex = factory_.insertionIndex(elem);
                const auto& mappedIndex = mapper.index(elem);
                for(std::size_t i = 0; i < entriesPerIndex; ++i)
                    reorderedArray[mappedIndex*entriesPerIndex + i] =
                        (*dataArray)[insertionIndex*entriesPerIndex + i];
            }
            ++dataArray;
        }
        // Force clearing memory
        std::vector<std::string>().swap(dataNames_);
        std::vector<std::vector<double> >().swap(cellData_);
    }
}

template<class Grid>
void HDF5Reader<Grid>::readGrid(const std::string& dataSetId)
{
    if ( isRoot_ ){
        // read and insert vertices
        std::ostringstream dataSetName;
        dataSetName << dataSetId << "/vertices";

        auto myDataSet = H5Dopen(fileId_, dataSetName.str().c_str(), H5P_DEFAULT);
        auto memorySpace = H5Dget_space(myDataSet);
        const int ndims =  H5Sget_simple_extent_ndims(memorySpace);
        assert(ndims==2);
        hsize_t dims[ndims];
        H5Sget_simple_extent_dims(memorySpace, dims, NULL);
        if(dims[1] != Grid::dimensionworld)
            DUNE_THROW(RangeError, "dimensionworld in HDF5 file is "<<dims[1]<<" expected "<<Grid::dimensionworld);

        std::vector<double> myData(dims[0]*dims[1]);

        auto status = H5Dread(myDataSet, H5T_NATIVE_DOUBLE, H5S_ALL,
                              H5S_ALL, H5P_DEFAULT, myData.data());
        status = H5Dclose(myDataSet);
        for(auto vertex = myData.begin(), vertexEnd = myData.end();
            vertex != vertexEnd;)
        {
            auto dataEnd = vertex + Grid::dimensionworld;
            Dune::FieldVector<typename Grid::ctype, Grid::dimensionworld> data;
            std::copy(vertex, dataEnd, data.begin());
            factory_.insertVertex(data);
            vertex = dataEnd;
        }

        // read and insert Elements
        // Currently only triangles supported!
        dataSetName.clear();
        dataSetName.str("");
        dataSetName << dataSetId << "/topology";
        std::cout<<" Reading "<<dataSetName.str().c_str()<<std::endl;
        myDataSet = H5Dopen(fileId_, dataSetName.str().c_str(), H5P_DEFAULT);
        memorySpace = H5Dget_space(myDataSet);
        const int ndims1 =  H5Sget_simple_extent_ndims(memorySpace);
        hsize_t dims1[ndims];
        H5Sget_simple_extent_dims(memorySpace, dims1, NULL);
        if( ndims1 != 1 )
             DUNE_THROW(RangeError, "Expected one dimensional array for topology");
        if( dims1[0]%3 != 0 )
            DUNE_THROW(RangeError, "Expected triangle topology but size of "
                       << "data space mismatches.");
        myData.resize(dims1[0]);
        status = H5Dread(myDataSet, H5T_NATIVE_DOUBLE, H5S_ALL,
                              H5S_ALL, H5P_DEFAULT, myData.data());
        status = H5Dclose(myDataSet);

        for(auto element = myData.begin(), elementEnd = myData.end();
            element != elementEnd;)
        {
            auto dataEnd = element + 3;
            std::vector<unsigned int> vertices(element, dataEnd);
            factory_.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension),
                                   vertices);
            element = dataEnd;
        }
        grid_.reset(factory_.createGrid());
    }
    else
    {
        grid_.reset(factory_.createGrid());
    }
}

template<class Grid>
std::string HDF5Reader<Grid>::readStringAttribute(hsize_t objectId,
                                                  const char* name)
{

    auto attributeId = H5Aopen(objectId, name, H5P_DEFAULT);
    if(attributeId < 0)
        return std::string();
    auto fileType = H5Aget_type(attributeId);
    auto space = H5Aget_space(attributeId);
    // Create memory type
    auto stringDataType = H5Tcopy (H5T_C_S1);
    H5Tset_size (stringDataType, H5T_VARIABLE);
    // Read Data
    //hsize_t dims[1] = { 1 };
    const char* data[1];
    auto status = H5Aread(attributeId, stringDataType, data);
    std::string result(data[0]);
    status = H5Dvlen_reclaim (stringDataType, space, H5P_DEFAULT, data);
    status = H5Aclose(attributeId);
    return result;
}

template<class Grid>
void HDF5Reader<Grid>::readAllCellData(const std::string& dataSetName)
{
    dataNames_.clear();
    cellData_.clear();
    if ( isRoot_ )
    {
        auto groupId = H5Gopen(fileId_, dataSetName.c_str(), H5P_DEFAULT);
        hsize_t numObjects;
        auto status = H5Gget_num_objs(groupId, &numObjects);
        for(hsize_t i = 0; i < numObjects; ++i)
        {
            char name[256];
            auto nameLength = H5Gget_objname_by_idx(groupId, i, name, 256);
            if ( nameLength <= 0)
                DUNE_THROW(RangeError, "Could not determine name for timestep data");
            std::cout << "data name is "<<name<<" length="<<nameLength;
            auto myDataSet = H5Dopen(groupId, name, H5P_DEFAULT);
            // Check for attribute name that overwrites the name of the data array
            if( H5Aexists(myDataSet, "name") > 0)
            {
                auto nameFromAttribute = readStringAttribute(myDataSet, "name");
                dataNames_.push_back(nameFromAttribute);
            }
            else
            {
                dataNames_.push_back(name);
            }
            auto memorySpace = H5Dget_space(myDataSet);
            const int ndims =  H5Sget_simple_extent_ndims(memorySpace);
            hsize_t dims[ndims];
            H5Sget_simple_extent_dims(memorySpace, dims, NULL);
            std::size_t size = dims[0];
            for(int i=1; i < ndims; ++i)
                size *= dims[i];
            cellData_.emplace_back(size);
            auto status = H5Dread(myDataSet, H5T_NATIVE_DOUBLE, H5S_ALL,
                                  H5S_ALL, H5P_DEFAULT, cellData_.back().data());
            H5Dclose(myDataSet);
        }
        H5Gclose(groupId);
    }
}
} // end namespace Dumux
} // end namespace Dune
#endif // DUMUX_HDF5GRIDREADER_HH
