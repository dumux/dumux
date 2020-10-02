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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup FacetCoupling
 * \brief Contains the grid manager class that creates the grids in the context
 *        of hybrid-dimensional coupled models, where the (n-1)-dimensional
 *        domains live on the element facets of the n-dimensional domains.
 *        Also, it allows to extract a grid data object containing parameters
 *        passed to elements and/or boundary segments. All grids are constructed
 *        from a single grid file.
 */
#ifndef DUMUX_FACETCOUPLING_GRID_MANAGER_HH
#define DUMUX_FACETCOUPLING_GRID_MANAGER_HH

#include <cassert>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/io/grid/griddata.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/typetraits/utility.hh>

#include "gmshreader.hh"

namespace Dumux {
namespace FCGridManagerChecks {

    // The grid creator and grid data classes provided below
    // require that the grids passed as template arguments are
    // ordered in descending grid dimension and they must have
    // the same world dimension. The following helper structs
    // verify these conditions and are reused below.
    // evaluates if dim is descending or dimworld is equal for two grids
    template<bool checkDimWorld, typename G1, typename G2>
    static constexpr bool evalCondition()
    {
        return checkDimWorld ? int(G1::dimensionworld) == int(G2::dimensionworld)
                             : int(G1::dimension) > int(G2::dimension);
    }

    // helper structs to evaluate conditions on the grid
    template<bool checkDimWorld, typename... Gs> struct FulfillConditions;
    template<bool checkDimWorld, typename G1, typename... Gs>
    struct FulfillConditions<checkDimWorld, G1, Gs...>
    {
        using G2 = typename std::tuple_element_t<0, std::tuple<Gs...>>;
        static constexpr bool value = evalCondition<checkDimWorld, G1, G2>() && FulfillConditions<checkDimWorld, Gs...>::value;
    };
    template<bool checkDimWorld, typename G1, typename G2>
    struct FulfillConditions<checkDimWorld, G1, G2> { static constexpr bool value = evalCondition<checkDimWorld, G1, G2>(); };
}

/*!
 * \ingroup FacetCoupling
 * \brief Grid data object to store element and boundary segment markers
 *        for all grids of the hierarchy.
 *
 * \tparam Grids the types of the grid hierarchy
 * \note Grids must be ordered in descending grid dimension
 */
template<typename... Grids>
class FacetCouplingGridDataWrapper
{
    // make sure all grids have the same world dimension and are ordered in descending dimension
    static_assert(FCGridManagerChecks::FulfillConditions<false, Grids...>::value, "All grids must have the same world dimension!");
    static_assert(FCGridManagerChecks::FulfillConditions<true, Grids...>::value, "Grids must be ordered w.r.t the dimension in descending order!");

    //! determine the number of involved grids
    static constexpr std::size_t numGrids = sizeof...(Grids);
    //! the i-th grid type
    template<std::size_t id> using Grid = typename std::tuple_element_t<id, std::tuple<Grids...>>;
    //! shared ptr to the i-th grid data type
    template<std::size_t id> using GridDataPtr = std::shared_ptr< GridData<Grid<id>> >;
    //! the intersection type of the i-th grid
    template<std::size_t id> using Intersection = typename Grid<id>::LeafGridView::Intersection;

public:
    //! export the i-th grid data type
    template<std::size_t id> using GridData = GridData<Grid<id>>;

    //! set the grid data object for the i-th grid
    template<std::size_t id>
    void setGridData(GridData<id>&& gridData)
    {
        std::get<id>(gridDataPtrTuple_) = std::make_shared<GridData<id>>( std::move(gridData) );
    }

    //! return the grid data for a specific grid
    template<std::size_t id>
    std::shared_ptr<const GridData<id>> getSubDomainGridData() const
    { return std::get<id>(gridDataPtrTuple_); }

    //! Returns domain marker of an element
    template<std::size_t id>
    int getElementDomainMarker(const typename Grid<id>::template Codim<0>::Entity& element) const
    { return std::get<id>(gridDataPtrTuple_)->getElementDomainMarker(element); }

    //! Returns the boundary marker of an intersection
    template<std::size_t id>
    int getBoundaryDomainMarker(const typename Grid<id>::LeafGridView::Intersection& is) const
    { return std::get<id>(gridDataPtrTuple_)->getBoundaryDomainMarker(is); }

    //! Returns the boundary marker for a given bounday segment index
    template<std::size_t id>
    int getBoundaryDomainMarker(int boundarySegmentIndex) const
    { return std::get<id>(gridDataPtrTuple_)->getBoundaryDomainMarker(boundarySegmentIndex); }

    //! Returns true if an intersection was inserted during grid creation
    template<std::size_t id>
    bool wasInserted(const Intersection<id>& intersection) const
    { return std::get<id>(gridDataPtrTuple_)->wasInserted(intersection); }

private:
    //! store a shared pointer to a grid data object for each grid
    using GridDataPtrTuple = typename makeFromIndexedType<std::tuple, GridDataPtr, std::make_index_sequence<numGrids>>::type;
    GridDataPtrTuple gridDataPtrTuple_;
};

/*!
 * \ingroup FacetCoupling
 * \brief Contains the embeddings between grids with codimension one
 *        among the grid hierarchy. All these embedments are given in
 *        insertion indices as they are read directly from the grid.
 *        Therefore, this class furthermore allows access to the insertion
 *        indices of entities. Additionally, it gives access to the grid
 *        views of the different grids on the hierarchy.
 *
 * \tparam Grids the types of the grid hierarchy
 * \note Grids must be ordered in descending grid dimension
 */
template<typename... Grids>
class FacetCouplingEmbeddings
{
    // make sure all grids have the same world dimension and are ordered in descending dimension
    static_assert(FCGridManagerChecks::FulfillConditions<false, Grids...>::value, "All grids must have the same world dimension!");
    static_assert(FCGridManagerChecks::FulfillConditions<true, Grids...>::value, "Grids must be ordered w.r.t the dimension in descending order!");

    //! the i-th grid type
    template<std::size_t id> using Grid = typename std::tuple_element_t<id, std::tuple<Grids...>>;
    //! the i-th grid factory type
    template<std::size_t id> using GridFactory = typename Dune::GridFactory<Grid<id>>;

    //! we use the bulk grid's index type here
    using GIType = typename IndexTraits< typename Grid<0>::LeafGridView >::GridIndex;
    //! the map type to store embedment data
    using EmbedmentMap = std::unordered_map<GIType, std::vector<GIType>>;

public:
    //! export the i-th grid view type
    template<std::size_t id> using GridView = typename Grid<id>::LeafGridView;

    //! export the number of created grids
    static constexpr std::size_t numGrids = sizeof...(Grids);
    //! export the grid id of the bulk grid (descending grid dim -> always zero!)
    static constexpr int bulkGridId = 0;
    //! state the dimension of the highest-dimensional grid
    static constexpr int bulkDim = Grid<bulkGridId>::dimension;

    //! export the bulk grid type
    using BulkGridView = GridView<bulkGridId>;
    //! export the type used for indices
    using GridIndexType = GIType;

    //! return reference to the i-th grid view
    template<std::size_t id>
    const GridView<id>& gridView() const
    { return *std::get<id>(gridViewPtrTuple_); }

    //! return the insertion index of an entity of the i-th grid
    template<std::size_t id, class Entity>
    GridIndexType insertionIndex(const Entity& entity) const
    { return std::get<id>(gridFactoryPtrTuple_)->insertionIndex(entity); }

    //! Returns the insertion indices of the entities embedded in given element
    template<std::size_t id>
    typename std::unordered_map< GridIndexType, std::vector<GridIndexType> >::mapped_type
    embeddedEntityIndices(const typename Grid<id>::template Codim<0>::Entity& element) const
    {
        const auto& map = embeddedEntityMaps_[id];
        auto it = map.find( std::get<id>(gridFactoryPtrTuple_)->insertionIndex(element) );
        if (it != map.end()) return it->second;
        else return typename std::unordered_map< GridIndexType, std::vector<GridIndexType> >::mapped_type();
    }

    //! Returns the insertion indices of the entities in which the element is embedded
    template<std::size_t id>
    typename std::unordered_map< GridIndexType, std::vector<GridIndexType> >::mapped_type
    adjoinedEntityIndices(const typename Grid<id>::template Codim<0>::Entity& element) const
    {
        const auto& map = adjoinedEntityMaps_[id];
        auto it = map.find( std::get<id>(gridFactoryPtrTuple_)->insertionIndex(element) );
        if (it != map.end()) return it->second;
        else return typename std::unordered_map< GridIndexType, std::vector<GridIndexType> >::mapped_type();
    }

    //! Returns const reference to maps of the embedded entities
    const std::unordered_map< GridIndexType, std::vector<GridIndexType> >& embeddedEntityMap(std::size_t id) const
    { assert(id < numGrids); return embeddedEntityMaps_[id]; }

    //! Returns non-const reference to maps of the embedded entities
    std::unordered_map< GridIndexType, std::vector<GridIndexType> >& embeddedEntityMap(std::size_t id)
    { assert(id < numGrids); return embeddedEntityMaps_[id]; }

    //! Returns const reference to the maps of the adjoined entities of dimension d+1
    const std::unordered_map< GridIndexType, std::vector<GridIndexType> >& adjoinedEntityMap(std::size_t id) const
    { assert(id < numGrids); return adjoinedEntityMaps_[id]; }

    //! Returns non-const reference to the maps of the adjoined entities of dimension d+1
    std::unordered_map< GridIndexType, std::vector<GridIndexType> >& adjoinedEntityMap(std::size_t id)
    { assert(id < numGrids); return adjoinedEntityMaps_[id]; }

    //! Returns the hierachy's insertion indices that make up the grid for the given id
    const std::vector<GridIndexType>& gridHierarchyIndices(std::size_t id) const
    { assert(id < numGrids); return gridVertexIndices_[id]; }

    //! Returns the number of vertices contained in the entire grid hierarch
    std::size_t numVerticesInHierarchy() const
    { return numVerticesInHierarchy_; }

    /*!
     * \brief Sets the required data for a specific grid on the hierarchy.
     * \param gridPtr shared pointer to this (id-th) grid
     * \param gridFactoryPtr shared pointer to this (id-th) grid factory
     * \param embeddedEntityMap map containing the lower-dimensional entities for this grid's elements
     * \param adjoinedEntityMap map containing the (d+1)-dimensional elements in which this grid's elements are embedded in
     * \param gridVertexIndices The hierachy's insertion indices that make up this grid
     * \param numVerticesInHierarchy Total number of vertices in entire grid hierarchy
     */
    template<std::size_t id>
    void setData( std::shared_ptr<Grid<id>> gridPtr,
                  std::shared_ptr<GridFactory<id>> gridFactoryPtr,
                  EmbedmentMap&& embeddedEntityMap,
                  EmbedmentMap&& adjoinedEntityMap,
                  std::vector<GridIndexType>&& gridVertexIndices,
                  std::size_t numVerticesInHierarchy )
    {
        std::get<id>(gridViewPtrTuple_) = std::make_shared<GridView<id>>(gridPtr->leafGridView());
        std::get<id>(gridFactoryPtrTuple_) = gridFactoryPtr;
        embeddedEntityMaps_[id] = std::move(embeddedEntityMap);
        adjoinedEntityMaps_[id] = std::move(adjoinedEntityMap);
        gridVertexIndices_[id] = std::move(gridVertexIndices);
        numVerticesInHierarchy_ = numVerticesInHierarchy;
    }

private:
    //! data on connectivity between the grids
    std::array<EmbedmentMap, numGrids> embeddedEntityMaps_;
    std::array<EmbedmentMap, numGrids> adjoinedEntityMaps_;

    //! Contains the hierarchy insertion indices that make up a lower-dimensional grid
    std::size_t numVerticesInHierarchy_;
    std::array<std::vector<GridIndexType>, numGrids> gridVertexIndices_;

    //! tuple to store the grids
    using Indices = std::make_index_sequence<numGrids>;
    template<std::size_t id> using GridViewPtr = std::shared_ptr<GridView<id>>;
    using GridPtrTuple = typename makeFromIndexedType<std::tuple, GridViewPtr, Indices>::type;
    GridPtrTuple gridViewPtrTuple_;

    //! tuple to store the grid grid factories
    template<std::size_t id> using GridFactoryPtr = std::shared_ptr< Dune::GridFactory<Grid<id>> >;
    using GridFactoryPtrTuple = typename makeFromIndexedType<std::tuple, GridFactoryPtr, Indices>::type;
    GridFactoryPtrTuple gridFactoryPtrTuple_;
};

/*!
 * \ingroup FacetCoupling
 * \brief Creates the grids in the context of hybrid-dimensional coupled models,
 *        where the (n-1)-dimensional domains live on the element facets of the
 *        n-dimensional domains.
 *
 * \tparam Grids the types of the grid hierarchy
 * \note Grids must be ordered in descending grid dimension
 */
template<typename... Grids>
class FacetCouplingGridManager
{
    // make sure all grids have the same world dimension and are ordered in descending dimension
    static_assert(FCGridManagerChecks::FulfillConditions<false, Grids...>::value, "All grids must have the same world dimension!");
    static_assert(FCGridManagerChecks::FulfillConditions<true, Grids...>::value, "Grids must be ordered w.r.t the dimension in descending order!");

    // we use a wrapper class for the grid data containing the data on all grids
    using GridDataWrapper = FacetCouplingGridDataWrapper<Grids...>;
public:
    //! export the i-th grid type
    template<std::size_t id> using Grid = typename std::tuple_element_t<id, std::tuple<Grids...>>;
    //! export the i-th grid pointer type
    template<std::size_t id> using GridPtr = typename std::shared_ptr< Grid<id> >;

    //! export the number of created grids
    static constexpr std::size_t numGrids = sizeof...(Grids);
    //! export the grid id of the bulk grid (descending grid dim -> always zero!)
    static constexpr int bulkGridId = 0;

    //! export the grid data (wrapper) type, i.e. parameters/markers
    using GridData = GridDataWrapper;
    //! export the type storing the embeddings
    using Embeddings = FacetCouplingEmbeddings<Grids...>;

    //! returns the i-th grid
    template<std::size_t id>
    const Grid<id>& grid() const
    { return *std::get<id>(gridPtrTuple_); }

    //! return a pointer to the grid data object
    std::shared_ptr<const GridData> getGridData() const
    {
        if (!enableEntityMarkers_)
            DUNE_THROW(Dune::IOError, "No grid data available");
        return gridDataPtr_;
    }

    //! return a pointer to the object containing embeddings
    std::shared_ptr<const Embeddings> getEmbeddings() const
    { return embeddingsPtr_; }

    //! creates the grids from a file given in parameter tree
    void init(const std::string& paramGroup = "")
    {
        // reset the grid & embedding data
        gridDataPtr_ = std::make_shared<GridData>();
        embeddingsPtr_ = std::make_shared<Embeddings>();

        // get filename and determine grid file extension
        const auto fileName = getParamFromGroup<std::string>(paramGroup, "Grid.File");
        const auto ext = getFileExtension(fileName);

        // get some parameters
        const bool verbose = getParamFromGroup<bool>(paramGroup, "Grid.Verbosity", false);
        const bool domainMarkers = getParamFromGroup<bool>(paramGroup, "Grid.DomainMarkers", false);
        const bool boundarySegments = getParamFromGroup<bool>(paramGroup, "Grid.BoundarySegments", false);

        // forward to the corresponding reader
        if (ext == "msh")
        {
            const auto thresh = getParamFromGroup<std::size_t>(paramGroup, "Grid.GmshPhysicalEntityThreshold", 0);
            FacetCouplingGmshReader<Grid<bulkGridId>, numGrids> gmshReader;
            gmshReader.read(fileName, (boundarySegments ? thresh : 0), verbose);
            passDataFromReader(gmshReader, domainMarkers, boundarySegments);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Reader for grid files of type ." + ext);

        // find out if entity markers are active
        enableEntityMarkers_ = domainMarkers || boundarySegments;
    }

    //! Distributes the grid on all processes of a parallel computation
    void loadBalance()
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridPtrTuple_)), [&](const auto id)
        {
            std::get<id>(this->gridPtrTuple_)->loadBalance();
        });
    }

protected:
    //! return non-const reference to i-th grid
    template<std::size_t id>
    Grid<id>& grid_()
    { return *std::get<id>(gridPtrTuple_); }

    //! return non-const pointer to the object containing embeddings
    std::shared_ptr<Embeddings> getEmbeddings_()
    { return embeddingsPtr_; }

private:
    //! Returns the filename extension of a given filename
    static std::string getFileExtension(const std::string& fileName)
    {
        const auto pos = fileName.rfind('.', fileName.length());
        if (pos != std::string::npos)
            return(fileName.substr(pos+1, fileName.length() - pos));
        else
            DUNE_THROW(Dune::IOError, "Please provide an extension for your grid file ('"<< fileName << "')!");
    }

    //! Creates the grids using the data in a mesh file reader
    template<typename MeshFileReader>
    void passDataFromReader(MeshFileReader& reader, bool domainMarkers, bool boundarySegments)
    {
        const auto& vertices = reader.gridVertices();

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridPtrTuple_)), [&](const auto id)
        {
            using GridFactory = Dune::GridFactory<Grid<id>>;
            auto factoryPtr = std::make_shared<GridFactory>();

            // insert grid vertices
            for (const auto idx : reader.vertexIndices(id))
                factoryPtr->insertVertex(vertices[idx]);

            // insert elements
            for (const auto& e : reader.elementData(id))
                factoryPtr->insertElement(e.gt, e.cornerIndices);

            // insert boundary segments
            if (boundarySegments)
                for (const auto& segment : reader.boundarySegmentData(id))
                    factoryPtr->insertBoundarySegment(segment);

            // make grid
            auto gridPtr = std::shared_ptr<Grid<id>>(factoryPtr->createGrid());

            // maybe create and set grid data object
            if (domainMarkers || boundarySegments)
            {
                typename GridDataWrapper::template GridData<id> gridData( gridPtr,
                                                                          factoryPtr,
                                                                          std::move(reader.elementMarkerMap(id)),
                                                                          std::move(reader.boundaryMarkerMap(id)) );
                gridDataPtr_->template setGridData<id>( std::move(gridData) );
            }

            // copy the embeddings
            embeddingsPtr_->template setData<id>( gridPtr,
                                                  factoryPtr,
                                                  std::move(reader.embeddedEntityMap(id)),
                                                  std::move(reader.adjoinedEntityMap(id)),
                                                  std::move(reader.vertexIndices(id)),
                                                  vertices.size() );

            // set the grid pointer
            std::get<id>(gridPtrTuple_) = gridPtr;
        });
    }

    //! tuple to store the grids
    using Indices = std::make_index_sequence<numGrids>;
    using GridPtrTuple = typename makeFromIndexedType<std::tuple, GridPtr, Indices>::type;
    GridPtrTuple gridPtrTuple_;

    //! grid data, i.e. parameters and markers
    bool enableEntityMarkers_;
    std::shared_ptr<GridData> gridDataPtr_;

    //! data on embeddings
    std::shared_ptr<Embeddings> embeddingsPtr_;
};

} // end namespace Dumux

#endif
