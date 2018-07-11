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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
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
#include <dumux/common/typetraits/utility.hh>

#include "gmshreader.hh"

namespace Dumux {
namespace Detail {

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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
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
    static_assert(Detail::FulfillConditions<false, Grids...>::value, "All grids must have the same world dimension!");
    static_assert(Detail::FulfillConditions<true, Grids...>::value, "Grids must be ordered w.r.t the dimension in descending order!");

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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
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
    static_assert(Detail::FulfillConditions<false, Grids...>::value, "All grids must have the same world dimension!");
    static_assert(Detail::FulfillConditions<true, Grids...>::value, "Grids must be ordered w.r.t the dimension in descending order!");

    // we use a wrapper class for the grid data containing the data on all grids
    using GridDataWrapper = FacetCouplingGridDataWrapper<Grids...>;
public:
    //! export the i-th grid type
    template<std::size_t id> using Grid = typename std::tuple_element_t<id, std::tuple<Grids...>>;
    //! export the i-th grid pointer type
    template<std::size_t id> using GridPtr = typename std::shared_ptr< Grid<id> >;
    //! export the i-th Dune::GridFactory type
    template<std::size_t id> using GridFactory = Dune::GridFactory< Grid<id> >;
    //! export the i-th Dune::GridFactory type
    template<std::size_t id> using GridFactoryPtr = std::shared_ptr< Dune::GridFactory<Grid<id>> >;

    //! export the number of created grids
    static constexpr std::size_t numGrids = sizeof...(Grids);
    //! export the grid id of the bulk grid (descending grid dim -> always zero!)
    static constexpr int bulkGridId = 0;
    //! state the dimension of the highest-dimensional grid
    static constexpr int bulkDim = Grid<bulkGridId>::dimension;

    //! export the bulk grid type
    using BulkGrid = Grid<bulkGridId>;
    //! we use the bulk grid's index type during read etc.
    using IndexType = typename BulkGrid::LeafGridView::IndexSet::IndexType;
    //! export the grid data (wrapper) type, i.e. parameters/markers
    using GridData = GridDataWrapper;

    //! returns the i-th grid
    template<std::size_t id>
    const Grid<id>& grid() const
    { return *std::get<id>(gridPtrTuple_); }

    //! returns the i-th grid factory
    template<std::size_t id>
    const GridFactory<id>& gridFactory() const
    { return *std::get<id>(gridFactoryPtrTuple_); }

    //! return a pointer to the grid data object
    std::shared_ptr<const GridData> getGridData() const
    {
        if (!enableEntityMarkers_)
            DUNE_THROW(Dune::IOError, "No grid data available");
        return gridDataPtr_;
    }

    //! Returns the insertion indices of the entities embedded in given element
    template<std::size_t id>
    typename std::unordered_map< IndexType, std::vector<IndexType> >::mapped_type
    embeddedEntityIndices(const typename Grid<id>::template Codim<0>::Entity& element) const
    {
        const auto& map = embeddedEntityMaps_[id];
        auto it = map.find( gridFactory<id>().insertionIndex(element) );
        if (it != map.end()) return it->second;
        else return typename std::unordered_map< IndexType, std::vector<IndexType> >::mapped_type();
    }

    //! Returns the insertion indices of the entities in which the element is embedded
    template<std::size_t id>
    typename std::unordered_map< IndexType, std::vector<IndexType> >::mapped_type
    embedmentEntityIndices(const typename Grid<id>::template Codim<0>::Entity& element) const
    {
        const auto& map = embedmentMaps_[id];
        auto it = map.find( gridFactory<id>().insertionIndex(element) );
        if (it != map.end()) return it->second;
        else return typename std::unordered_map< IndexType, std::vector<IndexType> >::mapped_type();
    }

    //! Returns the maps of the embedded entities
    const std::unordered_map< IndexType, std::vector<IndexType> >& embeddedEntityMap(std::size_t id) const
    { assert(id < numGrids); return embeddedEntityMaps_[id]; }

    std::unordered_map< IndexType, std::vector<IndexType> >& embeddedEntityMap(std::size_t id)
    { assert(id < numGrids); return embeddedEntityMaps_[id]; }

    //! Returns the maps of the embedments
    const std::unordered_map< IndexType, std::vector<IndexType> >& embedmentMap(std::size_t id) const
    { assert(id < numGrids); return embedmentMaps_[id]; }

    std::unordered_map< IndexType, std::vector<IndexType> >& embedmentMap(std::size_t id)
    { assert(id < numGrids); return embedmentMaps_[id]; }

    //! Returns the hierachy's insertion indices that make up the grid for the given id
    const std::vector<IndexType>& lowDimVertexIndices(std::size_t id) const
    { assert(id > 0 && id < numGrids); return lowDimGridVertexIndices_[id-1]; }

    //! creates the grids from a file given in parameter tree
    void init(const std::string& paramGroup = "")
    {
        // reset the grid data
        gridDataPtr_ = std::make_shared<GridData>();

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
            FacetCouplingGmshReader<BulkGrid, IndexType, numGrids> gmshReader;
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
        const auto& bulkGridVertices = reader.bulkGridVertices();
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridPtrTuple_)), [&](const auto id)
        {
            auto factoryPtr = std::make_shared<GridFactory<id>>();
            std::get<id>(this->gridFactoryPtrTuple_) = factoryPtr;

            // insert grid vertices
            if (id == 0)
                for (const auto& v : bulkGridVertices)
                    factoryPtr->insertVertex(v);
            else
            {
                for (const auto idx : reader.lowDimVertexIndices(id))
                    factoryPtr->insertVertex(bulkGridVertices[idx]);
                std::swap(lowDimGridVertexIndices_[id-1], reader.lowDimVertexIndices(id));
            }

            // insert elements
            for (const auto& e : reader.elementData(id))
                factoryPtr->insertElement(e.gt, e.cornerIndices);

            // insert boundary segments
            if (boundarySegments)
                for (const auto& segment : reader.boundarySegmentData(id))
                    factoryPtr->insertBoundarySegment(segment);

            // make grid
            std::get<id>(gridPtrTuple_) = std::shared_ptr<Grid<id>>(factoryPtr->createGrid());

            // maybe create and set grid data object
            if (domainMarkers || boundarySegments)
            {
                typename GridDataWrapper::template GridData<id> gridData( std::get<id>(gridPtrTuple_),
                                                                          factoryPtr,
                                                                          std::move(reader.elementMarkerMap(id)),
                                                                          std::move(reader.boundaryMarkerMap(id)) );
                gridDataPtr_->template setGridData<id>( std::move(gridData) );
            }

            // copy the embedments
            std::swap(embeddedEntityMaps_[id], reader.embeddedEntityMap(id));
            std::swap(embedmentMaps_[id], reader.embedmentMap(id));
        });
    }

    using Indices = std::make_index_sequence<numGrids>;
    //! tuple to store the grids
    using GridPtrTuple = typename makeFromIndexedType<std::tuple, GridPtr, Indices>::type;
    GridPtrTuple gridPtrTuple_;

    //! tuple to store the grid grid factories
    using GridFactoryPtrTuple = typename makeFromIndexedType<std::tuple, GridFactoryPtr, Indices>::type;
    GridFactoryPtrTuple gridFactoryPtrTuple_;

    //! data on connectivity between the grids
    std::array< std::unordered_map< IndexType, std::vector<IndexType> >, numGrids > embeddedEntityMaps_;
    std::array< std::unordered_map< IndexType, std::vector<IndexType> >, numGrids > embedmentMaps_;

    //! maps a low-dim vertex insertion index to the bulk grid vertex insertion index
    std::array< std::vector<IndexType>, numGrids-1 > lowDimGridVertexIndices_;

    //! grid data, i.e. parameters and markers
    bool enableEntityMarkers_;
    std::shared_ptr<GridData> gridDataPtr_;
};

} // end namespace Dumux

#endif // DUMUX_FACETCOUPLING_GRID_CREATOR_HH
