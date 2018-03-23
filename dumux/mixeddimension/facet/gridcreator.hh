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
 * \ingroup FacetCoupling
 * \brief Contains the grid creator class that creates the grids in the context
 *        of hybrid-dimensional coupled models, where the (n-1)-dimensional
 *        domains live on the element facets of the n-dimensional domains.
 *        All grids are constructed from a single grid file.
 */
#ifndef DUMUX_FACETCOUPLING_GRID_CREATOR_HH
#define DUMUX_FACETCOUPLING_GRID_CREATOR_HH

#include <cassert>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/utility.hh>

#include "gmshreader.hh"

namespace Dumux
{

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
class FacetCouplingGridCreator
{
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

    // make sure all grids have the same world dimension and are ordered in descending dimension
    static_assert(FulfillConditions<false, Grids...>::value, "All grids must have the same world dimension!");
    static_assert(FulfillConditions<true, Grids...>::value, "Grids must be ordered w.r.t the dimension in descending order!");

public:
    //! export the i-th grid type
    template<std::size_t id> using Grid = typename std::tuple_element_t<id, std::tuple<Grids...>>;
    //! export the i-th grid pointer type
    template<std::size_t id> using GridPtr = typename std::shared_ptr<Grid<id>>;
    //! export the i-th Dune::GridFactory type
    template<std::size_t id> using GridFactory = Dune::GridFactory<Grid<id>>;

    //! returns the i-th grid
    template<std::size_t id> Grid<id>& grid() { return *std::get<id>(gridPtrTuple_); }
    template<std::size_t id> const Grid<id>& grid() const { return *std::get<id>(gridPtrTuple_); }
    //! returns the i-th grid factory
    template<std::size_t id> GridFactory<id>& gridFactory() { return std::get<id>(gridFactoryTuple_); }
    template<std::size_t id> const GridFactory<id>& gridFactory() const { return std::get<id>(gridFactoryTuple_); }

    //! Encapsulates types used by this class
    struct Traits
    {
        //! export the number of created grids
        static constexpr std::size_t numGrids = sizeof...(Grids);
        //! export the bulk grid type
        using BulkGrid = Grid<0>;
        //! we use the bulk grid's index type during read etc.
        using IndexType = typename BulkGrid::LeafGridView::IndexSet::IndexType;
        //! maps to n-dimensional elements the set of (n-1)-dimensional elements embedded in it
        using EmbeddedEntityMap = std::map< IndexType, std::vector<IndexType> >;
        //! maps to m-dimensional elements the set of (n+1)-dimensional elements in which they are embedded
        using EmbedmentMap = std::map< IndexType, std::vector<IndexType> >;
        //! maps to each element a domain marker
        using ElementToDomainMarkerMap = std::vector<int>;
        //! maps to each boundary segment a boundary marker
        using BoundarySegmentToMarkerMap = std::vector<int>;
    };

    //! Returns domain marker of an element
    template<std::size_t id>
    typename Traits::ElementToDomainMarkerMap::value_type
    getElementDomainMarker(const typename Grid<id>::template Codim<0>::Entity& element)
    { return elementMarkerMaps_[id][gridFactory<id>().insertionIndex(element)]; }

    //! Returns the boundary marker of an intersection
    template<std::size_t id>
    typename Traits::BoundarySegmentToMarkerMap::value_type
    getBoundaryDomainMarker(const typename Grid<id>::LeafGridView::Intersection& is)
    {
        // this should only be called for intersections that were inserted
        assert(gridFactory<id>().wasInserted(is) && "Can't obtain boundary markers for intersections that weren't inserted!");
        return boundaryMarkerMaps_[id][gridFactory<id>().insertionIndex(is)];
    }

    //! Returns the insertion indices of the entities embedded in given element
    template<std::size_t id>
    const typename Traits::EmbeddedEntityMap::mapped_type
    embeddedEntityIndices(const typename Grid<id>::template Codim<0>::Entity& element)
    {
        const auto& map = embeddedEntityMaps_[id];
        auto it = map.find( gridFactory<id>().insertionIndex(element) );
        if (it != map.end()) return it->second;
        else return typename Traits::EmbeddedEntityMap::mapped_type();
    }

    //! Returns the insertion indices of the entities in which the element is embedded
    template<std::size_t id>
    const typename Traits::EmbedmentMap::mapped_type
    embedmentEntityIndices(const typename Grid<id>::template Codim<0>::Entity& element)
    {
        const auto& map = embedmentMaps_[id];
        auto it = map.find( gridFactory<id>().insertionIndex(element) );
        if (it != map.end()) return it->second;
        else return typename Traits::EmbedmentMap::mapped_type();
    }

    //! Returns the maps of the embedded entities
    typename Traits::EmbeddedEntityMap& embeddedEntityMap(std::size_t id)
    { assert(id < numGrids); return embeddedEntityMaps_[id]; }

    //! Returns the maps of the embedments
    typename Traits::EmbedmentMap& embedmentMap(std::size_t id)
    { assert(id < numGrids); return embedmentMaps_[id]; }

    //! creates the grids from a given grid file using the parameter tree
    void makeGrids(const std::string& fileName, const std::string& paramGroup = "")
    {
        const auto ext = getFileExtension(fileName);

        //! get some parameters
        const bool verbose = getParamFromGroup<bool>(paramGroup, "Grid.Verbosity", false);
        const bool domainMarkers = getParamFromGroup<bool>(paramGroup, "Grid.DomainMarkers", false);
        const bool boundarySegments = getParamFromGroup<bool>(paramGroup, "Grid.BoundarySegments", false);

        //! forward to the corresponding reader
        if (ext == "msh")
        {
            const auto thresh = getParamFromGroup<std::size_t>(paramGroup, "Grid.GmshPhysicalEntityThreshold", 0);
            FacetCouplingGmshReader<Traits> gmshReader;
            gmshReader.read(fileName, (boundarySegments ? thresh : 0), verbose);
            passDataFromReader(gmshReader, paramGroup, domainMarkers, boundarySegments);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Reader for grid files of type ." + ext);
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
    void passDataFromReader(MeshFileReader& reader, const std::string& paramGroup, bool domainMarkers, bool boundarySegments)
    {
        const auto& bulkGridVertices = reader.bulkGridVertices();
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridPtrTuple_)), [&](const auto id)
        {
            auto& factory = this->gridFactory<id>();

            // insert grid vertices
            if (id == 0)
                for (const auto& v : bulkGridVertices)
                    factory.insertVertex(v);
            else
                for (const auto idx : reader.vertexIndices(id))
                    factory.insertVertex(bulkGridVertices[idx]);

            // insert elements
            const auto& elements = reader.elementData(id);
            for (const auto& e : elements)
                factory.insertElement(e.gt, e.cornerIndices);

            // insert boundary segments
            if (boundarySegments)
            {
                const auto& segments = reader.boundarySegmentData(id);
                for (const auto& bs : segments)
                    factory.insertBoundarySegment(bs);
                std::swap(boundaryMarkerMaps_[id], reader.boundaryMarkerMap(id));
            }

            // maybe copy domain marker map
            if (domainMarkers)
                std::swap(elementMarkerMaps_[id], reader.elementMarkerMap(id));

            // copy the embedments
            std::swap(embeddedEntityMaps_[id], reader.embeddedEntityMap(id));
            std::swap(embedmentMaps_[id], reader.embedmentMap(id));

            // make grid
            std::get<id>(gridPtrTuple_) = std::shared_ptr<Grid<id>>(factory.createGrid());
        });
    }

    using Indices = std::make_index_sequence<Traits::numGrids>;
    //! tuple to store the grids
    using GridPtrTuple = typename makeFromIndexedType<std::tuple, GridPtr, Indices>::type;
    GridPtrTuple gridPtrTuple_;

    //! tuple to store the grid grid factories
    using GridFactoryTuble = typename makeFromIndexedType<std::tuple, GridFactory, Indices>::type;
    GridFactoryTuble gridFactoryTuple_;

    //! data on connectivity between the grids
    std::array< typename Traits::EmbeddedEntityMap, Traits::numGrids > embeddedEntityMaps_;
    std::array< typename Traits::EmbedmentMap, Traits::numGrids > embedmentMaps_;

    //! data on domain and boundary markers
    std::array< typename Traits::ElementToDomainMarkerMap, Traits::numGrids > elementMarkerMaps_;
    std::array< typename Traits::BoundarySegmentToMarkerMap, Traits::numGrids > boundaryMarkerMaps_;
};

} // end namespace Dumux

#endif // DUMUX_FACETCOUPLING_GRID_CREATOR_HH
