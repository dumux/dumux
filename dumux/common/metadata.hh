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
 * \ingroup Common
 * \brief The infrastructure to retrieve metadata information.
 */
#ifndef DUMUX_METADATA_HH
#define DUMUX_METADATA_HH

#include <iostream>
#include <list>
#include <sstream>
#include <unordered_map>
#include <fstream>
#include <functional>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <dune/common/concept.hh>
#include <dune/grid/common/gridview.hh>

#include "dumux/io/json/json.hpp"

#include "dumux/common/properties/propertysystem.hh"
#include "dumux/common/typetraits/utility.hh"
#include "dumux/common/typetraits/isvalid.hh"

#include "dumux/assembly/fvassembler.hh"
#include "dumux/assembly/diffmethod.hh"
#include "dumux/discretization/basegridgeometry.hh"
#include "dumux/discretization/fvgridvariables.hh"

namespace Dumux::MetaData {

namespace Concept {

//! Concept of GridGeometry
struct GridGeometry
{
    template<class GG>
    auto require(const GG& gg) -> decltype(
        gg.isPeriodic(),
        gg.numScv(),
        gg.numScvf(),
        gg.numBoundaryScvf(),
        gg.numDofs(),
        GG::discMethod
    );
};

//! Concept of GridVariables
struct GridVariables
{
    template<class GV>
    auto require(const GV& gv) -> decltype(
        Dune::Concept::requireType<typename GV::GridVolumeVariables>(),
        Dune::Concept::requireType<typename GV::VolumeVariables>(),
        Dune::Concept::requireType<typename GV::GridFluxVariablesCache>()
    );
};

//! Concept of GridView
struct GridView
{
    template<class GV>
    auto require(const GV& gv) -> decltype(
        Dune::Concept::requireBaseOf<Dune::GridView<typename GV::Traits>, GV>()
    );
};

} // end namespace Concept

namespace Detail {

    std::string removeNamespace(std::string&& s)
    {
        std::size_t last = s.find_last_of("::");

        if(last != std::string::npos)
            s.erase(0, last+1);

        return std::move(s);
    }

    template<class TTagTuple, class Collector>
    void collectTypeTagsFromTuple(Collector& collector, int depth=0, int parentBranch=-1)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<std::tuple_size_v<TTagTuple>>{},  [&](auto i)
        {
            using type = typename std::tuple_element<i, TTagTuple>::type;
            collector.push_back(std::tuple<int, int, std::string>{depth, parentBranch, removeNamespace(Dune::className<type>())});
            if constexpr (Dumux::Properties::Detail::hasParentTypeTag<type>(int{}))
                collectTypeTagsFromTuple<typename type::InheritsFrom>(collector, int{depth+1}, i);
        });
    }

} // end namespace Detail

/*!
 * \ingroup Common
 * \brief Class to collect metadata
 * \todo Doc me!
 */
class Metadata {

    using JsonTree = nlohmann::json;

public:
    /*!
     * \brief Get the json tree
     */
    JsonTree& getTree()
    {
        return tree_;
    }

    const JsonTree& getTree() const
    {
        return tree_;
    }

    /*!
     * \brief returns the object with id of the json tree
     */
    auto& operator[] (std::string id)
    { return getTree()[id]; }

    template <class T>
    static std::string className(const T& c, bool hideTemplates)
    {
        return hideTemplates ? hideTemplateArguments(Dune::className(c)) : Dune::className(c);
    }

    template <class T>
    static std::string className(bool hideTemplates)
    {
        return hideTemplates ? hideTemplateArguments(Dune::className<T>()) : Dune::className<T>();
    }

    static std::string hideTemplateArguments(std::string&& s)
    {
        std::size_t first = s.find("<");
        std::size_t last = s.find_last_of(">");

        if(first != std::string::npos && last != std::string::npos)
            s.replace(first, last-first+1, "<...>");

        s.erase(std::unique(std::begin(s), std::end(s),
                [](unsigned char a, unsigned char b){return std::isspace(a) && std::isspace(b);}), std::end(s));

        return std::move(s);
    }

private:
    JsonTree tree_;
};

//! prints json tree
template<class Collector>
void print(const Collector& collector)
{
    std::cout << collector.getTree().dump(4) << std::endl;
}

template<class Collector, class TypeTag, DiffMethod diffmethod, bool isImplicit>
void collectMetaData(Collector& collector, const FVAssembler<TypeTag, diffmethod, isImplicit>& a, bool hideTemplates = true)
{
    auto& obj = collector["Assembler"];
    obj["Type"] = Metadata::className(a, hideTemplates);
    obj["Stationary"] = a.isStationaryProblem();
}

template<class Collector, class GridGeometry>
auto collectMetaData(Collector& collector, const GridGeometry& gg, bool hideTemplates = true)
-> typename std::enable_if_t<Dune::models<Concept::GridGeometry, GridGeometry>()>
{
    auto& obj = collector["GridGeometry"];
    obj["Type"] = Metadata::className(gg, hideTemplates);
    obj["IsPeriodic"] = gg.isPeriodic();
    obj["DiscretisationMethod"] = toString(GridGeometry::discMethod);
    obj["NumScvs"] = gg.numScv();
    obj["NumScvfs"] = gg.numScvf();
    obj["SumBoundaryScvfs"] = gg.numBoundaryScvf();
    obj["NumDofs"] = gg.numDofs();
}

template<class Collector, class GridVariables>
auto collectMetaData(Collector& collector, const GridVariables& gv, bool hideTemplates = true)
-> typename std::enable_if_t<Dune::models<Concept::GridVariables, GridVariables>()>
{
    auto& obj = collector["GridVariables"];
    obj["Type"] = Metadata::className(gv, hideTemplates);
    obj["GridVolumeVariables"]["Type"] = Metadata::className<typename GridVariables::GridVolumeVariables>(hideTemplates);
    obj["VolumeVariables"]["Type"] = Metadata::className<typename GridVariables::VolumeVariables>(hideTemplates);
    obj["GridFluxVariablesCache"]["Type"] = Metadata::className<typename GridVariables::GridFluxVariablesCache>(hideTemplates);
}

template<class Collector, class GridView>
auto collectMetaData(Collector& collector, const GridView& gridView, bool hideTemplates = true)
-> typename std::enable_if_t<Dune::models<Concept::GridView, GridView>()>
{
    auto& obj = collector["GridView"];
    obj["Type"] = Metadata::className(gridView, hideTemplates);
    obj["dimension"] = GridView::dimension;
    obj["dimensionWorld"] = GridView::dimensionworld;
    obj["conforming"] = GridView::conforming;
    //obj["Grid"]["Type"] = Metadata::className(gridView.grid(), hideTemplates);
    for(int codim = 0; codim < GridView::dimension; ++codim)
       obj["numEntities"]["codim " + std::to_string(codim) ] = gridView.size(codim);

    //TODO parallel runs, i.e. overlapSize() etc.
}

template<class TypeTag, class Collector>
auto collectTypeTags(Collector& collector)
{
    auto& obj = collector["TTags"];
    obj = nlohmann::json::array();
    Detail::collectTypeTagsFromTuple<std::tuple<TypeTag>>(obj);
}

} // end namespace Dumux

#endif
