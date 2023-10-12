// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief The infrastructure to retrieve metadata information.
 */
#ifndef DUMUX_COMMON_METADATA_HH
#define DUMUX_COMMON_METADATA_HH

#include <iostream>
#include <list>
#include <sstream>
#include <unordered_map>
#include <fstream>
#include <functional>
#include <string>
#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <dune/common/concept.hh>
#include <dune/common/classname.hh>
#include <dune/grid/common/gridview.hh>

#include <dumux/io/json.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/common/typetraits/isvalid.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/fvgridvariables.hh>

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
 * \ingroup Core
 * \brief Class to collect metadata
 */
class Collector
{

    using JsonTree = Dumux::Json::JsonTree;

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
     * \brief Merges two trees by overwriting existing values
     */
    void merge(const Collector& collector)
    {
        this->getTree().merge_patch(collector.getTree());
    }

    /*!
     * \brief Append data from another collector
     * \param collector The json collector from which data is taken
     * \param convertToArrays Convert non-array types to array which allows appending data
     */
    void append(const Collector& collector, bool convertToArrays = false)
    {
        const auto& tree = collector.getTree();
        for (const auto& [key, values] : tree.items())
        {
            auto& dataAtKey = this->getTree()[key];
            if(dataAtKey.is_array())
            {
                if(values.is_array())
                    dataAtKey.insert(dataAtKey.end(), values.begin(), values.end());
                else
                    dataAtKey.push_back(values);
            }
            else if(dataAtKey.is_null())
            {
                dataAtKey = values;
            }
            else if(convertToArrays)
            {
                // convert to array and append data
                auto val(dataAtKey);
                dataAtKey = JsonTree::array({val});
                if(values.is_array())
                    dataAtKey.insert(dataAtKey.end(), values.begin(), values.end());
                else
                    dataAtKey.push_back(values);
            }
            else
                DUNE_THROW(Dune::InvalidStateException, "Unclear how to append data without conversion to array!");
        }
    }

    /*!
     * \brief returns the object with id of the json tree
     */
    auto& operator[] (const std::string& id)
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

//! convenience function to check if file exists
bool jsonFileExists(const std::string& fileName)
{
    std::ifstream infile(fileName + ".json");
    return infile.good();
}

//! reads a json file into a tree
template<class Collector>
void readJsonFile(Collector& collector, const std::string& fileName)
{
    std::ifstream i(fileName + ".json");
    i >> collector.getTree();
}

//! writes a json tree to file
template<class Collector>
void writeJsonFile(const Collector& collector, const std::string& fileName)
{
    std::ofstream o(fileName + ".json");
    o << std::setw(4) << collector.getTree() << std::endl;
}

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
    obj["Type"] = Collector::className(a, hideTemplates);
    obj["Stationary"] = a.isStationaryProblem();
}

template<class Collector, class GridGeometry>
auto collectMetaData(Collector& collector, const GridGeometry& gg, bool hideTemplates = true)
-> typename std::enable_if_t<Dune::models<Concept::GridGeometry, GridGeometry>()>
{
    auto& obj = collector["GridGeometry"];
    obj["Type"] = Collector::className(gg, hideTemplates);
    obj["IsPeriodic"] = gg.isPeriodic();
    obj["DiscretizationMethod"] = GridGeometry::discMethod.name();
    obj["NumScvs"] = gg.numScv();
    obj["NumScvfs"] = gg.numScvf();
    obj["NumBoundaryScvfs"] = gg.numBoundaryScvf();
    obj["NumDofs"] = gg.numDofs();
}

template<class Collector, class GridVariables>
auto collectMetaData(Collector& collector, const GridVariables& gv, bool hideTemplates = true)
-> typename std::enable_if_t<Dune::models<Concept::GridVariables, GridVariables>()>
{
    auto& obj = collector["GridVariables"];
    obj["Type"] = Collector::className(gv, hideTemplates);
    obj["GridVolumeVariables"]["Type"] = Collector::template className<typename GridVariables::GridVolumeVariables>(hideTemplates);
    obj["VolumeVariables"]["Type"] = Collector::template className<typename GridVariables::VolumeVariables>(hideTemplates);
    obj["GridFluxVariablesCache"]["Type"] = Collector::template className<typename GridVariables::GridFluxVariablesCache>(hideTemplates);
}

template<class Collector, class GridView>
auto collectMetaData(Collector& collector, const GridView& gridView, bool hideTemplates = true)
-> typename std::enable_if_t<Dune::models<Concept::GridView, GridView>()>
{
    auto& obj = collector["GridView"];
    obj["Type"] = Collector::className(gridView, hideTemplates);
    obj["dimension"] = GridView::dimension;
    obj["dimensionWorld"] = GridView::dimensionworld;
    obj["conforming"] = GridView::conforming;
    //obj["Grid"]["Type"] = Collector::className(gridView.grid(), hideTemplates);
    for(int codim = 0; codim < GridView::dimension; ++codim)
       obj["numEntities"]["codim " + std::to_string(codim) ] = gridView.size(codim);

    // TODO parallel runs, i.e. overlapSize() etc.
}

template<class TypeTag, class Collector>
auto collectTypeTags(Collector& collector)
{
    auto& obj = collector["TTags"];
    obj = Dumux::Json::JsonTree::array();
    Detail::collectTypeTagsFromTuple<std::tuple<TypeTag>>(obj);
}

} // end namespace Dumux::MetaData

#endif
