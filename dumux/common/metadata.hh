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

#include <nlohmann/json.hpp>

#include <dumux/common/typetraits/isvalid.hh>

#include "dumux/assembly/fvassembler.hh"
#include "dumux/assembly/diffmethod.hh"
#include "dumux/discretization/basegridgeometry.hh"
#include "dumux/discretization/fvgridvariables.hh"

namespace Dumux {

namespace Detail {

//! Helper to determine whether a given type inherits from BaseGridGeometry
struct isGridGeometry
{
    template<typename ...Args>
    void isConstructable(const BaseGridGeometry<Args...>&)
    {}

    template<class GridGeometry>
    auto operator()(GridGeometry&& gg)
    -> decltype(isConstructable(gg))
    {}
};

//! Helper to determine whether a given type inherits from FVGridVariables
struct isGridVariables
{
    template<typename ...Args>
    void isConstructable(const FVGridVariables<Args...>&)
    {}

    template<class GridVariables>
    auto operator()(GridVariables&& gv)
    -> decltype(isConstructable(gv))
    {}
};

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

    template <class T>
    static std::string className(const T& c, bool hideTemplates)
    {
        return hideTemplates ? hideTemplateArguments(Dune::className(c)) : Dune::className(c);
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
    auto& obj = collector.getTree()["Assembler"];
    obj["Type"] = Metadata::className(a, hideTemplates);
    obj["Stationary"] = a.isStationaryProblem();
}

template<class Collector, class GridGeometry>
auto collectMetaData(Collector& collector, const GridGeometry& gg, bool hideTemplates = true)
-> typename std::enable_if_t<decltype(isValid(Detail::isGridGeometry())(gg))::value, void>
{
    auto& obj = collector.getTree()["GridGeometry"];
    obj["Type"] = Metadata::className(gg, hideTemplates);
    obj["GridView"]["Type"] = Metadata::className(gg.gridView(), hideTemplates);
    obj["IsPeriodic"] = gg.isPeriodic();
    obj["DiscretisationMethod"] = toString(GridGeometry::discMethod);
    obj["MaxElementStencilSize"] = GridGeometry::maxElementStencilSize;
    obj["NumScvs"] = gg.numScv();
    obj["NumScvfs"] = gg.numScvf();
    obj["SumBoundaryScvfs"] = gg.numBoundaryScvf();
    obj["NumDofs"] = gg.numDofs();
}

template<class Collector, class GridVariables>
auto collectMetaData(Collector& collector, const GridVariables& gv, bool hideTemplates = true)
-> typename std::enable_if_t<decltype(isValid(Detail::isGridVariables())(gv))::value, void>
{
    auto& obj = collector.getTree()["GridVariables"];
    obj["Type"] = Metadata::className(gv, hideTemplates);
}

} // end namespace Dumux

#endif
