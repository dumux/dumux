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

    template<class TypeTag, DiffMethod diffmethod, bool isImplicit>
    static void collectMetaData(const FVAssembler<TypeTag, diffmethod, isImplicit>& a)
    {
        auto& obj = getTree()["Assembler"];
        obj["Type"] = Dune::className(a);
        obj["Stationary"] = a.isStationaryProblem() ? "true" : "false";
    }

    template<class GridGeometry>
    static auto collectMetaData(const GridGeometry& gg)
    -> typename std::enable_if_t<decltype(isValid(Detail::isGridGeometry())(gg))::value, void>
    {
        auto& obj = getTree()["GridGeometry"];
        obj["Type"] = Dune::className(gg);
        obj["GridView"] = Dune::className(gg.gridView());
        obj["IsPeriodic"] = gg.isPeriodic() ? "true" : "false";
    }

    template<class GridVariables>
    static auto collectMetaData(const GridVariables& gv)
    -> typename std::enable_if_t<decltype(isValid(Detail::isGridVariables())(gv))::value, void>
    {
        auto& obj = getTree()["GridVariables"];
        obj["Type"] = Dune::className(gv);
    }

    //! prints json tree
    static void print()
    {
        std::cout << getTree().dump(4) << std::endl;
    }

    /*!
     * \brief Get the json tree
     */
    static JsonTree& getTree()
    {
        static JsonTree parser_;
        return parser_;
    }

};

} // end namespace Dumux

#endif
