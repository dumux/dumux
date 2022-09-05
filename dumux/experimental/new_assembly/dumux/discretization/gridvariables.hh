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
 * \brief Variables class for PDEs defined on computational grids.
 */
#ifndef DUMUX_DISCRETIZATION_GRID_VARIABLES_HH
#define DUMUX_DISCRETIZATION_GRID_VARIABLES_HH

#include <memory>
#include <utility>

#include <dumux/experimental/new_assembly/dumux/common/defaultvariables.hh>
#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>
#include <dumux/experimental/new_assembly/dumux/common/indexstrategies.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Base class for grid variables.
 * \tparam GG The grid geometry type
 * \tparam X The coefficient vector for the degrees of freedom
 * \tparam IS The indexing strategy for accessing the coefficient vector
 */
template<typename GG, typename X, typename IS>
class GridVariables : public DefaultVariables<X>
{
    using ParentType = DefaultVariables<X>;

public:
    using GridGeometry = GG;
    using IndexStrategy = IS;

    template<typename... Args>
    GridVariables(std::shared_ptr<const GridGeometry> gridGeometry,
                  std::shared_ptr<const IndexStrategy> indexStrategy,
                  Args&&... args)
    : ParentType(std::forward<Args>(args)...)
    , gridGeometry_(gridGeometry)
    , indexStrategy_(indexStrategy)
    {}

    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    template<Concepts::MultiIndex MI> requires(
        Concepts::IndexStrategy<IS, MI>)
    decltype(auto) getDofIndex(const MI& multiIndex) const
    { return (*indexStrategy_)[multiIndex]; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    std::shared_ptr<const IndexStrategy> indexStrategy_;
};

} // namespace Dumux

#endif
