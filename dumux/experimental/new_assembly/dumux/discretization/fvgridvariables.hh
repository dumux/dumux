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
 * \brief Grid variables for finite-volume schemes.
 */
#ifndef DUMUX_DISCRETIZATION_FV_GRID_VARIABLES_HH
#define DUMUX_DISCRETIZATION_FV_GRID_VARIABLES_HH

#include <memory>
#include <utility>

#include <dune/istl/bvector.hh>

#include <dumux/common/enumerate.hh>
#include <dumux/common/reservedblockvector.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/gridvariables.hh>

namespace Dumux {

template<typename GridVariables>
class FVGridVariablesLocalView
{
    using Model = typename GridVariables::Model;
    using VolumeVariables = typename Model::VolumeVariables;

    using GridGeometry = typename GridVariables::GridGeometry;
    using LocalGridGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    FVGridVariablesLocalView(const GridVariables& gridVariables)
    : gridVariables_(gridVariables)
    {}

    void bind(const LocalGridGeometry& localGridGeometry) &
    {
        // TODO: local ids of entities within local view?
        for (const auto& scv : scvs(localGridGeometry))
        {
            // add inside VolVars
        }

        for (const auto& scv : outsideScvs(localGridGeometry))
        {
            // add outside ovlvars
        }

        for (const auto& scvf : boundaryScvfs(localGridGeometry))
            if (problem.isDirichlet())
            {
                // add boundary VolVars
            }
    }

private:
    const GridVariables& gridVariables_;
};

/*!
 * \ingroup Discretization
 * \brief Grid variables for finite-volume schemes.
 * \tparam M The model to be solved
 * \tparam GG The grid geometry type
 * \tparam X The coefficient vector for the degrees of freedom
 * \tparam IS The indexing strategy for accessing the coefficient vector
 */
template<typename M, typename GG, typename X, typename IS>
class FVGridVariables : public GridVariables<GG, X, IS>
{
    using ParentType = GridVariables<GG, X, IS>;

public:
    using Model = M;
    using ParentType::ParentType;

    template<typename... Args>
    FVGridVariables(std::shared_ptr<const Model> model, Args&&... args)
    : ParentType(std::forward<Args>(args)...)
    , model_(model)
    {}

private:
    std::shared_ptr<const Model> model_;
};

} // namespace Dumux

#endif
