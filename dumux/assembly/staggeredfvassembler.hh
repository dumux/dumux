// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \brief A linear system assembler (residual and Jacobian) for staggered finite volume schemes
 */
#ifndef DUMUX_STAGGERED_FV_ASSEMBLER_HH
#define DUMUX_STAGGERED_FV_ASSEMBLER_HH

#include <type_traits>

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/method.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/staggeredtraits.hh>
#include <dumux/multidomain/staggeredcouplingmanager.hh>

#include "diffmethod.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \brief A linear system assembler (residual and Jacobian) for staggered finite volume schemes.
 *        This is basically just a wrapper for the MultiDomainFVAssembler which simplifies the set-up
 *        of uncoupled problems using the staggered scheme.
 * \tparam TypeTag the TypeTag
 * \tparam diffMethod the differentiation method to residual compute derivatives
 * \tparam isImplicit if to use an implicit or explicit time discretization
 */
template<class TypeTag, DiffMethod diffMethod, bool isImplicit = true>
class StaggeredFVAssembler: public MultiDomainFVAssembler<StaggeredMultiDomainTraits<TypeTag, TypeTag>,
                                                          StaggeredCouplingManager<StaggeredMultiDomainTraits<TypeTag, TypeTag>>,
                                                          diffMethod>
{
    using ParentType = MultiDomainFVAssembler<StaggeredMultiDomainTraits<TypeTag, TypeTag>,
                                              StaggeredCouplingManager<StaggeredMultiDomainTraits<TypeTag, TypeTag>>,
                                              diffMethod>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using TimeLoop = TimeLoopBase<GetPropType<TypeTag, Properties::Scalar>>;

public:
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using CouplingManager = typename ParentType::CouplingManager;

    //! The constructor for stationary problems
    StaggeredFVAssembler(std::shared_ptr<const Problem> problem,
                         std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                         std::shared_ptr<GridVariables> gridVariables)
    : ParentType(std::make_tuple(problem, problem),
                 std::make_tuple(fvGridGeometry->cellCenterFVGridGeometryPtr(), fvGridGeometry->faceFVGridGeometryPtr()),
                 std::make_tuple(gridVariables->cellCenterGridVariablesPtr(), gridVariables->faceGridVariablesPtr()),
                 std::make_shared<CouplingManager>())
    {
        static_assert(isImplicit, "Explicit assembler for stationary problem doesn't make sense!");
        this->couplingManager_->setSubProblems(std::make_tuple(problem, problem));
    }

    //! The constructor for instationary problems
    StaggeredFVAssembler(std::shared_ptr<const Problem> problem,
                         std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                         std::shared_ptr<GridVariables> gridVariables,
                         std::shared_ptr<const TimeLoop> timeLoop)
    : ParentType(std::make_tuple(problem, problem),
                 std::make_tuple(fvGridGeometry->cellCenterFVGridGeometryPtr(), fvGridGeometry->faceFVGridGeometryPtr()),
                 std::make_tuple(gridVariables->cellCenterGridVariablesPtr(), gridVariables->faceGridVariablesPtr()),
                 std::make_shared<CouplingManager>(),
                 timeLoop)
    {
        static_assert(isImplicit, "Explicit assembler for stationary problem doesn't make sense!");
        this->couplingManager_->setSubProblems(std::make_tuple(problem, problem));
    }


    auto& gridVariables()
    { return ParentType::gridVariables(Dune::index_constant<0>()); }

    const auto& gridVariables() const
    { return ParentType::gridVariables(Dune::index_constant<0>()); }

    //! The global finite volume geometry
    const FVGridGeometry& fvGridGeometry() const
    { return ParentType::fvGridGeometry(Dune::index_constant<0>()).actualfvGridGeometry(); }

};

} // namespace Dumux

#endif
