// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \brief A linear system assembler (residual and Jacobian) for staggered finite volume schemes
 */
#ifndef DUMUX_STAGGERED_FV_ASSEMBLER_HH
#define DUMUX_STAGGERED_FV_ASSEMBLER_HH

#include <tuple>
#include <memory>

#include <dune/common/indices.hh>

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
class StaggeredFVAssembler : public MultiDomainFVAssembler<StaggeredMultiDomainTraits<TypeTag, TypeTag>,
                                                           StaggeredCouplingManager<StaggeredMultiDomainTraits<TypeTag, TypeTag>>,
                                                           diffMethod, isImplicit>
{
    using ParentType = MultiDomainFVAssembler<StaggeredMultiDomainTraits<TypeTag, TypeTag>,
                                              StaggeredCouplingManager<StaggeredMultiDomainTraits<TypeTag, TypeTag>>,
                                              diffMethod, isImplicit>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using TimeLoop = TimeLoopBase<GetPropType<TypeTag, Properties::Scalar>>;

public:
    using typename ParentType::ResidualType;
    using typename ParentType::SolutionVector;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using CouplingManager = typename ParentType::CouplingManager;

    //! The constructor for stationary problems
    StaggeredFVAssembler(std::shared_ptr<const Problem> problem,
                         std::shared_ptr<const GridGeometry> gridGeometry,
                         std::shared_ptr<GridVariables> gridVariables)
    : ParentType(std::make_tuple(problem, problem),
                 std::make_tuple(gridGeometry->faceFVGridGeometryPtr(), gridGeometry->cellCenterFVGridGeometryPtr()),
                 std::make_tuple(gridVariables->faceGridVariablesPtr(), gridVariables->cellCenterGridVariablesPtr()),
                 std::make_shared<CouplingManager>())
    {
        static_assert(isImplicit, "Explicit assembler for stationary problem doesn't make sense!");
        this->couplingManager_->setSubProblems(std::make_tuple(problem, problem));
    }

    //! The constructor for time-dependent problems
    StaggeredFVAssembler(std::shared_ptr<const Problem> problem,
                         std::shared_ptr<const GridGeometry> gridGeometry,
                         std::shared_ptr<GridVariables> gridVariables,
                         std::shared_ptr<const TimeLoop> timeLoop,
                         const SolutionVector& prevSol)
    : ParentType(std::make_tuple(problem, problem),
                 std::make_tuple(gridGeometry->faceFVGridGeometryPtr(), gridGeometry->cellCenterFVGridGeometryPtr()),
                 std::make_tuple(gridVariables->faceGridVariablesPtr(), gridVariables->cellCenterGridVariablesPtr()),
                 std::make_shared<CouplingManager>(),
                 timeLoop,
                 prevSol)
    {
        this->couplingManager_->setSubProblems(std::make_tuple(problem, problem));
    }

    auto& gridVariables()
    { return ParentType::gridVariables(Dune::index_constant<0>()); }

    const auto& gridVariables() const
    { return ParentType::gridVariables(Dune::index_constant<0>()); }

    const GridGeometry& gridGeometry() const
    { return ParentType::gridGeometry(Dune::index_constant<0>()).actualGridGeometry(); }

};

} // namespace Dumux

#endif
