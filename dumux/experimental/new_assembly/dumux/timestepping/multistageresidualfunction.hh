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
 * \ingroup Assembly
 * \brief Class that exposes the residual of a multi-stage PDE as a linearizable function.
 */
#ifndef DUMUX_TIME_STEPPING_MULTI_STAGE_RESIDUAL_FUNCTION_HH
#define DUMUX_TIME_STEPPING_MULTI_STAGE_RESIDUAL_FUNCTION_HH

#include <memory>
#include <utility>
#include <functional>

#include <dumux/experimental/new_assembly/dumux/common/variables.hh>
#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>
#include <dumux/experimental/new_assembly/dumux/common/linearization.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/operatorweights.hh>

#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

#include <dumux/experimental/new_assembly/dumux/timestepping/multistagevariables.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/multistageparams.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<typename Vars>
decltype(auto) operatorWeights(const Vars&)
{ return Dumux::OperatorWeights<Variables::ScalarType<Vars>>{}; }

} // namespace Detail
#endif // DOXYGEN


namespace Concepts {

template<typename T, typename Jacobian, typename Residual>
concept StageAssembler = requires  {
    typename T::Variables;
    Variables<typename T::Variables>;
} and requires(const T& t, Jacobian& jacobian, Residual& res, const typename T::Variables& vars) {
    { t.assembleResidual(res, vars, Detail::operatorWeights(vars)) };
    { t.assembleJacobianAndResidual(jacobian, res, vars, Detail::operatorWeights(vars)) };
};

} // namespace Concepts


/*!
 * \ingroup Assembly
 * \brief Class that exposes the residual of a multi-stage PDE as linearizable function.
 * \tparam Residual The range type of the PDE to hold the residual
 * \tparam Jacobian The type used for the Jacobian matrix of the PDE system
 * \tparam Assembler Object that assembles the residual and the Jacobian (per stage).
 * \tparam Derivative An object to represent a derivative. Must be constructible from
 *                    a reference to a Jacobian. The default is a reference to the Jacobian.
 * \note This implementation wraps the given variables in `MultiStageVariables`,
 *       which additionaly expose the variables of previous stages.
 */
template<LinearSystem::Interoperable Residual,
         LinearSystem::Interoperable Jacobian,
         Concepts::StageAssembler<Jacobian, Residual> Assembler,
         typename Derivative = Auto>
class MultiStageResidualFunction
{
    using Variables = Dumux::MultiStageVariablesView<typename Assembler::Variables>;

    static constexpr bool defaultDeriv = std::same_as<Derivative, Auto>;
    using DerivativeType = std::conditional_t<defaultDeriv, Jacobian, Derivative>;
    using DerivativeStorage = std::conditional_t<defaultDeriv, std::reference_wrapper<Jacobian>, Derivative>;
    static_assert(std::constructible_from<DerivativeStorage, Jacobian&>,
                  "Derivative type is required to be constructible from a reference to a Jacobian");

public:
    using Domain = Variables;
    using Range = Residual;
    using Linearization = Dumux::Linearization<DerivativeType, Residual>;

    // copies are dangerous because they would operate on the same residual/jacobian
    MultiStageResidualFunction& operator=(const MultiStageResidualFunction&) = delete;
    MultiStageResidualFunction(const MultiStageResidualFunction&) = delete;

    // moves are fine because the ownership is transferred
    MultiStageResidualFunction& operator=(MultiStageResidualFunction&&) = default;
    MultiStageResidualFunction(MultiStageResidualFunction&&) = default;

    //! Constructor with externally managed jacobian/residual
    MultiStageResidualFunction(std::shared_ptr<const Assembler> assembler,
                               std::shared_ptr<Jacobian> jacobian,
                               std::shared_ptr<Residual> residual)
    : assembler_(std::move(assembler))
    , jacobian_(std::move(jacobian))
    , residual_(std::move(residual))
    {}

    //! Constructor taking ownership over jacobian/residual
    MultiStageResidualFunction(std::shared_ptr<const Assembler> assembler,
                               Jacobian&& jacobian,
                               Residual&& residual)
    : assembler_(std::move(assembler))
    , jacobian_(std::make_shared<Jacobian>(std::move(jacobian)))
    , residual_(std::make_shared<Residual>(std::move(residual)))
    {}

    const Residual& evaluateAt(const Variables& vars)
    {
        LinearSystem::fill(*residual_, 0.0);
        assemblePreviousStagesResiduals_(vars);
        assembleStageResidual_(vars, vars.numStages() - 1);
        return *residual_;
    }

    Linearization linearizeAt(const Variables& vars)
    {
        LinearSystem::fill(*jacobian_, 0.0);
        LinearSystem::fill(*residual_, 0.0);

        const auto& params = vars.stageParams();
        const auto curStageIdx = vars.numStages() - 1;

        assemblePreviousStagesResiduals_(vars);
        assembler_->assembleJacobianAndResidual(
            *jacobian_,
            *residual_,
            vars,
            makeOperatorWeights_(params, curStageIdx)
        );

        derivative_ = std::make_unique<DerivativeStorage>(*jacobian_);
        return {*derivative_, *residual_};
    }

private:
    void assemblePreviousStagesResiduals_(const Variables& vars)
    {
        for (std::size_t i = 0; i < vars.numStages() - 1; ++i)
            assembleStageResidual_(vars, i);
    }

    void assembleStageResidual_(const Variables& vars, const std::size_t stageIdx)
    {
        assembler_->assembleResidual(
            *residual_,
            vars.stageVariables(stageIdx),
            makeOperatorWeights_(vars.stageParams(), stageIdx)
        );
    }

    template<typename S>
    auto makeOperatorWeights_(const Dumux::MultiStageParams<S>& params, std::size_t stageIdx) const
    {
        using Scalar = Dumux::Variables::ScalarType<Variables>;
        OperatorWeights<Scalar> weights;
        if (!params.skipTemporal(stageIdx))
            weights.temporalWeight = params.temporalWeight(stageIdx);
        if (!params.skipSpatial(stageIdx))
            weights.spatialWeight = params.spatialWeight(stageIdx);
        return weights;
    }

    std::shared_ptr<const Assembler> assembler_;
    std::shared_ptr<Jacobian> jacobian_;
    std::shared_ptr<Range> residual_;
    std::unique_ptr<DerivativeStorage> derivative_{nullptr};
};

} // namespace Dumux

#endif
