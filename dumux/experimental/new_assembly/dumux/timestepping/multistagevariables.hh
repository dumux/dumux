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
 * \ingroup TimeStepping
 * \brief Variables class usable in the context of multi-stage time integration.
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_VARIABLES_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_VARIABLES_HH

#include <memory>
#include <utility>
#include <cassert>
#include <type_traits>
#include <concepts>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/variables.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/multistageparams.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<typename V> requires(std::copy_constructible<std::decay_t<V>>)
class MultiStageVariablesImpl
{
    using SV = std::decay_t<V>;
    using Scalar = Dumux::Variables::ScalarType<SV>;

public:
    using Dofs = typename SV::Dofs;
    using TimeLevel = typename SV::TimeLevel;

    using StageVariables = SV;
    using StageParams = MultiStageParams<Scalar>;

    template<std::convertible_to<V> _V>
    MultiStageVariablesImpl(_V&& vars)
    : vars_(std::forward<_V>(vars))
    {}

    decltype(auto) timeLevel() const
    { return vars_.timeLevel(); }

    decltype(auto) dofs() const
    { return vars_.dofs(); }

    void update(const Dofs& dofs)
    { vars_.update(dofs); }

    void update(const Dofs& dofs, const TimeLevel& timeLevel)
    { vars_.update(dofs, timeLevel); }

    void updateTime(const TimeLevel& timeLevel)
    { vars_.updateTime(timeLevel); }

    std::size_t numStages() const
    { return prevStageVariables_.size() + 1; }

    void setStageParams(std::shared_ptr<const StageParams> p)
    {
        assert(p->size() == numStages() + 1);
        params_ = p;
        prevStageVariables_.push_back(vars_);
        updateTime_(vars_, params_->size() - 1);
    }

    void clearStages()
    {
        params_ = nullptr;
        prevStageVariables_.clear();
    }

    const StageVariables& stageVariables(std::size_t i) const
    {
        assert(i <= numStages());
        return i < numStages() - 1 ? prevStageVariables_[i] : vars_;
    }

    const StageParams& stageParams() const
    {
        assert(params_);
        return *params_;
    }

    explicit(false) operator StageVariables&() { return vars_; }
    explicit(false) operator const StageVariables&() const { return vars_; }

private:
    void updateTime_(StageVariables& vars, std::size_t stageIdx)
    {
        assert(params_->size() > stageIdx);
        vars.updateTime(
            TimeLevel{
                params_->timeAtStage(stageIdx),
                params_->timeAtStage(0),
                params_->timeStepFraction(stageIdx)
            }
        );
    }

    V vars_;
    std::vector<StageVariables> prevStageVariables_;
    std::shared_ptr<const StageParams> params_{nullptr};
};

template<typename T>
concept DerivedFromMultiStageVariablesImpl
    = std::derived_from<T, MultiStageVariablesImpl<typename T::StageVariables>>
    or std::derived_from<T, MultiStageVariablesImpl<typename T::StageVariables&>>;

} // namespace Detail
#endif // DOXYGEN


/*!
 * \ingroup TimeStepping
 * \brief Variables class usable in the context of multi-stage time integration.
 */
template<Concepts::TimeDependentVariables V>
class MultiStageVariables
: public Detail::MultiStageVariablesImpl<std::decay_t<V>>
{
    using ParentType = Detail::MultiStageVariablesImpl<std::decay_t<V>>;

public:
    using typename ParentType::StageVariables;

    template<typename... Args> requires(
        std::constructible_from<StageVariables, Args...>)
    MultiStageVariables(Args&&... args)
    : ParentType(StageVariables{std::forward<Args>(args)...})
    {}

    MultiStageVariables(StageVariables&& vars)
    : ParentType(std::move(vars))
    {}
};


/*!
 * \ingroup TimeStepping
 * \brief Non-owning variables class usable in the context of multi-stage time integration.
 * \note This implementation does not take ownership, but is solely a view on the given
 *       variables object. Thus, its lifetime is bound to the one provided upon construction.
 */
template<Concepts::TimeDependentVariables V>
class MultiStageVariablesView
: public Detail::MultiStageVariablesImpl<V&>
{
    using ParentType = Detail::MultiStageVariablesImpl<V&>;

public:
    MultiStageVariablesView(V& vars)
    : ParentType(vars)
    {}
};


//! Register `MultiStageVariablesView` as view
template<typename V>
struct Traits::IsView<MultiStageVariablesView<V>>
: public std::true_type {};

} // namespace Dumux

#endif
