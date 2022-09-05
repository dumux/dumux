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
 * \brief Class that exposes the residual of a PDE as linearizable function.
 */
#ifndef DUMUX_ASSEMBLY_RESIDUAL_FUNCTION_HH
#define DUMUX_ASSEMBLY_RESIDUAL_FUNCTION_HH

#include <memory>
#include <utility>
#include <functional>
#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>
#include <dumux/experimental/new_assembly/dumux/common/linearization.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

namespace Dumux {
namespace Concepts {

template<typename T, typename J, typename R>
concept Assembler = requires(const T& t, J& jacobian, R& res) {
    typename T::Variables;

    { t.assembleResidual(res, std::declval<const typename T::Variables&>()) };
    {
        t.assembleJacobianAndResidual(jacobian,
                                      res,
                                      std::declval<const typename T::Variables&>())
    };
};

} // namespace Concepts

/*!
 * \ingroup Assembly
 * \brief Class that exposes the residual of a PDE as linearizable function.
 * \tparam Residual The range type of the PDE to hold the residual
 * \tparam Jacobian The type used for the Jacobian matrix of the PDE system
 * \tparam Assembler Object that assembles the residual and the Jacobian
 * \tparam Derivative An object to represent a derivative. Must be constructible from
 *                    a reference to a Jacobian. The default is a reference to the Jacobian.
 */
template<LinearSystem::Interoperable Residual,
         LinearSystem::Interoperable Jacobian,
         Concepts::Assembler<Jacobian, Residual> Assembler,
         typename Derivative = Auto>
class ResidualFunction
{
    using Variables = typename Assembler::Variables;

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
    ResidualFunction& operator=(const ResidualFunction&) = delete;
    ResidualFunction(const ResidualFunction&) = delete;

    // moves are fine because the ownership is transferred
    ResidualFunction& operator=(ResidualFunction&&) = default;
    ResidualFunction(ResidualFunction&&) = default;

    //! Constructor with externally managed jacobian/residual
    ResidualFunction(std::shared_ptr<const Assembler> assembler,
                     std::shared_ptr<Jacobian> jacobian,
                     std::shared_ptr<Residual> residual)
    : assembler_(std::move(assembler))
    , jacobian_(std::move(jacobian))
    , residual_(std::move(residual))
    {}

    //! Constructor taking ownership over jacobian/residual
    ResidualFunction(std::shared_ptr<const Assembler> assembler,
                     Jacobian&& jacobian,
                     Residual&& residual)
    : assembler_(std::move(assembler))
    , jacobian_(std::make_shared<Jacobian>(std::move(jacobian)))
    , residual_(std::make_shared<Residual>(std::move(residual)))
    {}

    //! Evaluate the residual at the given evaluation point
    const Residual& evaluateAt(const Variables& vars)
    {
        LinearSystem::fill(*residual_, 0.0);
        this->assembler_->assembleResidual(
            *residual_,
            vars
        );
        return *residual_;
    }

    //! Linearize the residual function at the given evaluation point
    Linearization linearizeAt(const Variables& vars)
    {
        LinearSystem::fill(*residual_, 0.0);
        LinearSystem::fill(*jacobian_, 0.0);
        assembler_->assembleJacobianAndResidual(
            *jacobian_,
            *residual_,
            vars
        );
        derivative_ = std::make_unique<DerivativeStorage>(*jacobian_);
        return {*derivative_, *residual_};
    }

private:
    std::shared_ptr<const Assembler> assembler_;
    std::shared_ptr<Jacobian> jacobian_;
    std::shared_ptr<Range> residual_;
    std::unique_ptr<DerivativeStorage> derivative_{nullptr};
};

} // namespace Dumux

#endif
