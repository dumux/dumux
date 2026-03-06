// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief free functions for the evaluation of primary/secondary variables inside elements.
 */
#ifndef DUMUX_DISCRETIZATION_EVALUATION_HH
#define DUMUX_DISCRETIZATION_EVALUATION_HH

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/typetraits.hh>
#include <dune/common/math.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/tag.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::CVFE::Detail {

//! \brief Return the number of components of a scalar or vector-valued quantity.
template<class T>
inline constexpr int numComponents()
{
    if constexpr (Dune::IsNumber<T>::value)
        return 1;
    else
        return T::dimension;
}

//! \brief Expression wrapping a function with an operator tag.
template<class F, class Op>
struct Expr {
    const F& f;
    static constexpr Op op{};
};

template<class T>
struct IsExpr : std::false_type {};

template<class F, class Op>
struct IsExpr<Expr<F, Op>> : std::true_type {};

}

namespace Dumux::EvaluationOperators {
/*!
 * \ingroup Discretization
 * \brief Identity operator
 */
struct Identity : public Utility::Tag<Identity> {};

/*!
 * \ingroup Discretization
 * \brief Grad operator
 */
struct Gradient : public Utility::Tag<Gradient> {};
}

namespace Dumux::CVFE {

//! \brief Wrap as gradient expression.
template<class F>
Detail::Expr<F, Dumux::EvaluationOperators::Gradient> grad(const F& f) { return {f}; }

//! \brief Wrap as identity expression.
template<class F>
Detail::Expr<F, Dumux::EvaluationOperators::Identity> id(const F& f) { return {f}; }

//! \brief Create an identity expression from a function.
template<class F>
Detail::Expr<F, Dumux::EvaluationOperators::Identity> expr(const F& f) { return id(f); }

//! \brief Forward an already wrapped expression unchanged.
template<class F, class Op>
Detail::Expr<F, Op> expr(const Detail::Expr<F, Op>& e) { return e; }

//! \brief Evaluate a function at an interpolation point.
template<class ElementDiscretization, class ElementVariables, class IpData, class F,
         std::enable_if_t<!Detail::IsExpr<std::remove_cvref_t<F>>::value, int> = 0>
auto eval(const ElementDiscretization& elementDisc,
          const ElementVariables& elemVars,
          const IpData& ipData,
          const F& f)
{
    return eval(elementDisc, elemVars, ipData, expr(f));
}

/*! \brief Evaluate the identity expression at an interpolation point.
 *  Interpolates by using shape values of a local basis from `cache(elemVars, ipData)`.
 */
template<class ElementDiscretization, class ElementVariables, class IpData, class F>
auto eval(const ElementDiscretization& elementDisc,
          const ElementVariables& elemVars,
          const IpData& ipData,
          Detail::Expr<F,Dumux::EvaluationOperators::Identity> e)
{
    using Variables = typename ElementVariables::Variables;
    static_assert(std::is_invocable_v<const F&, const Variables&>,
                    "Function must be callable with ElementVariables::Variables");

    using FReturnType = std::remove_cvref_t<std::invoke_result_t<const F&, const Variables&>>;

    FReturnType result(0.0);
    const auto& ipCache = cache(elemVars, ipData);
    const auto& shapeValues = ipCache.shapeValues();
    for (const auto& localDof : localDofs(elementDisc))
        result += shapeValues[localDof.index()][0] * e.f(elemVars[localDof]);

    return result;
}

/*! \brief Evaluate the gradient expression at an interpolation point.
 *  Interpolates by using shape function gradients of a local basis from `cache(elemVars, ipData)`.
 *  Returns a component-wise gradient matrix.
 */
template<class ElementDiscretization, class ElementVariables, class IpData, class F>
auto eval(const ElementDiscretization& elementDisc,
          const ElementVariables& elemVars,
          const IpData& ipData,
          Detail::Expr<F,Dumux::EvaluationOperators::Gradient> e)
{
    using Variables = typename ElementVariables::Variables;
    static_assert(std::is_invocable_v<const F&, const Variables&>,
                    "Function must be callable with ElementVariables::Variables");

    using FReturnType = std::remove_cvref_t<std::invoke_result_t<const F&, const Variables&>>;
    using GlobalPosition = typename ElementDiscretization::Element::Geometry::GlobalCoordinate;
    using Scalar = typename GlobalPosition::value_type;
    using Gradients = Dune::FieldMatrix<Scalar, Detail::numComponents<FReturnType>(), GlobalPosition::dimension>;

    Gradients gradients(0.0);
    const auto& ipCache = cache(elemVars, ipData);
    for (const auto& localDof : localDofs(elementDisc))
    {
        const auto& vals = e.f(elemVars[localDof]);
        if constexpr (!Dune::IsNumber<FReturnType>::value)
            for (int comp = 0; comp < Detail::numComponents<FReturnType>(); ++comp)
                gradients[comp].axpy(vals[comp], ipCache.gradN(localDof.index()));
        else
            gradients[0].axpy(vals, ipCache.gradN(localDof.index()));
    }

    return gradients;
}

//! \brief Evaluate multiple expressions.
template<class ElementDiscretization, class ElementVariables, class IpData, class... Args>
requires (sizeof...(Args) > 1)
auto eval(const ElementDiscretization& elementDisc,
          const ElementVariables& elemVars,
          const IpData& ipData,
          const Args&... args)
{
    return std::tuple{
        eval(elementDisc, elemVars, ipData, expr(args))...
    };
}

} // namespace Dumux::CVFE

#endif
