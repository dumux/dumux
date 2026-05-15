// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Interpolate a given function into the function space of the discretization
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_INTERPOLATE_HH
#define DUMUX_DISCRETIZATION_CVFE_INTERPOLATE_HH

#include <cassert>
#include <optional>
#include <ranges>
#include <type_traits>

#include <dune/grid/common/rangegenerators.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/common/tag.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/projection/l2_projection.hh>

namespace Dumux::CVFE::InterpolationPolicy {

/*!
 * \ingroup Discretization
 * \brief Interpolation policy that evaluates a function at the dof positions.
 */
struct DofPositionEvaluation : public Utility::Tag<DofPositionEvaluation> {};

/*!
 * \ingroup Discretization
 * \brief Interpolation policy that uses an L2 projection.
 */
struct L2Projection : public Utility::Tag<L2Projection>
{
    std::size_t maxIterations = 100;
    double residualReduction = 1e-13;
    int verbosity = 0;
};


} // end namespace Dumux::CVFE::InterpolationPolicy

namespace Dumux::CVFE {

namespace Detail {

template<class ElementDiscretization, class Function, class ValueType>
concept InterpolatableAtElemDiscIpData = requires(
    const ElementDiscretization& eDisc,
    const std::ranges::range_value_t<decltype(localDofs(eDisc))>& ld,
    Function& f,
    ValueType& v)
{
    { ipData(eDisc, ld) } -> Concept::LocalDofIpData;
    v = f(eDisc, ipData(eDisc, ld));
};

template<class ElementDiscretization, class Function, class ValueType>
concept InterpolatableAtGlobalPos = requires(
    typename ElementDiscretization::Element::Geometry::GlobalCoordinate globalPos,
    Function& f,
    ValueType& v)
{ v = f(globalPos); };

template<class ElementDiscretization, class Function, class ValueType>
concept InterpolatableAtElemDiscGlobalPos = requires(
    ElementDiscretization& eDisc,
    typename ElementDiscretization::Element::Geometry::GlobalCoordinate globalPos,
    Function& f,
    ValueType& v)
{ v = f(eDisc, globalPos); };

/*!
 * \brief Wraps a function callable as f(globalPos) or f(elemDisc, globalPos)
 *        into a local-function interface: bind(element) and operator()(localCoord).
 *
 * This allows both signatures to be used uniformly in e.g. L2 projection,
 * where the function is evaluated at quadrature points in local coordinates.
 */
template<class GridDiscretization, class Function>
class LocalFunctionAdapter
{
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using Element = typename ElementDiscretization::Element;
    using LocalCoord = typename Element::Geometry::LocalCoordinate;
    using GlobalPos = typename Element::Geometry::GlobalCoordinate;

public:
    LocalFunctionAdapter(const GridDiscretization& gridDisc, Function f)
    : function_(std::move(f))
    , elemDisc_(localView(gridDisc))
    {}

    void bind(const Element& element)
    {
        element_.emplace(element);
        elemDisc_.bind(element);
    }

    auto operator()(const LocalCoord& localCoord)
    {
        const GlobalPos globalPos = elemDisc_.elementGeometry().global(localCoord);
        if constexpr (requires(ElementDiscretization& elemDisc , GlobalPos g, std::decay_t<Function>& f){ f(elemDisc, g); })
            return function_(elemDisc_, globalPos);
        else
            return function_(globalPos);
    }

private:
    Function function_;
    ElementDiscretization elemDisc_;
    std::optional<Element> element_;
};

template<class GridDiscretization, class Function>
auto makeLocalFunctionAdapter(const GridDiscretization& gridDisc, Function&& f)
{
    return LocalFunctionAdapter<GridDiscretization, std::decay_t<Function>>(
        gridDisc, std::forward<Function>(f));
}

} // end namespace Detail

/*!
 * \ingroup Discretization
 * \brief Interpolate a function into a coefficient vector by evaluating it at the dof positions.
 *
 * Accepts functions callable as function(elemDisc, ipData), function(elemDisc, globalPos), or function(globalPos).
 */
template<class GridDiscretization, class CoefficientVector, class Function>
void interpolate(const GridDiscretization& gridDisc,
                 CoefficientVector& coeffs,
                 Function&& function,
                 InterpolationPolicy::DofPositionEvaluation = {})
requires (Detail::InterpolatableAtElemDiscIpData<typename GridDiscretization::LocalView, Function, typename CoefficientVector::value_type>
       || Detail::InterpolatableAtElemDiscGlobalPos<typename GridDiscretization::LocalView, Function, typename CoefficientVector::value_type>
       || Detail::InterpolatableAtGlobalPos<typename GridDiscretization::LocalView, Function, typename CoefficientVector::value_type>)
{
    using LocalView = typename GridDiscretization::LocalView;
    using ValueType = typename CoefficientVector::value_type;
    auto elemDisc = localView(gridDisc);
    for (const auto& element : elements(gridDisc.gridView()))
    {
        elemDisc.bind(element);
        for (const auto& localDof : localDofs(elemDisc))
        {
            if constexpr (Detail::InterpolatableAtElemDiscIpData<LocalView, Function, ValueType>)
                coeffs[localDof.dofIndex()] = function(elemDisc, ipData(elemDisc, localDof));
            else if constexpr (Detail::InterpolatableAtElemDiscGlobalPos<LocalView, Function, ValueType>)
                coeffs[localDof.dofIndex()] = function(elemDisc, ipData(elemDisc, localDof).global());
            else
                coeffs[localDof.dofIndex()] = function(ipData(elemDisc, localDof).global());
        }
    }
}

/*!
 * \ingroup Discretization
 * \brief Interpolate a function into a coefficient vector using an L2 projection.
 *
 * The function must be callable as function(globalPos) or function(elemDisc, globalPos),
 * and is wrapped into a local-function interface for use in L2 projection.
 */
template<class GridDiscretization, class CoefficientVector, class Function>
void interpolate(const GridDiscretization& gridDisc,
                 CoefficientVector& coeffs,
                 Function&& function,
                 InterpolationPolicy::L2Projection policy)
requires (Detail::InterpolatableAtGlobalPos<typename GridDiscretization::LocalView, Function, typename CoefficientVector::value_type>
       || Detail::InterpolatableAtElemDiscGlobalPos<typename GridDiscretization::LocalView, Function, typename CoefficientVector::value_type>)
{
    using Basis = FEBasisFromCVFEGridDiscretization<GridDiscretization>;
    const Basis basis(gridDisc);
    using L2Proj = Dumux::L2Projection<Basis>;

    typename L2Proj::Params params;
    params.maxIterations = policy.maxIterations;
    params.residualReduction = policy.residualReduction;
    params.verbosity = policy.verbosity;

    auto localFunction = Detail::makeLocalFunctionAdapter(gridDisc, std::forward<Function>(function));
    const L2Proj projection(basis);
    const auto result = projection.project(localFunction, params);
    assert(coeffs.size() == result.size());
    using ResultValue = typename decltype(result)::value_type; // FieldVector<Scalar, numEq>
    constexpr int numEq = ResultValue::dimension;
    for (std::size_t i = 0; i < result.size(); ++i)
        if constexpr (numEq == 1)
            coeffs[i] = result[i][0];
        else
            coeffs[i] = result[i];
}

} // end namespace Dumux::CVFE

#endif
