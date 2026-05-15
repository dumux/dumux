// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief L2-projections of analytic functions into a given function space
 */
#ifndef DUMUX_DISCRETIZATION_L2_PROJECTION_HH
#define DUMUX_DISCRETIZATION_L2_PROJECTION_HH

#include <array>
#include <optional>
#include <type_traits>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/typetraits.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#ifdef HAVE_DUNE_FUNCTIONS
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#endif
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/parallel/parallel_for.hh>

namespace Dumux {

namespace Detail {

/*!
 * \brief Create a local function from a given function.
 */
template<class Function, class GridView>
auto makeLocalFunction(Function&& f, const GridView& gridView)
{
    using Element = typename GridView::template Codim<0>::Entity;
    using LocalCoord = typename Element::Geometry::LocalCoordinate;
    constexpr bool hasLocalInterface =
        requires(std::decay_t<Function>& lf, const Element& e, const LocalCoord& x)
        {
            lf.bind(e);
            lf(x);
        };

    // If f already provides the local-function interface, use it directly.
    // This must be checked first to avoid incorrectly re-wrapping it via dune-functions.
    if constexpr (hasLocalInterface)
        return std::forward<Function>(f);
#ifdef HAVE_DUNE_FUNCTIONS
    else
        return localFunction(Dune::Functions::makeGridViewFunction(std::forward<Function>(f), gridView));
#else
    else
        DUNE_THROW(Dune::InvalidStateException, "Function must provide bind(element) and operator()(localCoord), "
                                                "or dune-functions must be available to wrap a global f(globalPos) callable.");
#endif
}

} // end namespace Detail

/*!
 * \ingroup Discretization
 * \brief Wraps a CVFE grid discretization to expose the FE basis interface expected by L2Projection.
 *
 * Wraps a CVFE GridDiscretization to expose the minimal interface expected by L2Projection:
 * - size(), gridView(), localView()
 * - localView provides bind(element), tree().finiteElement(), index(i)
 */
template<class GridDiscretization>
class FEBasisFromCVFEGridDiscretization
{
    using GV = typename GridDiscretization::GridView;
    using Element = typename GV::template Codim<0>::Entity;
    using FE = typename GridDiscretization::FeCache::FiniteElementType;

public:
    using GridView = GV;

    struct LocalTree
    {
        using FiniteElement = FE;
        const FE* fe_ = nullptr;
        const FiniteElement& finiteElement() const { return *fe_; }
    };

    struct LocalView
    {
        using Tree = LocalTree;

        explicit LocalView(const GridDiscretization& gg) : gg_(gg) {}

        void bind(const Element& element)
        {
            element_ = element;
            tree_.fe_ = &gg_.feCache().get(element.type());
        }

        const Tree& tree() const { return tree_; }

        std::size_t index(std::size_t index) const
        {
            const auto& localKey = tree_.fe_->localCoefficients().localKey(index);
            // TODO: Currently we assume that this is the default dof mapping when having multiple dofs per sub-entity
            // meaning that they are localKey.index() larger than 0
            return gg_.dofMapper().subIndex(*element_, localKey.subEntity(), localKey.codim()) + localKey.index();
        }

    private:
        const GridDiscretization& gg_;
        std::optional<Element> element_;
        Tree tree_;
    };

    explicit FEBasisFromCVFEGridDiscretization(const GridDiscretization& gg) : gg_(gg) {}

    std::size_t size() const { return gg_.numDofs(); }
    const GridView& gridView() const { return gg_.gridView(); }
    LocalView localView() const { return LocalView(gg_); }

private:
    const GridDiscretization& gg_;
};

template <class FEBasis>
class L2Projection
{
    static constexpr int dim = FEBasis::GridView::dimension;
    using FiniteElement = typename FEBasis::LocalView::Tree::FiniteElement;
    using Scalar = typename FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType;
    using ShapeValue = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
    static_assert(ShapeValue::dimension == 1, "Only scalar-valued shape functions are supported for L2 projection.");
    using Matrix = Dune::BCRSMatrix<Scalar>;

public:
    template<int numEq = 1>
    using CoefficientVector = Dune::BlockVector<Dune::FieldVector<Scalar, numEq>>;

    //! Parameters that can be passed to project()
    struct Params
    {
        std::size_t maxIterations{100};
        Scalar residualReduction{1e-13};
        int verbosity{0};
    };

    L2Projection(const FEBasis& feBasis)
    : feBasis_(feBasis)
    , solver_()
    {
        solver_.setMatrix(std::make_shared<Matrix>(createMassMatrix_(feBasis)));
    }

    template <class Function>
    auto project(Function&& function, const Params& params = Params{}) const
    {
        auto localFunc = Detail::makeLocalFunction(std::forward<Function>(function), feBasis_.gridView());
        using LocalPosition = typename FEBasis::GridView::template Codim<0>::Entity::Geometry::LocalCoordinate;
        using ReturnType = std::invoke_result_t<Function, LocalPosition>;
        constexpr int numEq = []() {
            if constexpr (Dune::IsNumber<ReturnType>::value) return 1;
            else return ReturnType::dimension;
        }();
        using CoeffVec = CoefficientVector<numEq>;

        const auto numDofs = feBasis_.size();
        // assemble one RHS per equation component
        std::array<Dune::BlockVector<Scalar>, numEq> rhs;
        for (auto& r : rhs) { r.resize(numDofs); r = 0.0; }

        // assemble right hand side
        auto localView = feBasis_.localView();
        for (const auto& element : elements(feBasis_.gridView()))
        {
            localView.bind(element);
            localFunc.bind(element);

            const auto& localFiniteElement = localView.tree().finiteElement();
            const int order = dim*localFiniteElement.localBasis().order();
            const auto& quad = Dune::QuadratureRules<Scalar, dim>::rule(element.type(), order);
            const auto geometry = element.geometry();

            for (auto&& qp : quad)
            {
                const auto weight = qp.weight();
                const auto ie = geometry.integrationElement(qp.position());

                std::vector<ShapeValue> shapeValues;
                localFiniteElement.localBasis().evaluateFunction(qp.position(), shapeValues);
                const auto functionValue = localFunc(qp.position());

                for (int i = 0; i < localFiniteElement.localBasis().size(); ++i)
                {
                    const auto globalI = localView.index(i);
                    const auto w = ie*weight*shapeValues[i][0];
                    if constexpr (numEq == 1)
                        rhs[0][globalI] += w*functionValue;
                    else
                        for (int compIdx = 0; compIdx < numEq; compIdx++)
                            rhs[compIdx][globalI] += w*functionValue[compIdx];
                }
            }
        }

        // solve numEq independent systems in parallel (same matrix, different RHS)
        CoeffVec coeffs(numDofs);
        Dumux::parallelFor(numEq, [&](std::size_t compIdx)
        {
            // each thread needs its own solver copy (not thread-safe to share)
            auto solver = solver_;
            Dune::ParameterTree solverParams;
            solverParams["maxit"] = std::to_string(params.maxIterations);
            solverParams["reduction"] = std::to_string(params.residualReduction);
            solverParams["verbose"] = std::to_string(params.verbosity);
            solver.setParams(solverParams);

            Dune::BlockVector<Scalar> sol(numDofs); sol = 0.0;
            solver.solve(sol, rhs[compIdx]);

            for (std::size_t i = 0; i < numDofs; ++i)
                coeffs[i][compIdx] = sol[i];
        });

        return coeffs;
    }

private:
    Matrix createMassMatrix_(const FEBasis& feBasis) const
    {
        Matrix massMatrix;

        auto pattern = getFEJacobianPattern(feBasis);
        pattern.exportIdx(massMatrix);
        massMatrix = 0.0;

        auto localView = feBasis.localView();
        for (const auto& element : elements(feBasis.gridView()))
        {
            localView.bind(element);

            const auto& localFiniteElement = localView.tree().finiteElement();
            const int order = 2*dim*localFiniteElement.localBasis().order();
            const auto& quad = Dune::QuadratureRules<Scalar, dim>::rule(element.type(), order);
            const auto geometry = element.geometry();

            for (auto&& qp : quad)
            {
                const auto weight = qp.weight();
                const auto ie = geometry.integrationElement(qp.position());

                std::vector<ShapeValue> shapeValues;
                localFiniteElement.localBasis().evaluateFunction(qp.position(), shapeValues);

                for (int i = 0; i < localFiniteElement.localBasis().size(); ++i)
                {
                    const auto globalI = localView.index(i);
                    massMatrix[globalI][globalI] += ie*weight*shapeValues[i]*shapeValues[i];

                    for (int j = i+1; j < localFiniteElement.localBasis().size(); ++j)
                    {
                        const auto globalJ = localView.index(j);
                        const auto value = ie*weight*shapeValues[i]*shapeValues[j];
                        massMatrix[globalI][globalJ] += value;
                        massMatrix[globalJ][globalI] += value;
                    }
                }
            }
        }

        return massMatrix;
    }

    const FEBasis& feBasis_;
    SSORCGIstlSolver<
        SeqLinearSolverTraits, LinearAlgebraTraits<Matrix, Dune::BlockVector<Scalar>>
    > solver_;
};

} // end namespace Dumux

#endif
