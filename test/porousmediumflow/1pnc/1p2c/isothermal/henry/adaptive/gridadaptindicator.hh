// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A solute-concentration indicator for adaptive refinement of the Henry benchmark.
 *
 * Same role as dumux/porousmediumflow/2p/gridadaptindicator.hh (element-wise
 * neighbor-delta of a scalar field, global-extrema-normalized refine/coarsen bounds,
 * 2:1-balance neighbor walk), keyed on the solute mass fraction X^solute
 * (FluidSystem::soluteIdx, the same primary-variable index problem.hh already indexes
 * PrimaryVariables/NumEqVector with directly) instead of saturation. Unlike the 2p
 * indicator, there is no phase mass/state to worry about (OnePNC's X^solute is a plain
 * primary variable that always exists), so no special-casing is needed for that.
 */
#ifndef DUMUX_HENRY_ADAPTIVE_GRIDADAPTINDICATOR_HH
#define DUMUX_HENRY_ADAPTIVE_GRIDADAPTINDICATOR_HH

#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>

namespace Dumux {

/*!
 * \brief Concentration-jump indicator for grid adaptation of the Henry benchmark.
 * \tparam TypeTag the problem TypeTag (HenryFahsBenchmarkTest or HenryFahsCase2BenchmarkTest)
 *
 * \note The primary-variable index of the solute mass fraction is taken as a constructor
 *       argument (FluidSystem::soluteIdx, see problem.hh, which already indexes
 *       PrimaryVariables/NumEqVector with it directly) rather than derived from Indices.
 */
template<class TypeTag>
class HenryGridAdaptIndicator
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

public:
    /*!
     * \brief The constructor
     * \param gridGeometry the finite volume grid geometry
     * \param soluteEqIdx primary-variable index of the solute mass fraction (FluidSystem::soluteIdx)
     * \param paramGroup the parameter group to read Adaptive.MinLevel/MaxLevel from
     */
    HenryGridAdaptIndicator(std::shared_ptr<const GridGeometry> gridGeometry,
                            int soluteEqIdx,
                            const std::string& paramGroup = "")
    : gridGeometry_(gridGeometry)
    , soluteIdx_(soluteEqIdx)
    , refineBound_(std::numeric_limits<Scalar>::max())
    , coarsenBound_(std::numeric_limits<Scalar>::lowest())
    , maxSoluteDelta_(gridGeometry_->gridView().size(0), 0.0)
    , minLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MinLevel", 0))
    , maxLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MaxLevel", 0))
    {}

    void setMinLevel(std::size_t minLevel) { minLevel_ = minLevel; }
    void setMaxLevel(std::size_t maxLevel) { maxLevel_ = maxLevel; }
    void setLevels(std::size_t minLevel, std::size_t maxLevel)
    { minLevel_ = minLevel; maxLevel_ = maxLevel; }

    /*!
     * \brief Calculate the refine/coarsen indicator for each element, based on the maximum
     *        jump in X^solute across element interfaces, normalized by the global range of
     *        X^solute.
     *
     * \param sol the current solution
     * \param refineTol elements whose max neighbor-delta exceeds refineTol*globalRange(X^solute)
     *                  are marked for refinement (subject to Adaptive.MaxLevel)
     * \param coarsenTol elements whose max neighbor-delta is below coarsenTol*globalRange(X^solute)
     *                   are marked for coarsening (subject to Adaptive.MinLevel)
     */
    void calculate(const SolutionVector& sol,
                   Scalar refineTol = 0.05,
                   Scalar coarsenTol = 0.001)
    {
        refineBound_ = std::numeric_limits<Scalar>::max();
        coarsenBound_ = std::numeric_limits<Scalar>::lowest();
        maxSoluteDelta_.assign(gridGeometry_->gridView().size(0), 0.0);

        if (minLevel_ > maxLevel_)
            DUNE_THROW(Dune::InvalidStateException, "Adaptive.MinLevel must not exceed Adaptive.MaxLevel");
        else if (minLevel_ == maxLevel_)
            return; // adaptivity disabled (default: MinLevel = MaxLevel = 0)

        if (coarsenTol > refineTol)
            DUNE_THROW(Dune::InvalidStateException, "Refine tolerance must be higher than coarsen tolerance");

        const auto& gridView = gridGeometry_->gridView();

        Scalar globalMax = std::numeric_limits<Scalar>::lowest();
        Scalar globalMin = std::numeric_limits<Scalar>::max();

        for (const auto& element : elements(gridView))
        {
            const auto globalIdxI = gridGeometry_->elementMapper().index(element);

            const auto geometry = element.geometry();
            const auto elemSol = elementSolution(element, sol, *gridGeometry_);
            const Scalar cI = evalSolution(element, geometry, *gridGeometry_, elemSol, geometry.center())[soluteIdx_];

            using std::min; using std::max;
            globalMin = min(cI, globalMin);
            globalMax = max(cI, globalMax);

            for (const auto& intersection : intersections(gridView, element))
            {
                if (!intersection.neighbor())
                    continue;

                const auto outside = intersection.outside();
                const auto globalIdxJ = gridGeometry_->elementMapper().index(outside);

                // visit each interior facet only once: elementMapper() gives every leaf
                // element a unique index regardless of level, so comparing indices alone
                // is enough
                if (globalIdxI < globalIdxJ)
                {
                    const auto outsideGeometry = outside.geometry();
                    const auto elemSolJ = elementSolution(outside, sol, *gridGeometry_);
                    const Scalar cJ = evalSolution(outside, outsideGeometry, *gridGeometry_, elemSolJ, outsideGeometry.center())[soluteIdx_];

                    using std::abs;
                    const Scalar localDelta = abs(cI - cJ);
                    maxSoluteDelta_[globalIdxI] = max(maxSoluteDelta_[globalIdxI], localDelta);
                    maxSoluteDelta_[globalIdxJ] = max(maxSoluteDelta_[globalIdxJ], localDelta);
                }
            }
        }

        const Scalar globalDelta = globalMax - globalMin;
        refineBound_ = refineTol * globalDelta;
        coarsenBound_ = coarsenTol * globalDelta;

        // 2:1-balance: ensure any element whose refinement would create a level jump larger
        // than one between neighbors is refined too.
        for (const auto& element : elements(gridView))
            if (this->operator()(element) > 0)
                checkNeighborsRefine_(element);
    }

    /*!
     * \brief function call operator
     * \return  1 if the element should be refined, -1 if it should be coarsened, 0 otherwise
     */
    int operator() (const Element& element) const
    {
        const auto idx = gridGeometry_->elementMapper().index(element);
        if (element.hasFather() && maxSoluteDelta_[idx] < coarsenBound_)
            return -1;
        else if (element.level() < maxLevel_ && maxSoluteDelta_[idx] > refineBound_)
            return 1;
        else
            return 0;
    }

    //! Per-element indicator value computed by the last calculate() call.
    const std::vector<Scalar>& values() const { return maxSoluteDelta_; }
    Scalar refineBound() const { return refineBound_; }
    Scalar coarsenBound() const { return coarsenBound_; }

private:

    // element is about to refine by one level. Force any neighbor that's already
    // more than one level coarser (will be two levels after refinement) to refine too, then recurse into that neighbor
    // since forcing it can create the same violation one ring further out.
    bool checkNeighborsRefine_(const Element& element, std::size_t level = 1)
    {
        for (const auto& intersection : intersections(gridGeometry_->gridView(), element))
        {
            if (!intersection.neighbor())
                continue;

            const auto outside = intersection.outside();
            if (outside.level() < maxLevel_ && outside.level() < element.level())
            {
                maxSoluteDelta_[gridGeometry_->elementMapper().index(outside)] = std::numeric_limits<Scalar>::max();
                if (level < maxLevel_)
                    checkNeighborsRefine_(outside, level + 1);
            }
        }
        return true;
    }

    std::shared_ptr<const GridGeometry> gridGeometry_;
    int soluteIdx_;
    Scalar refineBound_;
    Scalar coarsenBound_;
    std::vector<Scalar> maxSoluteDelta_;
    std::size_t minLevel_;
    std::size_t maxLevel_;
};

} // end namespace Dumux

#endif
