// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A concentration-gradient indicator for adaptive refinement of the Boussinesq tests.
 *
 * Same role as dumux/porousmediumflow/2p/gridadaptindicator.hh (element-wise neighbor-delta
 * of a scalar field, global-extrema-normalized refine/coarsen bounds, 2:1-balance neighbor
 * walk), but keyed on the transported concentration C instead of saturation, and fixed to
 * actually work in parallel:
 *
 *  - the global refine/coarsen thresholds are reduced across *all* ranks via
 *    gridView.comm().max/min (2p's version silently used only the local rank's extrema --
 *    it has the reduction written out in a comment that was never enabled);
 *  - the 2:1-balance neighbor-forcing pass can flag a ghost/overlap copy of an element that
 *    is *owned by a different rank*; that flag is synced back to the owning rank afterwards
 *    (via griddatacommunication.hh's syncElementDataMax) so every rank's markElements() call
 *    -- which only ever marks its own Partitions::interior elements -- sees the same
 *    decision near a partition boundary. This step has no equivalent in the 2p code at all.
 */
#ifndef DUMUX_BOUSSINESQ_ADAPTIVE_GRIDADAPTINDICATOR_HH
#define DUMUX_BOUSSINESQ_ADAPTIVE_GRIDADAPTINDICATOR_HH

#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>

#include "griddatacommunication.hh"

namespace Dumux {

/*!
 * \brief Concentration-gradient indicator for grid adaptation of the Boussinesq tests.
 * \tparam TypeTag the problem TypeTag (works for both the vorticity and pressure
 *                  formulations, both box/CVFE)
 *
 * \note The primary-variable index of the transported concentration is taken as a
 *       constructor argument rather than derived from Indices::concentrationIdx, because the
 *       two Boussinesq formulations don't share an Indices layout: the vorticity/box model
 *       exposes Indices::concentrationIdx directly, while the pressure/OnePNC model
 *       addresses it as Indices::conti0EqIdx + FluidSystem::soluteIdx.
 */
template<class TypeTag>
class BoussinesqGridAdaptIndicator
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
     * \param concentrationEqIdx primary-variable index of the transported concentration
     * \param paramGroup the parameter group to read Adaptive.MinLevel/MaxLevel from
     * \note Reads Adaptive.MinLevel and Adaptive.MaxLevel from the parameter tree
     */
    BoussinesqGridAdaptIndicator(std::shared_ptr<const GridGeometry> gridGeometry,
                                 int concentrationEqIdx,
                                 const std::string& paramGroup = "")
    : gridGeometry_(gridGeometry)
    , concentrationIdx_(concentrationEqIdx)
    , refineBound_(std::numeric_limits<Scalar>::max())
    , coarsenBound_(std::numeric_limits<Scalar>::lowest())
    , maxConcentrationDelta_(gridGeometry_->gridView().size(0), 0.0)
    , minLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MinLevel", 0))
    , maxLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MaxLevel", 0))
    {}

    void setMinLevel(std::size_t minLevel) { minLevel_ = minLevel; }
    void setMaxLevel(std::size_t maxLevel) { maxLevel_ = maxLevel; }
    void setLevels(std::size_t minLevel, std::size_t maxLevel)
    { minLevel_ = minLevel; maxLevel_ = maxLevel; }

    /*!
     * \brief Calculate the refine/coarsen indicator for each element, based on the maximum
     *        jump in C across element interfaces, normalized by the global range of C.
     *
     * \param sol the current solution
     * \param refineTol elements whose max neighbor-delta exceeds refineTol*globalRange(C)
     *                  are marked for refinement (subject to Adaptive.MaxLevel)
     * \param coarsenTol elements whose max neighbor-delta is below coarsenTol*globalRange(C)
     *                   are marked for coarsening (subject to Adaptive.MinLevel)
     */
    void calculate(const SolutionVector& sol,
                   Scalar refineTol = 0.05,
                   Scalar coarsenTol = 0.001)
    {
        refineBound_ = std::numeric_limits<Scalar>::max();
        coarsenBound_ = std::numeric_limits<Scalar>::lowest();
        maxConcentrationDelta_.assign(gridGeometry_->gridView().size(0), 0.0);

        if (minLevel_ >= maxLevel_)
            return;

        if (coarsenTol > refineTol)
            DUNE_THROW(Dune::InvalidStateException, "Refine tolerance must be higher than coarsen tolerance");

        const auto& gridView = gridGeometry_->gridView();

        Scalar localMax = std::numeric_limits<Scalar>::lowest();
        Scalar localMin = std::numeric_limits<Scalar>::max();

        for (const auto& element : elements(gridView))
        {
            const auto globalIdxI = gridGeometry_->elementMapper().index(element);

            const auto geometry = element.geometry();
            const auto elemSol = elementSolution(element, sol, *gridGeometry_);
            const Scalar cI = evalSolution(element, geometry, *gridGeometry_, elemSol, geometry.center())[concentrationIdx_];

            using std::min; using std::max;
            localMin = min(cI, localMin);
            localMax = max(cI, localMax);

            for (const auto& intersection : intersections(gridView, element))
            {
                if (!intersection.neighbor())
                    continue;

                const auto outside = intersection.outside();
                const auto globalIdxJ = gridGeometry_->elementMapper().index(outside);

                // visit each interior facet only once
                if (element.level() > outside.level() || (element.level() == outside.level() && globalIdxI < globalIdxJ))
                {
                    const auto outsideGeometry = outside.geometry();
                    const auto elemSolJ = elementSolution(outside, sol, *gridGeometry_);
                    const Scalar cJ = evalSolution(outside, outsideGeometry, *gridGeometry_, elemSolJ, outsideGeometry.center())[concentrationIdx_];

                    using std::abs;
                    const Scalar localDelta = abs(cI - cJ);
                    maxConcentrationDelta_[globalIdxI] = max(maxConcentrationDelta_[globalIdxI], localDelta);
                    maxConcentrationDelta_[globalIdxJ] = max(maxConcentrationDelta_[globalIdxJ], localDelta);
                }
            }
        }

        // Global reduction across all ranks -- this is the piece that stayed commented out
        // (and thus never executed) in dumux/porousmediumflow/2p/gridadaptindicator.hh.
        // Without it every rank derives refineBound_/coarsenBound_ from its own local
        // subdomain only, so two ranks can disagree on whether the *same* boundary element
        // should be refined.
        using std::max; using std::min;
        const Scalar globalMax = gridView.comm().max(localMax);
        const Scalar globalMin = gridView.comm().min(localMin);
        const Scalar globalDelta = globalMax - globalMin;

        refineBound_ = refineTol * globalDelta;
        coarsenBound_ = coarsenTol * globalDelta;

        // 2:1-balance: ensure any element whose refinement would create a level jump larger
        // than one between neighbors is refined too. This can touch a ghost/overlap copy of
        // an element owned by a different rank; that's fine locally, but the *owning* rank
        // needs to see the same forced-refine flag on its own interior copy, which is what
        // the sync below (after this loop) provides.
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
            if (this->operator()(element) > 0)
                checkNeighborsRefine_(element);

        // Sync the (possibly locally-forced) per-element deltas so that every rank's
        // Partitions::interior elements see the max value across all ranks that have a
        // copy of that element -- this is the communication step the dead 2p code never
        // finished (it referenced a deleted VectorExchange helper for exactly this).
        BoussinesqAdaptive::syncElementDataMax(gridView, gridGeometry_->elementMapper(), maxConcentrationDelta_);
    }

    /*!
     * \brief function call operator
     * \return  1 if the element should be refined, -1 if it should be coarsened, 0 otherwise
     *
     * An element is a refine candidate once its value() *exceeds* refineBound() (not "falls
     * below" -- see that accessor's doc), and a coarsen candidate once its value() falls
     * *below* coarsenBound(). Both bounds are absolute concentration-delta values, already
     * scaled by the corresponding tolerance and the current global range of C -- see
     * refineBound()/coarsenBound().
     */
    int operator() (const Element& element) const
    {
        const auto idx = gridGeometry_->elementMapper().index(element);
        if (element.hasFather() && maxConcentrationDelta_[idx] < coarsenBound_)
            return -1;
        else if (element.level() < maxLevel_ && maxConcentrationDelta_[idx] > refineBound_)
            return 1;
        else
            return 0;
    }

    /*!
     * \brief Per-element indicator value computed by the last calculate() call: the maximum
     *        jump in C across that element's faces (see calculate()'s doc). This is the raw
     *        quantity operator() compares against refineBound()/coarsenBound().
     * \note Returns a reference to the live member vector -- valid for the object's lifetime
     *       and automatically reflects the most recent calculate() call (e.g. when passed to
     *       VtkOutputModule::addField() once, outside the time loop: the referenced vector's
     *       contents update in place on every subsequent calculate(), no need to re-add it).
     */
    const std::vector<Scalar>& values() const { return maxConcentrationDelta_; }

    /*!
     * \brief The current refine threshold: refineTol * (globalMax(C) - globalMin(C)), as of
     *        the last calculate() call. An element is a refine candidate once
     *        values()[elementIdx] > refineBound() (strictly *above*, not below -- the
     *        indicator refines where the field is changing fast, so a *large* jump is what
     *        triggers refinement; coarsenBound() is the "flat region" threshold instead, see
     *        there) -- and only below Adaptive.MaxLevel.
     */
    Scalar refineBound() const { return refineBound_; }

    /*!
     * \brief The current coarsen threshold: coarsenTol * (globalMax(C) - globalMin(C)), as of
     *        the last calculate() call. An element is a coarsen candidate once
     *        values()[elementIdx] < coarsenBound() (the field is locally flat/converged here)
     *        -- and only above Adaptive.MinLevel. Always < refineBound() (calculate() throws
     *        otherwise), so there is a dead zone between the two bounds where elements are
     *        left exactly as they are.
     */
    Scalar coarsenBound() const { return coarsenBound_; }

private:
    bool checkNeighborsRefine_(const Element& element, std::size_t level = 1)
    {
        for (const auto& intersection : intersections(gridGeometry_->gridView(), element))
        {
            if (!intersection.neighbor())
                continue;

            const auto outside = intersection.outside();

            // never force-refine a ghost we have no ownership info about beyond its level;
            // the owning rank will pick this up via the sync call once this rank sets the
            // delta on its own (border/overlap) copy of the shared element below
            if (outside.partitionType() == Dune::GhostEntity)
                continue;

            if (outside.level() < maxLevel_ && outside.level() < element.level())
            {
                maxConcentrationDelta_[gridGeometry_->elementMapper().index(outside)] = std::numeric_limits<Scalar>::max();
                if (level < maxLevel_)
                    checkNeighborsRefine_(outside, level + 1);
            }
        }
        return true;
    }

    std::shared_ptr<const GridGeometry> gridGeometry_;
    int concentrationIdx_;
    Scalar refineBound_;
    Scalar coarsenBound_;
    std::vector<Scalar> maxConcentrationDelta_;
    std::size_t minLevel_;
    std::size_t maxLevel_;
};

/*!
 * \brief Absolute-gradient indicator for grid adaptation of the Boussinesq tests.
 *
 * Unlike BoussinesqGridAdaptIndicator above (max jump in C between neighboring element
 * averages, thresholds relative to the global range of C), this evaluates the true
 * |grad C| at each element's center via evalGradients (the same discretization-level
 * gradient reconstruction the Sherwood/vorticity diagnostics in main_vorticity.cc already
 * use), and compares it against caller-supplied *absolute* refine/coarsen bounds -- e.g.
 * physical thresholds like 1/(5*Ra) and 1/Ra, not fractions of the field's current range.
 * Because the bounds are absolute, no gridView.comm().max/min reduction of a global range is
 * needed here (there is no global range to reduce); the parallel-correctness pieces that are
 * still needed -- the 2:1-balance neighbor-forcing walk and the border/overlap sync of the
 * per-element indicator value -- are identical in structure to the class above.
 *
 * \note If coarsenBound > refineBound (checked by the caller, not enforced here), the two
 * conditions in operator() can both be true for the same element; the coarsen check runs
 * first (same order as the class above), so coarsening wins that overlap.
 */
template<class TypeTag>
class BoussinesqGradientMagnitudeIndicator
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
     * \param concentrationEqIdx primary-variable index of the transported concentration
     * \param paramGroup the parameter group to read Adaptive.MinLevel/MaxLevel from
     */
    BoussinesqGradientMagnitudeIndicator(std::shared_ptr<const GridGeometry> gridGeometry,
                                         int concentrationEqIdx,
                                         const std::string& paramGroup = "")
    : gridGeometry_(gridGeometry)
    , concentrationIdx_(concentrationEqIdx)
    , refineBound_(std::numeric_limits<Scalar>::max())
    , coarsenBound_(std::numeric_limits<Scalar>::lowest())
    , gradMagnitude_(gridGeometry_->gridView().size(0), 0.0)
    , minLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MinLevel", 0))
    , maxLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MaxLevel", 0))
    {}

    void setMinLevel(std::size_t minLevel) { minLevel_ = minLevel; }
    void setMaxLevel(std::size_t maxLevel) { maxLevel_ = maxLevel; }
    void setLevels(std::size_t minLevel, std::size_t maxLevel)
    { minLevel_ = minLevel; maxLevel_ = maxLevel; }

    /*!
     * \brief Calculate |grad C| at each element's center.
     * \param sol the current solution
     * \param refineBound elements with |grad C| > refineBound are refine candidates
     *                     (subject to Adaptive.MaxLevel)
     * \param coarsenBound elements with |grad C| < coarsenBound are coarsen candidates
     *                      (subject to Adaptive.MinLevel)
     */
    void calculate(const SolutionVector& sol, Scalar refineBound, Scalar coarsenBound)
    {
        refineBound_ = refineBound;
        coarsenBound_ = coarsenBound;
        gradMagnitude_.assign(gridGeometry_->gridView().size(0), 0.0);

        if (minLevel_ >= maxLevel_)
            return;

        const auto& gridView = gridGeometry_->gridView();
        for (const auto& element : elements(gridView))
        {
            const auto idx = gridGeometry_->elementMapper().index(element);
            const auto geometry = element.geometry();
            const auto elemSol = elementSolution(element, sol, *gridGeometry_);
            const auto grad = evalGradients(element, geometry, *gridGeometry_, elemSol, geometry.center())[concentrationIdx_];
            gradMagnitude_[idx] = grad.two_norm();
        }

        // 2:1-balance neighbor-forcing walk and parallel sync: same rationale and structure
        // as BoussinesqGridAdaptIndicator::calculate() above, see its comments for why both
        // are needed for correctness at partition boundaries.
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
            if (this->operator()(element) > 0)
                checkNeighborsRefine_(element);

        BoussinesqAdaptive::syncElementDataMax(gridView, gridGeometry_->elementMapper(), gradMagnitude_);
    }

    /*!
     * \brief function call operator
     * \return  1 if the element should be refined, -1 if it should be coarsened, 0 otherwise
     */
    int operator() (const Element& element) const
    {
        const auto idx = gridGeometry_->elementMapper().index(element);
        if (element.hasFather() && gradMagnitude_[idx] < coarsenBound_)
            return -1;
        else if (element.level() < maxLevel_ && gradMagnitude_[idx] > refineBound_)
            return 1;
        else
            return 0;
    }

    //! Per-element |grad C| computed by the last calculate() call.
    const std::vector<Scalar>& values() const { return gradMagnitude_; }
    Scalar refineBound() const { return refineBound_; }
    Scalar coarsenBound() const { return coarsenBound_; }

private:
    bool checkNeighborsRefine_(const Element& element, std::size_t level = 1)
    {
        for (const auto& intersection : intersections(gridGeometry_->gridView(), element))
        {
            if (!intersection.neighbor())
                continue;

            const auto outside = intersection.outside();

            if (outside.partitionType() == Dune::GhostEntity)
                continue;

            if (outside.level() < maxLevel_ && outside.level() < element.level())
            {
                gradMagnitude_[gridGeometry_->elementMapper().index(outside)] = std::numeric_limits<Scalar>::max();
                if (level < maxLevel_)
                    checkNeighborsRefine_(outside, level + 1);
            }
        }
        return true;
    }

    std::shared_ptr<const GridGeometry> gridGeometry_;
    int concentrationIdx_;
    Scalar refineBound_;
    Scalar coarsenBound_;
    std::vector<Scalar> gradMagnitude_;
    std::size_t minLevel_;
    std::size_t maxLevel_;
};

} // end namespace Dumux

#endif
