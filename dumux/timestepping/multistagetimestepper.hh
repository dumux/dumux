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
 * \brief A time stepper performing a single time step of a transient simulation
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_TIMESTEPPER_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_TIMESTEPPER_HH

#include <iterator>
#include <algorithm>
#include <memory>
#include <vector>

namespace Dumux {

//! forward declaration
template<class Scalar>
class MultiStageMethod;

//! Data object for the parameters of a given stage
template<class Scalar>
class MultiStageParams
{
    struct Params {
        Scalar alpha, betaDt, timeAtStage, dtFraction;
        bool skipTemporal, skipSpatial;
    };
public:
    //! Extract params for stage i from method m
    MultiStageParams(const MultiStageMethod<Scalar>& m, std::size_t i, const Scalar t, const Scalar dt)
    : size_(i+1)
    {
        params_.resize(size_);
        for (std::size_t k = 0; k < size_; ++k)
        {
            auto& p = params_[k];
            p.alpha = m.temporalWeight(i, k);
            p.betaDt = m.spatialWeight(i, k)*dt;
            p.timeAtStage = t + m.timeStepWeight(k)*dt;
            p.dtFraction = m.timeStepWeight(k);

            using std::abs;
            p.skipTemporal = (abs(p.alpha) < 1e-6);
            p.skipSpatial = (abs(p.betaDt) < 1e-6);
        }
    }

    std::size_t size () const
    { return size_; }

    //! weights of the temporal operator residual (\f$ \alpha_{ik} \f$)
    Scalar temporalWeight (std::size_t k) const
    { return params_[k].alpha; }

    //! weights of the spatial operator residual (\f$ \beta_{ik} \f$)
    Scalar spatialWeight (std::size_t k) const
    { return params_[k].betaDt; }

    //! the time at which we have to evaluate the operators
    Scalar timeAtStage (std::size_t k) const
    { return params_[k].timeAtStage; }

    //! the fraction of a time step corresponding to the k-th stage
    Scalar timeStepFraction (std::size_t k) const
    { return params_[k].dtFraction; }

    //! If \f$ \alpha_{ik} = 0\f$
    Scalar skipTemporal (std::size_t k) const
    { return params_[k].skipTemporal; }

    //! If \f$ \beta_{ik} = 0\f$
    Scalar skipSpatial (std::size_t k) const
    { return params_[k].skipSpatial; }

private:
    std::size_t size_;
    std::vector<Params> params_;
};

/*!
 * \brief Solution vectors for the multi-stage time stepping scheme
 */
template<class SolutionVector>
class MultiStageSolutions
{
public:
    MultiStageSolutions(std::size_t numStages)
    {
        // avoid resizing vectors after construction
        sols_.reserve(numStages);
        owned_.reserve(numStages);
    }

    MultiStageSolutions(const MultiStageSolutions&) = delete;
    MultiStageSolutions(MultiStageSolutions&&) = delete;
    MultiStageSolutions& operator=(const MultiStageSolutions&) = delete;
    MultiStageSolutions& operator=(MultiStageSolutions&&) = delete;

    ~MultiStageSolutions()
    {
        // destroy all owned solution resources
        for (int i = 0; i < sols_.size(); ++i)
            if (owned_[i])
                delete sols_[i];
    }

    //! push an existing solution
    void pushSolution(SolutionVector& sol)
    {
        sols_.push_back(&sol);
        owned_.push_back(false);
    }

    //! create an intermediate/temporary solution
    void createSolution(const SolutionVector& guess)
    {
        sols_.push_back(new SolutionVector(guess));
        owned_.push_back(true);
    }

    auto size() const
    { return sols_.size(); }

    SolutionVector& operator[] (std::size_t i)
    { return *sols_[i]; }

    const SolutionVector& operator[] (std::size_t i) const
    { return *sols_[i]; }

    SolutionVector& back()
    { return *sols_.back(); }

    const SolutionVector& back() const
    { return *sols_.back(); }

private:
    std::vector<SolutionVector*> sols_;
    std::vector<bool> owned_;
};


/*!
 * \brief Time stepping with a multi-stage method
 * \note We limit ourselves to "diagonally" implicit multi-stage methods where solving
 *       a stage can only depend on the values of the same stage and stages before
 *       but not future stages (which would require solving larger linear systems)
 */
template<class PDESolver>
class MultiStageTimeStepper
{
    using Scalar = typename PDESolver::Assembler::Scalar;
    using SolutionVector = typename PDESolver::Assembler::ResidualType;
    using Solutions = MultiStageSolutions<SolutionVector>;

public:

    /*!
     * \brief The constructor
     * \param pdeSolver Solver class for solving a PDE in each stage
     * \param msMethod The multi-stage method which is to be used for time integration
     * \todo TODO: Add time step control if the pde solver doesn't converge
     */
    MultiStageTimeStepper(std::shared_ptr<PDESolver> pdeSolver,
                          std::shared_ptr<const MultiStageMethod<Scalar>> msMethod)
    : pdeSolver_(pdeSolver)
    , msMethod_(msMethod)
    {}

    /*!
     * \brief Advance one time step of the given time loop
     * \param x The container in which to write the solution of the next time step.
     *          This is also used as the initial guess.
     * \param t The current time level
     * \param dt The time step size to be used
     * TODO: Add time step control if the pde solver doesn't converge
     */
    void step(SolutionVector& x, const Scalar t, const Scalar dt)
    {
        const auto numStages = msMethod_->numStages();
        Solutions solutions(numStages);

        // The solution of stage zero is oldOld (\f$ x^{n}\f$).
        solutions.createSolution(x);

        // make sure the assembler gets rid of stuff from the last time integration
        pdeSolver_->assembler().prepareTimeIntegration(numStages);

        // For each intermediate stage we solve for the stage solution.
        // In the last stage (stageIdx = numStages), we obtain the solution (\f$ x^{n+1}\f$)
        // and write it into the provided container
        for (auto stageIdx = 1UL; stageIdx <= numStages; ++stageIdx)
        {
            if (stageIdx == numStages)
            {
                solutions.pushSolution(x);
                solutions.back() = solutions[stageIdx-1];
            }
            else
                solutions.createSolution(solutions[stageIdx-1]);

            // extract parameters for this stage from the time stepping method
            auto stageParams = std::make_shared<MultiStageParams<Scalar>>(*msMethod_, stageIdx, t, dt);

            // prepare the assembler for this stage
            pdeSolver_->assembler().prepareStage(solutions, stageParams);

            // assemble & solve
            pdeSolver_->solve(solutions.back());
        }
    }

    /*!
     * \brief Set/change the time step method
     */
    void setMethod(std::shared_ptr<const MultiStageMethod<Scalar>> msMethod)
    { msMethod_ = msMethod; }

private:
    std::shared_ptr<PDESolver> pdeSolver_;
    std::shared_ptr<const MultiStageMethod<Scalar>> msMethod_;
};

} // end namespace Dumux

#endif
