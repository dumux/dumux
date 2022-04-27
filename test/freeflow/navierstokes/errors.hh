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
 * \ingroup NavierStokesTests
 * \copydoc Dumux::NavierStokesTestError
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ERRORS_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ERRORS_HH

#include <vector>
#include <cmath>
#include <fstream>
#include <dune/common/indices.hh>
#include <dune/common/fvector.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/io/format.hh>

namespace Dumux {

/*!
 * \brief Compute errors between an analytical solution and the numerical approximation
 */
template<class Problem, class Scalar = double>
class NavierStokesErrors
{
    using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Extrusion = Extrusion_t<GridGeometry>;
    using Indices = typename Problem::Indices;
    using Element = typename GridGeometry::LocalView::Element;
    using PrimaryVariables = std::decay_t<decltype(
        std::declval<Problem>().dirichlet(std::declval<Element>(), std::declval<SubControlVolume>())
    )>;

    static constexpr int dim = GridGeometry::GridView::dimension;

public:
    template<class SolutionVector>
    NavierStokesErrors(std::shared_ptr<const Problem> problem,
                       const SolutionVector& curSol,
                       Scalar time = 0.0)
    : problem_(problem)
    { calculateErrors_(curSol, time); }

    /*!
     * \brief Computes errors between an analytical solution and the numerical approximation
     *
     * \param curSol The current solution vector
     * \param time The current time
     */
    template<class SolutionVector>
    void update(const SolutionVector& curSol, Scalar time = 0.0)
    { calculateErrors_(curSol, time); }

    //! The (absolute) discrete l2 error
    const PrimaryVariables& l2Absolute() const { return l2Absolute_; }
    //! The relative discrete l2 error (relative to the discrete l2 norm of the reference solution)
    const PrimaryVariables& l2Relative() const { return l2Relative_; }
    //! The (absolute) discrete l-infinity error
    const PrimaryVariables& lInfAbsolute() const { return lInfAbsolute_; }
    //! The relative discrete l-infinity error (relative to the discrete loo norm of the reference solution)
    const PrimaryVariables& lInfRelative() const { return lInfRelative_; }

    //! Time corresponding to the error (returns 0 per default)
    Scalar time() const { return time_; }

private:
    template<class SolutionVector>
    void calculateErrors_(const SolutionVector& curSol, Scalar time)
    {
        // store time information
        time_ = time;

        // calculate helping variables
        Scalar totalVolume = 0.0;

        PrimaryVariables sumReference(0.0);
        PrimaryVariables sumError(0.0);
        PrimaryVariables maxReference(0.0);
        PrimaryVariables maxError(0.0);

        auto fvGeometry = localView(problem_->gridGeometry());
        for (const auto& element : elements(problem_->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                totalVolume += Extrusion::volume(scv);

                // compute the pressure errors
                const auto analyticalSolutionCellCenter
                    = problem_->analyticalSolution(scv.dofPosition(), time)[Indices::pressureIdx];
                const auto numericalSolutionCellCenter
                    = curSol[GridGeometry::cellCenterIdx()][scv.dofIndex()][Indices::pressureIdx - dim];

                const Scalar pError = absDiff_(analyticalSolutionCellCenter, numericalSolutionCellCenter);
                const Scalar pReference = absDiff_(analyticalSolutionCellCenter, 0.0);

                maxError[Indices::pressureIdx] = std::max(maxError[Indices::pressureIdx], pError);
                maxReference[Indices::pressureIdx] = std::max(maxReference[Indices::pressureIdx], pReference);
                sumError[Indices::pressureIdx] += pError * pError * Extrusion::volume(scv);
                sumReference[Indices::pressureIdx] += pReference * pReference * Extrusion::volume(scv);

                for (const auto& scvf : scvfs(fvGeometry))
                {
                    // compute the velocity errors
                    const auto velocityIndex = Indices::velocity(scvf.directionIndex());

                    const Scalar staggeredHalfVolume = Extrusion::volume(
                        typename GridGeometry::Traits::FaceSubControlVolume(
                            0.5*(scv.center() + scvf.center()), 0.5*Extrusion::volume(scv)
                        )
                    );

                    const auto analyticalSolutionFace = problem_->analyticalSolution(scvf.center(), time)[velocityIndex];
                    const auto numericalSolutionFace = curSol[GridGeometry::faceIdx()][scvf.dofIndex()][0];

                    const Scalar velError = absDiff_(analyticalSolutionFace, numericalSolutionFace);
                    const Scalar velReference = absDiff_(analyticalSolutionFace, 0.0);

                    maxError[velocityIndex] = std::max(maxError[velocityIndex], velError);
                    maxReference[velocityIndex] = std::max(maxReference[velocityIndex], velReference);
                    sumError[velocityIndex] += velError * velError * staggeredHalfVolume;
                    sumReference[velocityIndex] += velReference * velReference * staggeredHalfVolume;
                }
            }
        }

        // calculate errors
        for (int i = 0; i < PrimaryVariables::size(); ++i)
        {
            l2Absolute_[i] = std::sqrt(sumError[i] / totalVolume);
            l2Relative_[i] = std::sqrt(sumError[i] / sumReference[i]);
            lInfAbsolute_[i] = maxError[i];
            lInfRelative_[i] = maxError[i] / maxReference[i];
        }
    }

    template<class T>
    T absDiff_(const T& a, const T& b) const
    { using std::abs; return abs(a-b); }

    std::shared_ptr<const Problem> problem_;

    PrimaryVariables l2Absolute_;
    PrimaryVariables l2Relative_;
    PrimaryVariables lInfAbsolute_;
    PrimaryVariables lInfRelative_;
    Scalar time_;
};

template<class Problem, class SolutionVector>
NavierStokesErrors(std::shared_ptr<Problem>, SolutionVector&&)
-> NavierStokesErrors<Problem>;

template<class Problem, class SolutionVector, class Scalar>
NavierStokesErrors(std::shared_ptr<Problem>, SolutionVector&&, Scalar)
-> NavierStokesErrors<Problem, Scalar>;


/*!
 * \brief An error CSV file writer
 */
template<class Problem, class Scalar = double>
class NavierStokesErrorCSVWriter
{
    using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using Indices = typename Problem::Indices;

    static constexpr int dim = GridGeometry::GridView::dimension;

public:
    NavierStokesErrorCSVWriter(std::shared_ptr<const Problem> problem, const std::string& suffix = "")
    : name_(problem->name() + (suffix.empty() ? "_error" : "_error_" + suffix))
    {
        const int numCCDofs = problem->gridGeometry().numCellCenterDofs();
        const int numFaceDofs = problem->gridGeometry().numFaceDofs();

        // print auxiliary file with the number of dofs
        std::ofstream logFileDofs(name_ + "_dofs.csv", std::ios::trunc);
        logFileDofs << "cc dofs, face dofs, all dofs\n"
                    << numCCDofs << numFaceDofs << numCCDofs + numFaceDofs << "\n";

        // clear error file
        std::ofstream logFile(name_ + ".csv", std::ios::trunc);

        // write header
        logFile << "time";
        for (const std::string e : { "L2Abs", "L2Rel", "LinfAbs", "LinfRel" })
            printError_(logFile, errorNames_(e), "{:s}");
        logFile << "\n";
    }

    void printErrors(const NavierStokesErrors<Problem>& errors) const
    {
        // append to error file
        std::ofstream logFile(name_ + ".csv", std::ios::app);

        logFile << Fmt::format("{:.5e}", errors.time());
        printError_(logFile, errors.l2Absolute());
        printError_(logFile, errors.l2Relative());
        printError_(logFile, errors.lInfAbsolute());
        printError_(logFile, errors.lInfRelative());

        logFile << "\n";
    }

private:
    template<class Error>
    void printError_(std::ofstream& logFile, const Error& error, const std::string& format = "{:.5e}") const
    {
        logFile << Fmt::format(", " + format, error[Indices::pressureIdx]);
        for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << Fmt::format(", " + format, error[Indices::velocity(dirIdx)]);
    }

    std::vector<std::string> errorNames_(const std::string& e) const
    {
        if constexpr (dim == 1)
            return { e + "(u)", e + "(p)" };
        else if constexpr (dim == 2)
            return { e + "(u)", e + "(v)", e + "(p)" };
        else
            return { e + "(u)", e + "(v)", e + "(w)", e + "(p)" };
    }

    std::string name_;
};

template<class Problem>
NavierStokesErrorCSVWriter(std::shared_ptr<Problem>)
-> NavierStokesErrorCSVWriter<Problem>;

template<class Problem>
NavierStokesErrorCSVWriter(std::shared_ptr<Problem>, const std::string&)
-> NavierStokesErrorCSVWriter<Problem>;


/*!
 * \brief Append errors to the log file during a convergence test
 */
template<class Problem>
void convergenceTestAppendErrors(std::ofstream& logFile,
                                 std::shared_ptr<Problem> problem,
                                 const NavierStokesErrors<Problem>& errors)
{
    const auto numCCDofs = problem->gridGeometry().numCellCenterDofs();
    const auto numFaceDofs = problem->gridGeometry().numFaceDofs();

    logFile << Fmt::format(
        "[ConvergenceTest] numCCDofs = {} numFaceDofs = {}",
        numCCDofs, numFaceDofs
    );

    const auto print = [&](const auto& e, const std::string& name){
        using Indices = typename Problem::Indices;
        logFile << Fmt::format(
            " {0}(p) = {1} {0}(vx) = {2} {0}(vy) = {3}",
            name, e[Indices::pressureIdx], e[Indices::velocityXIdx], e[Indices::velocityYIdx]
        );
    };

    print(errors.l2Absolute(), "L2Abs");
    print(errors.l2Relative(), "L2Rel");
    print(errors.lInfAbsolute(), "LinfAbs");
    print(errors.lInfRelative(), "LinfRel");

    logFile << "\n";
}

/*!
 * \brief Append errors to the log file during a convergence test
 */
template<class Problem>
void convergenceTestAppendErrors(std::shared_ptr<Problem> problem,
                                 const NavierStokesErrors<Problem>& errors)
{
    const auto logFileName = problem->name() + ".log";
    std::ofstream logFile(logFileName, std::ios::app);
    convergenceTestAppendErrors(logFile, problem, errors);
}

namespace NavierStokesTest {

/*!
 * \brief Compute errors between an analytical solution and the numerical approximation
 */
template<class MomentumProblem, class MassProblem, class Scalar = double>
class Errors
{
    static constexpr int dim
        = std::decay_t<decltype(std::declval<MomentumProblem>().gridGeometry())>::GridView::dimension;

    using ErrorVector = Dune::FieldVector<Scalar, dim+1>;

public:
    template<class SolutionVector>
    Errors(std::shared_ptr<const MomentumProblem> momentumProblem,
           std::shared_ptr<const MassProblem> massProblem,
           const SolutionVector& curSol,
           Scalar time = 0.0)
    : momentumProblem_(momentumProblem)
    , massProblem_(massProblem)
    { calculateErrors_(curSol, time); }

    /*!
     * \brief Computes errors between an analytical solution and the numerical approximation
     *
     * \param curSol The current solution vector
     * \param time The current time
     */
    template<class SolutionVector>
    void update(const SolutionVector& curSol, Scalar time = 0.0)
    { calculateErrors_(curSol, time); }

    //! The (absolute) discrete l2 error
    const ErrorVector& l2Absolute() const { return l2Absolute_; }
    //! The relative discrete l2 error (relative to the discrete l2 norm of the reference solution)
    const ErrorVector& l2Relative() const { return l2Relative_; }
    //! The (absolute) discrete l-infinity error
    const ErrorVector& lInfAbsolute() const { return lInfAbsolute_; }
    //! The relative discrete l-infinity error (relative to the discrete loo norm of the reference solution)
    const ErrorVector& lInfRelative() const { return lInfRelative_; }

    //! Time corresponding to the error (returns 0 per default)
    Scalar time() const { return time_; }

private:
    template<class SolutionVector>
    void calculateErrors_(const SolutionVector& curSol, Scalar time)
    {
        // store time information
        time_ = time;

        // calculate helping variables
        Scalar totalVolume = 0.0;

        ErrorVector sumReference(0.0);
        ErrorVector sumError(0.0);
        ErrorVector maxReference(0.0);
        ErrorVector maxError(0.0);

        using namespace Dune::Indices;

        // velocity errors
        {
            auto fvGeometry = localView(momentumProblem_->gridGeometry());
            for (const auto& element : elements(momentumProblem_->gridGeometry().gridView()))
            {
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                {
                    using GridGeometry = std::decay_t<decltype(std::declval<MassProblem>().gridGeometry())>;
                    using Extrusion = Extrusion_t<GridGeometry>;

                    // compute the pressure errors
                    using Indices = typename MomentumProblem::Indices;
                    const auto velIdx = Indices::velocity(scv.dofAxis());
                    const auto analyticalSolution
                        = momentumProblem_->analyticalSolution(scv.dofPosition(), time)[velIdx];
                    const auto numericalSolution
                        = curSol[_0][scv.dofIndex()][0];

                    const Scalar vError = absDiff_(analyticalSolution, numericalSolution);
                    const Scalar vReference = absDiff_(analyticalSolution, 0.0);

                    maxError[velIdx+1] = std::max(maxError[velIdx+1], vError);
                    maxReference[velIdx+1] = std::max(maxReference[velIdx+1], vReference);
                    sumError[velIdx+1] += vError * vError * Extrusion::volume(scv);
                    sumReference[velIdx+1] += vReference * vReference * Extrusion::volume(scv);
                }
            }
        }

        // pressure errors
        {
            auto fvGeometry = localView(massProblem_->gridGeometry());
            for (const auto& element : elements(massProblem_->gridGeometry().gridView()))
            {
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                {
                    using GridGeometry = std::decay_t<decltype(std::declval<MomentumProblem>().gridGeometry())>;
                    using Extrusion = Extrusion_t<GridGeometry>;
                    totalVolume += Extrusion::volume(scv);

                    // compute the pressure errors
                    using Indices = typename MassProblem::Indices;
                    const auto analyticalSolution
                        = massProblem_->analyticalSolution(scv.dofPosition(), time)[Indices::pressureIdx];
                    const auto numericalSolution
                        = curSol[_1][scv.dofIndex()][Indices::pressureIdx];

                    const Scalar pError = absDiff_(analyticalSolution, numericalSolution);
                    const Scalar pReference = absDiff_(analyticalSolution, 0.0);

                    maxError[0] = std::max(maxError[0], pError);
                    maxReference[0] = std::max(maxReference[0], pReference);
                    sumError[0] += pError * pError * Extrusion::volume(scv);
                    sumReference[0] += pReference * pReference * Extrusion::volume(scv);
                }
            }
        }

        // calculate errors
        for (int i = 0; i < ErrorVector::size(); ++i)
        {
            l2Absolute_[i] = std::sqrt(sumError[i] / totalVolume);
            l2Relative_[i] = std::sqrt(sumError[i] / sumReference[i]);
            lInfAbsolute_[i] = maxError[i];
            lInfRelative_[i] = maxError[i] / maxReference[i];
        }
    }

    template<class T>
    T absDiff_(const T& a, const T& b) const
    { using std::abs; return abs(a-b); }

    std::shared_ptr<const MomentumProblem> momentumProblem_;
    std::shared_ptr<const MassProblem> massProblem_;

    ErrorVector l2Absolute_;
    ErrorVector l2Relative_;
    ErrorVector lInfAbsolute_;
    ErrorVector lInfRelative_;
    Scalar time_;
};

template<class MomentumProblem, class MassProblem, class SolutionVector>
Errors(std::shared_ptr<MomentumProblem>, std::shared_ptr<MassProblem>, SolutionVector&&)
-> Errors<MomentumProblem, MassProblem>;

template<class MomentumProblem, class MassProblem, class SolutionVector, class Scalar>
Errors(std::shared_ptr<MomentumProblem>, std::shared_ptr<MassProblem>, SolutionVector&&, Scalar)
-> Errors<MomentumProblem, MassProblem, Scalar>;


/*!
 * \brief An error CSV file writer
 */
template<class MomentumProblem, class MassProblem, class Scalar = double>
class ErrorCSVWriter
{
    static constexpr int dim
        = std::decay_t<decltype(std::declval<MomentumProblem>().gridGeometry())>::GridView::dimension;

public:
    ErrorCSVWriter(std::shared_ptr<const MomentumProblem> momentumProblem,
                   std::shared_ptr<const MassProblem> massProblem,
                   const std::string& suffix = "")
    : name_(massProblem->name() + (suffix.empty() ? "_error" : "_error_" + suffix))
    {
        const int numCCDofs = massProblem->gridGeometry().numDofs();
        const int numFaceDofs = momentumProblem->gridGeometry().numDofs();

        // print auxiliary file with the number of dofs
        std::ofstream logFileDofs(name_ + "_dofs.csv", std::ios::trunc);
        logFileDofs << "cc dofs, face dofs, all dofs\n"
                    << numCCDofs << numFaceDofs << numCCDofs + numFaceDofs << "\n";

        // clear error file
        std::ofstream logFile(name_ + ".csv", std::ios::trunc);

        // write header
        logFile << "time";
        using ErrorNames = std::vector<std::string>;
        for (const std::string e : { "L2Abs", "L2Rel", "LinfAbs", "LinfRel" })
            printError_(logFile, ErrorNames({ e + "(p)", e + "(u)", e + "(v)", e + "(w)" }), "{:s}");
        logFile << "\n";
    }

    void printErrors(const Errors<MomentumProblem, MassProblem>& errors) const
    {
        // append to error file
        std::ofstream logFile(name_ + ".csv", std::ios::app);

        logFile << Fmt::format("{:.5e}", errors.time());
        printError_(logFile, errors.l2Absolute());
        printError_(logFile, errors.l2Relative());
        printError_(logFile, errors.lInfAbsolute());
        printError_(logFile, errors.lInfRelative());

        logFile << "\n";
    }

private:
    template<class Error>
    void printError_(std::ofstream& logFile, const Error& error, const std::string& format = "{:.5e}") const
    {
        using MassIndices = typename MassProblem::Indices;
        using MomIndices = typename MomentumProblem::Indices;
        logFile << Fmt::format(", " + format, error[MassIndices::pressureIdx]);
        for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << Fmt::format(", " + format, error[MomIndices::velocity(dirIdx)]);
    }

    std::string name_;
};

template<class MomentumProblem, class MassProblem>
ErrorCSVWriter(std::shared_ptr<MomentumProblem>, std::shared_ptr<MassProblem>)
-> ErrorCSVWriter<MomentumProblem, MassProblem>;

template<class MomentumProblem, class MassProblem>
ErrorCSVWriter(std::shared_ptr<MomentumProblem>, std::shared_ptr<MassProblem>, const std::string&)
-> ErrorCSVWriter<MomentumProblem, MassProblem>;


/*!
 * \brief Append errors to the log file during a convergence test
 */
template<class MomentumProblem, class MassProblem>
void convergenceTestAppendErrors(std::ofstream& logFile,
                                 std::shared_ptr<MomentumProblem> momentumProblem,
                                 std::shared_ptr<MassProblem> massProblem,
                                 const Errors<MomentumProblem, MassProblem>& errors)
{
    const auto numCCDofs = massProblem->gridGeometry().numDofs();
    const auto numFaceDofs = momentumProblem->gridGeometry().numDofs();

    logFile << Fmt::format(
        "[ConvergenceTest] numCCDofs = {} numFaceDofs = {}",
        numCCDofs, numFaceDofs
    );

    const auto print = [&](const auto& e, const std::string& name){
        using MassIndices = typename MassProblem::Indices;
        using MomIndices = typename MomentumProblem::Indices;
        logFile << Fmt::format(
            " {0}(p) = {1} {0}(vx) = {2} {0}(vy) = {3}",
            name, e[MassIndices::pressureIdx], e[MomIndices::velocityXIdx], e[MomIndices::velocityYIdx]
        );
    };

    print(errors.l2Absolute(), "L2Abs");
    print(errors.l2Relative(), "L2Rel");
    print(errors.lInfAbsolute(), "LinfAbs");
    print(errors.lInfRelative(), "LinfRel");

    logFile << "\n";
}


/*!
 * \brief Append errors to the log file during a convergence test
 */
template<class MomentumProblem, class MassProblem>
void convergenceTestAppendErrors(std::shared_ptr<MomentumProblem> momentumProblem,
                                 std::shared_ptr<MassProblem> massProblem,
                                 const Errors<MomentumProblem, MassProblem>& errors)
{
    const auto logFileName = massProblem->name() + ".log";
    std::ofstream logFile(logFileName, std::ios::app);
    convergenceTestAppendErrors(logFile, momentumProblem, massProblem, errors);
}

} // end namespace NavierStokesTest

} // end namespace Dumux

#endif
