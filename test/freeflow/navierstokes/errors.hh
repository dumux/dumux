// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/io/format.hh>
#include <dumux/geometry/diameter.hh>

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
                totalVolume += Extrusion::volume(fvGeometry, scv);

                // compute the pressure errors
                const auto analyticalSolutionCellCenter
                    = problem_->analyticalSolution(scv.dofPosition(), time)[Indices::pressureIdx];
                const auto numericalSolutionCellCenter
                    = curSol[GridGeometry::cellCenterIdx()][scv.dofIndex()][Indices::pressureIdx - dim];

                const Scalar pError = absDiff_(analyticalSolutionCellCenter, numericalSolutionCellCenter);
                const Scalar pReference = absDiff_(analyticalSolutionCellCenter, 0.0);

                maxError[Indices::pressureIdx] = std::max(maxError[Indices::pressureIdx], pError);
                maxReference[Indices::pressureIdx] = std::max(maxReference[Indices::pressureIdx], pReference);
                sumError[Indices::pressureIdx] += pError * pError * Extrusion::volume(fvGeometry, scv);
                sumReference[Indices::pressureIdx] += pReference * pReference * Extrusion::volume(fvGeometry, scv);

                for (const auto& scvf : scvfs(fvGeometry))
                {
                    // compute the velocity errors
                    const auto velocityIndex = Indices::velocity(scvf.directionIndex());

                    const Scalar staggeredHalfVolume = Extrusion::volume(fvGeometry,
                        typename GridGeometry::Traits::FaceSubControlVolume(
                            0.5*(scv.center() + scvf.center()), 0.5*Extrusion::volume(fvGeometry, scv)
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
                    << numCCDofs << ", " << numFaceDofs << ", " << numCCDofs + numFaceDofs << "\n";

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
        logFile << Fmt::vformat(", " + format, Fmt::make_format_args(error[Indices::pressureIdx]));
        for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << Fmt::vformat(", " + format, Fmt::make_format_args(error[Indices::velocity(dirIdx)]));
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
 * \brief Compute errors between an analytical solution and the numerical approximation for momentum eq.
 */
template<class Problem, class Scalar = double>
class ErrorsSubProblem
{
    static constexpr bool isMomentumProblem = Problem::isMomentumProblem();
    static constexpr int dim
    = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>::GridView::dimension;

    using ErrorVector = typename std::conditional_t< isMomentumProblem,
                                                     Dune::FieldVector<Scalar, dim>,
                                                     Dune::FieldVector<Scalar, 1> >;

public:
    template<class SolutionVector>
    ErrorsSubProblem(std::shared_ptr<const Problem> problem,
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
    const ErrorVector& l2Absolute() const { return l2Absolute_; }
    //! The relative discrete l2 error (relative to the discrete l2 norm of the reference solution)
    const ErrorVector& l2Relative() const { return l2Relative_; }
    //! The (absolute) discrete l-infinity error
    const ErrorVector& lInfAbsolute() const { return lInfAbsolute_; }
    //! The relative discrete l-infinity error (relative to the discrete loo norm of the reference solution)
    const ErrorVector& lInfRelative() const { return lInfRelative_; }

    //! Time corresponding to the error (returns 0 per default)
    Scalar time() const { return time_; }

    //! Maximum diameter of primal grid elements
    Scalar hMax() const { return hMax_; }

    //! Volume of domain
    Scalar totalVolume() const { return totalVolume_; }

    //! Number of scvs that have been considered in error calculation
    Scalar numDofs() const { return numDofs_; }

private:
    template<class SolutionVector>
    void calculateErrors_(const SolutionVector& curSol, Scalar time)
    {
        // store time information
        time_ = time;

        // calculate helping variables
        totalVolume_ = 0.0;
        hMax_ = 0.0;
        numDofs_ = problem_->gridGeometry().numDofs();

        using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
        // We do not consider the overlapping Dofs, i.e. elements, in the errors
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::pq1bubble)
            numDofs_ -= problem_->gridGeometry().gridView().size(0);

        ErrorVector sumReference(0.0);
        ErrorVector sumError(0.0);
        ErrorVector maxReference(0.0);
        ErrorVector maxError(0.0);

        using namespace Dune::Indices;

        auto fvGeometry = localView(problem_->gridGeometry());
        for (const auto& element : elements(problem_->gridGeometry().gridView()))
        {
            hMax_ = std::max(hMax_, diameter(element.geometry()));
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                using Extrusion = Extrusion_t<GridGeometry>;

                // velocity errors
                if constexpr (isMomentumProblem)
                {
                    if constexpr (GridGeometry::discMethod == DiscretizationMethods::fcstaggered)
                    {
                        totalVolume_ += Extrusion::volume(fvGeometry, scv);
                        // compute the velocity errors
                        using Indices = typename Problem::Indices;
                        const auto velIdx = Indices::velocity(scv.dofAxis());
                        const auto analyticalSolution
                            = problem_->analyticalSolution(scv.dofPosition(), time)[velIdx];
                        const auto numericalSolution
                            = curSol[scv.dofIndex()][0];

                        const Scalar vError = absDiff_(analyticalSolution, numericalSolution);
                        const Scalar vReference = absDiff_(analyticalSolution, 0.0);

                        maxError[velIdx] = std::max(maxError[velIdx], vError);
                        maxReference[velIdx] = std::max(maxReference[velIdx], vReference);
                        sumError[velIdx] += vError * vError * Extrusion::volume(fvGeometry, scv);
                        sumReference[velIdx] += vReference * vReference * Extrusion::volume(fvGeometry, scv);
                    }
                    else if constexpr (GridGeometry::discMethod == DiscretizationMethods::fcdiamond
                                       || GridGeometry::discMethod == DiscretizationMethods::box)
                    {
                        totalVolume_ += Extrusion::volume(fvGeometry, scv);
                        for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
                        {
                            const auto analyticalSolution
                                = problem_->analyticalSolution(scv.dofPosition(), time)[dirIdx];
                            const auto numericalSolution
                                = curSol[scv.dofIndex()][dirIdx];

                            const Scalar vError = absDiff_(analyticalSolution, numericalSolution);
                            const Scalar vReference = absDiff_(analyticalSolution, 0.0);

                            maxError[dirIdx] = std::max(maxError[dirIdx], vError);
                            maxReference[dirIdx] = std::max(maxReference[dirIdx], vReference);
                            sumError[dirIdx] += vError * vError * Extrusion::volume(fvGeometry, scv);
                            sumReference[dirIdx] += vReference * vReference * Extrusion::volume(fvGeometry, scv);
                        }
                    }
                    else if constexpr (GridGeometry::discMethod == DiscretizationMethods::pq1bubble)
                    {
                        if(!scv.isOverlapping())
                        {
                            totalVolume_ += Extrusion::volume(fvGeometry, scv);
                            for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
                            {
                                const auto analyticalSolution
                                    = problem_->analyticalSolution(scv.dofPosition(), time)[dirIdx];
                                const auto numericalSolution
                                    = curSol[scv.dofIndex()][dirIdx];

                                const Scalar vError = absDiff_(analyticalSolution, numericalSolution);
                                const Scalar vReference = absDiff_(analyticalSolution, 0.0);

                                maxError[dirIdx] = std::max(maxError[dirIdx], vError);
                                maxReference[dirIdx] = std::max(maxReference[dirIdx], vReference);
                                sumError[dirIdx] += vError * vError * Extrusion::volume(fvGeometry, scv);
                                sumReference[dirIdx] += vReference * vReference * Extrusion::volume(fvGeometry, scv);
                            }
                        }
                    }
                    else
                        DUNE_THROW(Dune::InvalidStateException,
                                   "Unknown momentum discretization scheme in Navier-Stokes error calculation");
                }
                // pressure errors
                else
                {
                    totalVolume_ += Extrusion::volume(fvGeometry, scv);
                    // compute the pressure errors
                    using Indices = typename Problem::Indices;
                    const auto analyticalSolution
                        = problem_->analyticalSolution(scv.dofPosition(), time)[Indices::pressureIdx];
                    const auto numericalSolution
                        = curSol[scv.dofIndex()][Indices::pressureIdx];

                    const Scalar pError = absDiff_(analyticalSolution, numericalSolution);
                    const Scalar pReference = absDiff_(analyticalSolution, 0.0);

                    maxError[0] = std::max(maxError[0], pError);
                    maxReference[0] = std::max(maxReference[0], pReference);
                    sumError[0] += pError * pError * Extrusion::volume(fvGeometry, scv);
                    sumReference[0] += pReference * pReference * Extrusion::volume(fvGeometry, scv);
                }
            }
        }

        // calculate errors
        for (int i = 0; i < ErrorVector::size(); ++i)
        {
            l2Absolute_[i] = std::sqrt(sumError[i] / totalVolume_);
            l2Relative_[i] = std::sqrt(sumError[i] / sumReference[i]);
            lInfAbsolute_[i] = maxError[i];
            lInfRelative_[i] = maxError[i] / maxReference[i];
        }
    }

    template<class T>
    T absDiff_(const T& a, const T& b) const
    { using std::abs; return abs(a-b); }

    std::shared_ptr<const Problem> problem_;

    ErrorVector l2Absolute_;
    ErrorVector l2Relative_;
    ErrorVector lInfAbsolute_;
    ErrorVector lInfRelative_;
    Scalar time_;
    Scalar hMax_;
    Scalar totalVolume_;
    std::size_t numDofs_;
};

/*!
 * \brief Compute errors between an analytical solution and the numerical approximation
 */
template<class MomentumProblem, class MassProblem, class Scalar = double>
class Errors
{
    static constexpr int dim
        = std::decay_t<decltype(std::declval<MomentumProblem>().gridGeometry())>::GridView::dimension;

    using ErrorVector = Dune::FieldVector<Scalar, dim+1>;
    using NumSubProblemVector = Dune::FieldVector<Scalar, 2>;

public:
    template<class SolutionVector>
    Errors(std::shared_ptr<const MomentumProblem> momentumProblem,
           std::shared_ptr<const MassProblem> massProblem,
           const SolutionVector& curSol,
           Scalar time = 0.0)
    : momentumErrors_(momentumProblem, curSol[Dune::Indices::_0], time)
    , massErrors_(massProblem, curSol[Dune::Indices::_1], time)
    { update(curSol, time); }

    /*!
     * \brief Computes errors between an analytical solution and the numerical approximation
     *
     * \param curSol The current solution vector
     * \param time The current time
     */
    template<class SolutionVector>
    void update(const SolutionVector& curSol, Scalar time = 0.0)
    {
        momentumErrors_.update(curSol[Dune::Indices::_0], time);
        massErrors_.update(curSol[Dune::Indices::_1], time);

        time_ = time;

        const auto&  l2AbsoluteMomentum = momentumErrors_.l2Absolute();
        const auto&  l2RelativeMomentum = momentumErrors_.l2Relative();
        const auto&  lInfAbsoluteMomentum = momentumErrors_.lInfAbsolute();
        const auto&  lInfRelativeMomentum = momentumErrors_.lInfRelative();

        const auto&  l2AbsoluteMass = massErrors_.l2Absolute();
        const auto&  l2RelativeMass = massErrors_.l2Relative();
        const auto&  lInfAbsoluteMass = massErrors_.lInfAbsolute();
        const auto&  lInfRelativeMass = massErrors_.lInfRelative();

        l2Absolute_[0] = l2AbsoluteMass[0];
        l2Relative_[0] = l2RelativeMass[0];
        lInfAbsolute_[0] = lInfAbsoluteMass[0];
        lInfRelative_[0] = lInfRelativeMass[0];

        std::copy( l2AbsoluteMomentum.begin(), l2AbsoluteMomentum.end(), l2Absolute_.begin() + 1 );
        std::copy( l2RelativeMomentum.begin(), l2RelativeMomentum.end(), l2Relative_.begin() + 1 );
        std::copy( lInfAbsoluteMomentum.begin(), lInfAbsoluteMomentum.end(), lInfAbsolute_.begin() + 1 );
        std::copy( lInfRelativeMomentum.begin(), lInfRelativeMomentum.end(), lInfRelative_.begin() + 1 );

        hMax_[0] = massErrors_.hMax();
        totalVolume_[0] = massErrors_.totalVolume();
        numDofs_[0] = massErrors_.numDofs();
        hMax_[1] = momentumErrors_.hMax();
        totalVolume_[1] = momentumErrors_.totalVolume();
        numDofs_[1] = momentumErrors_.numDofs();
    }

    //! The (absolute) discrete l2 error
    const ErrorVector& l2Absolute() const { return l2Absolute_; }
    //! The relative discrete l2 error (relative to the discrete l2 norm of the reference solution)
    const ErrorVector& l2Relative() const { return l2Relative_; }
    //! The (absolute) discrete l-infinity error
    const ErrorVector& lInfAbsolute() const { return lInfAbsolute_; }
    //! The relative discrete l-infinity error (relative to the discrete loo norm of the reference solution)
    const ErrorVector& lInfRelative() const { return lInfRelative_; }

    //! Volume of scvs considered in error calculation
    const NumSubProblemVector& totalVolume() const { return totalVolume_; }
    //! Number of scvs considered in error calculation
    const NumSubProblemVector& numDofs() const { return numDofs_; }
    //! Maximum diameter of primal grid elements
    const NumSubProblemVector& hMax() const { return hMax_; }

    //! Time corresponding to the error (returns 0 per default)
    Scalar time() const { return time_; }

private:
    ErrorsSubProblem<MomentumProblem, Scalar> momentumErrors_;
    ErrorsSubProblem<MassProblem, Scalar> massErrors_;

    ErrorVector l2Absolute_;
    ErrorVector l2Relative_;
    ErrorVector lInfAbsolute_;
    ErrorVector lInfRelative_;

    NumSubProblemVector hMax_;
    NumSubProblemVector totalVolume_;
    NumSubProblemVector numDofs_;

    Scalar time_;
};

template<class Problem, class SolutionVector>
ErrorsSubProblem(std::shared_ptr<Problem>, SolutionVector&&)
-> ErrorsSubProblem<Problem>;

template<class Problem, class SolutionVector, class Scalar>
ErrorsSubProblem(std::shared_ptr<Problem>, SolutionVector&&, Scalar)
-> ErrorsSubProblem<Problem, Scalar>;

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
                    << numCCDofs << ", " << numFaceDofs << ", " << numCCDofs + numFaceDofs << "\n";

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
        logFile << Fmt::vformat(", " + format, Fmt::make_format_args(error[MassIndices::pressureIdx]));
        for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
            logFile << Fmt::vformat(", " + format, Fmt::make_format_args(error[MomIndices::velocity(dirIdx)+1]));
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

/*!
 * \brief Append errors to the log file during a convergence test
 */
template<class MomentumProblem>
void convergenceTestAppendErrorsMomentum(std::ofstream& logFile,
                                         std::shared_ptr<MomentumProblem> problem,
                                         const ErrorsSubProblem<MomentumProblem>& errors)
{
    const auto numFaceDofs = problem->gridGeometry().numDofs();

    logFile << Fmt::format(
        "[ConvergenceTest] numFaceDofs = {}",
        numFaceDofs
    );

    const auto print = [&](const auto& e, const std::string& name){
        using MomIndices = typename MomentumProblem::Indices;
        logFile << Fmt::format(
            " {0}(vx) = {1} {0}(vy) = {2}",
            name, e[MomIndices::velocityXIdx], e[MomIndices::velocityYIdx]
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
template<class MomentumProblem>
void convergenceTestAppendErrorsMomentum(std::shared_ptr<MomentumProblem> problem,
                                         const ErrorsSubProblem<MomentumProblem>& errors)
{
    const auto logFileName = problem->name() + ".log";
    std::ofstream logFile(logFileName, std::ios::app);
    convergenceTestAppendErrorsMomentum(logFile, problem, errors);
}

} // end namespace NavierStokesTest

} // end namespace Dumux

#endif
