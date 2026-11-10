// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \copydoc Dumux::NavierStokesTestError
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ERRORS_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ERRORS_HH

#include <vector>
#include <cmath>
#include <fstream>
#include <tuple>

#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/io/format.hh>
#include <dumux/geometry/diameter.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/common/typetraits/griddiscretization.hh>

#include <test/testing/convergencetest.hh>

namespace Dumux {

/*!
 * \brief Velocity labels used for printing errors
 */
std::vector<std::string> velocityLabels()
{ return {"(vx)", "(vy)", "(vz)"}; }

namespace NavierStokesTest {

/*!
 * \brief Compute errors between an analytical solution and the numerical approximation for momentum eq.
 */
template<class Problem, class Scalar = double>
class ErrorsSubProblem
{
    static constexpr bool isMomentumProblem = Problem::isMomentumProblem();
    static constexpr int dim
    = ProblemTraits<Problem>::GridGeometry::GridView::dimension;

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

        using GridGeometry = typename ProblemTraits<Problem>::GridGeometry;
        const auto& gridGeometry = Dumux::gridDiscretization(*problem_);

        numDofs_ = gridGeometry.numDofs();

        // We do not consider the overlapping Dofs, i.e. elements, in the errors
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::pq1bubble
                      || GridGeometry::discMethod == DiscretizationMethods::pq2)
            numDofs_ -= gridGeometry.gridView().size(0);

        ErrorVector sumReference(0.0);
        ErrorVector sumError(0.0);
        ErrorVector maxReference(0.0);
        ErrorVector maxError(0.0);

        using namespace Dune::Indices;

        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
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
                    else if constexpr (GridGeometry::discMethod == DiscretizationMethods::pq1bubble
                                       || GridGeometry::discMethod == DiscretizationMethods::pq2)
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
        = ProblemTraits<MomentumProblem>::GridGeometry::GridView::dimension;

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
        = ProblemTraits<MomentumProblem>::GridGeometry::GridView::dimension;

public:
    ErrorCSVWriter(std::shared_ptr<const MomentumProblem> momentumProblem,
                   std::shared_ptr<const MassProblem> massProblem,
                   const std::string& suffix = "")
    : name_(massProblem->name() + (suffix.empty() ? "_error" : "_error_" + suffix))
    {
        const auto& massGridGeometry = Dumux::gridDiscretization(*massProblem);

        const auto& momentumGridGeometry = Dumux::gridDiscretization(*momentumProblem);

        const int numCCDofs = massGridGeometry.numDofs();
        const int numFaceDofs = momentumGridGeometry.numDofs();

        // print auxiliary file with the number of dofs
        std::ofstream logFileDofs(name_ + "_dofs.csv", std::ios::trunc);
        logFileDofs << "cc dofs, face dofs, all dofs\n"
                    << numCCDofs << ", " << numFaceDofs << ", " << numCCDofs + numFaceDofs << "\n";

        // clear error file
        std::ofstream logFile(name_ + ".csv", std::ios::trunc);

        const auto velLabels = velocityLabels();

        // write header
        logFile << "time";
        using ErrorNames = std::vector<std::string>;
        for (const std::string e : { "L2Abs", "L2Rel", "LinfAbs", "LinfRel" })
            printError_(logFile, ErrorNames({ e + "(p)", e + velLabels[0], e + velLabels[1], e + velLabels[2] }), "{:s}");
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
 * \brief The L2 norm and the H1 norm of the error of a solution field, integrated with a quadrature rule of the given order
 * \return The integrated volume and the two norms
 */
template<class Problem, class GridVariables, class SolutionVector>
std::tuple<double, Dune::FieldVector<double, 2>> calculateL2AndH1Errors(const Problem& problem,
                                                                        const GridVariables& gridVariables,
                                                                        const SolutionVector& x,
                                                                        int order = 5)
{
    using GridGeometry = Dumux::GridDiscretization_t<GridVariables>;
    using Extrusion = Extrusion_t<GridGeometry>;
    double totalVolume = 0.0;
    Dune::FieldVector<double, 2> errors(0.0);
    const auto& gridDiscretization = Dumux::gridDiscretization(problem);
    auto fvGeometry = localView(gridDiscretization);

    auto curGridVars = [&]() -> decltype(auto)
    {
        if constexpr (requires { gridVariables.curGridVolVars(); })
            return gridVariables.curGridVolVars();
        else
            return gridVariables.curGridVars();
    };

    auto elemVolVars = localView(curGridVars());
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        fvGeometry.bind(element);
        const auto geometry = fvGeometry.elementGeometry();

        elemVolVars.bind(element, fvGeometry, x);
        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        const auto& quad = Dune::QuadratureRules<double, GridGeometry::GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            const auto& localPos = qp.position();
            const auto& qpVolumeWeight = qp.weight() * Extrusion::integrationElement(geometry, localPos);
            totalVolume += qpVolumeWeight;

            const auto& globalPos = geometry.global(localPos);
            const auto analyticalSolution = problem.analyticalSolution(globalPos);
            const auto numericalSolution = evalSolutionAtLocalPos(element, geometry, gridDiscretization, elemSol, localPos);
            const auto solDiff = numericalSolution - analyticalSolution;
            errors[0] += (solDiff * solDiff) * qpVolumeWeight;

            const auto gradAnalyticalSolution = problem.gradAnalyticalSolution(globalPos);
            const auto gradNumericalSolution = evalGradientsAtLocalPos(element, geometry, gridDiscretization, elemSol, localPos);
            const auto gradDiff = gradNumericalSolution - gradAnalyticalSolution;
            double gradDiffNorm = 0.0;
            for (int i = 0; i < gradDiff.size(); ++i)
                gradDiffNorm += gradDiff[i] * gradDiff[i];

            errors[1] += (solDiff * solDiff + gradDiffNorm) * qpVolumeWeight;
        }
    }
    errors[0] = std::sqrt(errors[0]);
    errors[1] = std::sqrt(errors[1]);

    return {totalVolume, errors};
}

/*!
 * \brief The L2 norm of the error of a solution field, integrated with a quadrature rule of the given order
 * \return The integrated volume and the norm
 */
template<class Problem, class GridVariables, class SolutionVector>
std::tuple<double, double> calculateL2Error(const Problem& problem,
                                            const GridVariables& gridVariables,
                                            const SolutionVector& x,
                                            int order = 5)
{
    using GridGeometry = Dumux::GridDiscretization_t<GridVariables>;
    using Extrusion = Extrusion_t<GridGeometry>;
    double totalVolume = 0.0;
    double error = 0.0;
    const auto& gridDiscretization = Dumux::gridDiscretization(problem);
    auto fvGeometry = localView(gridDiscretization);

    auto curGridVars = [&]() -> decltype(auto)
    {
        if constexpr (requires { gridVariables.curGridVolVars(); })
            return gridVariables.curGridVolVars();
        else
            return gridVariables.curGridVars();
    };

    auto elemVolVars = localView(curGridVars());
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        fvGeometry.bind(element);
        const auto geometry = fvGeometry.elementGeometry();

        elemVolVars.bind(element, fvGeometry, x);
        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        const auto& quad = Dune::QuadratureRules<double, GridGeometry::GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            const auto& localPos = qp.position();
            const auto& qpVolumeWeight = qp.weight() * Extrusion::integrationElement(geometry, localPos);
            totalVolume += qpVolumeWeight;

            const auto& globalPos = geometry.global(localPos);
            const auto analyticalSolution = problem.analyticalSolution(globalPos);
            const auto numericalSolution = evalSolutionAtLocalPos(element, geometry, gridDiscretization, elemSol, localPos);
            const auto solDiff = numericalSolution - analyticalSolution;
            error += (solDiff * solDiff) * qpVolumeWeight;
        }
    }
    error = std::sqrt(error);

    return {totalVolume, error};
}


/*!
 * \brief Write the mesh size and the discrete l2 and l-infinity errors for the convergence test driver
 * \note L2Abs(v) is the l2 error of the velocity vector
 */
template<class MomentumProblem, class MassProblem, class Scalar>
void writeConvergenceTestErrors(std::shared_ptr<MomentumProblem> momentumProblem,
                                std::shared_ptr<MassProblem> massProblem,
                                const Errors<MomentumProblem, MassProblem, Scalar>& errors)
{
    static constexpr int dim
        = ProblemTraits<MomentumProblem>::GridGeometry::GridView::dimension;
    using MassIndices = typename MassProblem::Indices;
    using MomIndices = typename MomentumProblem::Indices;

    const auto velLabels = velocityLabels();
    std::vector<std::pair<std::string, Scalar>> values;
    const auto add = [&](const auto& e, const std::string& name){
        values.emplace_back(name + "(p)", e[MassIndices::pressureIdx]);
        for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
            values.emplace_back(name + velLabels[dirIdx], e[MomIndices::velocity(dirIdx)+1]);
    };

    add(errors.l2Absolute(), "L2Abs");
    add(errors.l2Relative(), "L2Rel");
    add(errors.lInfAbsolute(), "LinfAbs");
    add(errors.lInfRelative(), "LinfRel");

    using std::sqrt;
    Scalar velocityError = 0.0;
    for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
        velocityError += errors.l2Absolute()[MomIndices::velocity(dirIdx)+1]*errors.l2Absolute()[MomIndices::velocity(dirIdx)+1];
    values.emplace_back("L2Abs(v)", sqrt(velocityError));

    Test::writeConvergenceTestErrors(massProblem->name(), errors.hMax()[0], values);
}

/*!
 * \brief Write the mesh size and the discrete l2 and l-infinity errors of the velocity for the convergence test driver
 */
template<class MomentumProblem, class Scalar>
void writeConvergenceTestErrorsMomentum(std::shared_ptr<MomentumProblem> problem,
                                        const ErrorsSubProblem<MomentumProblem, Scalar>& errors)
{
    static constexpr int dim
        = ProblemTraits<MomentumProblem>::GridGeometry::GridView::dimension;
    using MomIndices = typename MomentumProblem::Indices;

    const auto velLabels = velocityLabels();
    std::vector<std::pair<std::string, Scalar>> values;
    const auto add = [&](const auto& e, const std::string& name){
        for (int dirIdx = 0; dirIdx < dim; ++dirIdx)
            values.emplace_back(name + velLabels[dirIdx], e[MomIndices::velocity(dirIdx)]);
    };

    add(errors.l2Absolute(), "L2Abs");
    add(errors.l2Relative(), "L2Rel");
    add(errors.lInfAbsolute(), "LinfAbs");
    add(errors.lInfRelative(), "LinfRel");

    Test::writeConvergenceTestErrors(problem->name(), errors.hMax(), values);
}

} // end namespace NavierStokesTest

} // end namespace Dumux

#endif
