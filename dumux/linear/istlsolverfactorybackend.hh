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
 * \ingroup Linear
 * \brief Provides a generic linear solver based on the ISTL that chooses the
 *        solver and preconditioner at runtime
 */

#ifndef DUMUX_LINEAR_ISTL_SOLVERFACTORYBACKEND_HH
#define DUMUX_LINEAR_ISTL_SOLVERFACTORYBACKEND_HH

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>

// preconditioners
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>

// solvers
#include <dune/istl/solvers.hh>
#include <dune/istl/solverfactory.hh>

#include <dumux/linear/solver.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/parallelhelpers.hh>

#if DUNE_VERSION_NEWER_REV(DUNE_ISTL,2,7,1)

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A linear solver using the dune-istl solver factory
 *        to choose the solver and preconditioner at runtime.
 * \note the solvers are configured via the input file
 * \note requires Dune version 2.7.1 or newer
 */
template <class Matrix, class Vector, class GridGeometry>
class IstlSolverFactoryBackend : public LinearSolver
{
    using GridView = typename GridGeometry::GridView;
    using LinearSolverTraits =  Dumux::LinearSolverTraits<Matrix, Vector, GridGeometry>;
    using Grid = typename GridView::Grid;
    using LinearOperator = typename LinearSolverTraits::LinearOperator;
    using ScalarProduct = typename LinearSolverTraits::ScalarProduct;
    using VType = typename LinearSolverTraits::VType;
    using Comm = typename LinearSolverTraits::Comm;
    using BCRSMat = typename LinearSolverTraits::LinearOperator::matrix_type;
    using DofMapper = typename LinearSolverTraits::DofMapper;

public:
    //! translation table for solver parameters
    static std::vector<std::array<std::string,2> > dumuxToIstlSolverParams;

    /*!
     * \brief Construct the backend for the sequential case only
     *
     * \param paramGroup the parameter group for parameter lookup
     */
    IstlSolverFactoryBackend(const std::string& paramGroup = "")
    : paramGroup_(paramGroup)
    , firstCall_(true)
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
            DUNE_THROW(Dune::InvalidStateException, "Using sequential constructor for parallel run. Use signature with gridView and dofMapper!");

        reset();
    }

    /*!
     * \brief Construct the backend for parallel or sequential runs
     *
     * \param gridView the grid view for parallel communication via the grid
     * \param dofMapper an index mapper for dof entities
     * \param paramGroup the parameter group for parameter lookup
     */
    IstlSolverFactoryBackend(const GridView& gridView,
                             const DofMapper& dofMapper,
                             const std::string& paramGroup = "")
    : paramGroup_(paramGroup)
    , parallelHelper_(std::make_unique<ParallelISTLHelper<GridView, LinearSolverTraits>>(gridView, dofMapper))
    , firstCall_(true)
    {
        reset();
    }

    /*!
     * \brief Solve a linear system.
     *
     * \param A the matrix
     * \param x the seeked solution vector, containing the initial solution upon entry
     * \param b the right hand side vector
     */
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        std::shared_ptr<Comm> comm;
        std::shared_ptr<LinearOperator> linearOperator;
        std::shared_ptr<ScalarProduct> scalarProduct; // not used.

#if HAVE_MPI
        if constexpr (LinearSolverTraits::isParallel)
            prepareLinearAlgebraParallel<LinearSolverTraits>(A, b, comm, linearOperator, scalarProduct, *parallelHelper_, firstCall_);
        else
            prepareLinearAlgebraSequential<LinearSolverTraits>(A, comm, linearOperator, scalarProduct);
#else
        prepareLinearAlgebraSequential<LinearSolverTraits>(A, comm, linearOperator, scalarProduct);
#endif

        // construct solver
        auto solver = getSolverFromFactory_(linearOperator);

        // solve linear system
        solver->apply(x, b, result_);

        firstCall_ = false;
        return result_.converged;
    }

    //! reset the linear solver factory
    void reset()
    {
        resetDefaultParameters();
        convertParameterTree_(paramGroup_);
        checkMandatoryParameters_();
        Dune::initSolverFactories<typename LinearSolverTraits::LinearOperator>();
        name_ = params_.get<std::string>("preconditioner.type") + "-preconditioned " + params_.get<std::string>("type");
        if (params_.get<int>("verbose", 0) > 0)
            std::cout << "Initialized linear solver of type: " << name_ << std::endl;
    }

    //! reset some defaults for the solver parameters
    void resetDefaultParameters()
    {
        params_["restart"] = "10";
        params_["maxit"] = "250";
        params_["reduction"] = "1e-13";
        params_["verbose"] = "0";
        params_["preconditioner.iterations"] = "1";
        params_["preconditioner.relaxation"] = "1.0";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

    const std::string& name() const
    {
        return name_;
    }

private:

    void convertParameterTree_(const std::string& paramGroup="")
    {
        auto linearSolverGroups = getParamSubGroups("LinearSolver", paramGroup);
        if (linearSolverGroups.empty()) // no linear solver parameters were specified
            return;

        for (const auto& [dumuxKey, istlKey] : dumuxToIstlSolverParams)
        {
            for (const auto& group : linearSolverGroups)
            {
                const auto fullDumuxKey = group + "." + dumuxKey;
                const auto value = getParam<std::string>(fullDumuxKey, "");
                if (!value.empty())
                {
                    params_[istlKey] = value;
                    break; // skip groups with smaller depth in the tree
                }
            }
        }
    }

    void checkMandatoryParameters_()
    {
        if (!params_.hasKey("type"))
            DUNE_THROW(Dune::InvalidStateException, "Solver factory needs \"LinearSolver.Type\" parameter to select the solver");

        if (!params_.hasKey("preconditioner.type"))
            DUNE_THROW(Dune::InvalidStateException, "Solver factory needs \"LinearSolver.Preconditioner.Type\" parameter to select the preconditioner");
    }

    auto getSolverFromFactory_(std::shared_ptr<LinearOperator>& fop)
    {
        try { return Dune::getSolverFromFactory(fop, params_); }
        catch(Dune::Exception& e)
        {
            std::cerr << "Could not create solver with factory" << std::endl;
            std::cerr << e.what() << std::endl;
            throw e;
        }
    }

    const std::string paramGroup_;
    std::unique_ptr<ParallelISTLHelper<GridView, LinearSolverTraits>> parallelHelper_;
    bool firstCall_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::string name_;
};

//! translation table for solver parameters
template<class Matrix, class Vector, class Geometry>
std::vector<std::array<std::string, 2>>
IstlSolverFactoryBackend<Matrix, Vector, Geometry>::dumuxToIstlSolverParams =
{
    // solver params
    {"Verbosity", "verbose"},
    {"MaxIterations", "maxit"},
    {"ResidualReduction", "reduction"},
    {"Type", "type"},
    {"GMResRestart", "restart"}, // cycles before restarting
    {"Restart", "restart"}, // cycles before restarting
    {"MaxOrthogonalizationVectors", "mmax"},

    // preconditioner params
    // we expect parameters in the subgroup "Preconditioner"
    // there are some legacy param names that are still supported
    // subgroup parameters overwrite legacy parameters
    {"PreconditionerVerbosity", "preconditioner.verbosity"},
    {"Preconditioner.Verbosity", "preconditioner.verbosity"},
    {"PreconditionerType", "preconditioner.type"},
    {"Preconditioner.Type", "preconditioner.type"},
    {"PreconditionerIterations", "preconditioner.iterations"},
    {"Preconditioner.Iterations", "preconditioner.iterations"},
    {"PreconditionerRelaxation", "preconditioner.relaxation"},
    {"Preconditioner.Relaxation", "preconditioner.relaxation"},
    {"Preconditioner.ILUOrder", "preconditioner.n"},
    {"Preconditioner.ILUResort", "preconditioner.resort"},
    {"Preconditioner.AmgSmootherRelaxation", "preconditioner.smootherRelaxation"},
    {"Preconditioner.AmgSmootherIterations", "preconditioner.smootherIterations"},
    {"Preconditioner.AmgMaxLevel", "preconditioner.maxLevel"},
    {"Preconditioner.AmgCoarsenTarget", "preconditioner.coarsenTarget"},
    {"Preconditioner.AmgMinCoarseningRate", "preconditioner.minCoarseningRate"},
    {"Preconditioner.AmgProlongationDampingFactor", "preconditioner.prolongationDampingFactor"},
    {"Preconditioner.AmgAlpha", "preconditioner.alpha"},
    {"Preconditioner.AmgBeta", "preconditioner.beta"},
    {"Preconditioner.AmgAdditive", "preconditioner.additive"},
    {"Preconditioner.AmgGamma", "preconditioner.gamma"},
    {"Preconditioner.AmgPreSmoothingSteps", "preconditioner.preSteps"},
    {"Preconditioner.AmgPostSmoothingSteps", "preconditioner.postSteps"},
    {"Preconditioner.AmgCriterionSymmetric", "preconditioner.criterionSymmetric"},
    {"Preconditioner.AmgStrengthMeasure", "preconditioner.strengthMeasure"},
    {"Preconditioner.AmgDiagonalRowIndex", "preconditioner.diagonalRowIndex"},
    {"Preconditioner.AmgDefaultAggregationSizeMode", "preconditioner.defaultAggregationSizeMode"},
    {"Preconditioner.AmgDefaultAggregationDimension", "preconditioner.defaultAggregationDimension"},
    {"Preconditioner.AmgMaxAggregateDistance", "preconditioner.maxAggregateDistance"},
    {"Preconditioner.AmgMinAggregateSize", "preconditioner.minAggregateSize"},
    {"Preconditioner.AmgMaxAggregateSize", "preconditioner.maxAggregateSize"}
};

} // end namespace Dumux

#else
#warning "Generic dune-istl solver factory backend needs dune-istl >= 2.7.1!"
#endif // DUNE version check
#endif // header guard
