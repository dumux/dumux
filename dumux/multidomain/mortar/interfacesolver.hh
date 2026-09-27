// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Solver for the interface problem of a mortar-coupled model.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_INTERFACE_SOLVER_HH
#define DUMUX_MULTIDOMAIN_MORTAR_INTERFACE_SOLVER_HH

#include <cstddef>
#include <memory>
#include <string>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/solver.hh>
#include <dune/istl/solvers.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

#include "compatibility.hh"
#include "couplingmode.hh"
#include "interfaceoperator.hh"
#include "preconditioners.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief Solves the interface problem of a mortar-coupled model: the mortar datum for which
 *        the conjugate traces of the subdomains agree, by a Krylov method driven by the
 *        interface operator.
 *
 * The right-hand side is the residual of the subdomains solved with their own data, the
 * operator is the residual of the subdomains solved without it, and after the solve the
 * subdomains are solved once more with their data at the solution. In natural mode the
 * iterate is confined to the data compatible with floating subdomains. The solver is set up
 * for the coupling mode the model is in at construction.
 *
 * Parameters, read from the given group:
 * - `Mortar.Solver`: `cg` (default), `gmres` or `bicgstab`. The interface operator is
 *   symmetric where every subdomain reads back the adjoint of the data it is imposed, as a
 *   cell-centred subdomain does with data per trace cell; otherwise, e.g. with data per trace
 *   vertex read back per trace cell, it is not, which rules out conjugate gradients.
 * - `Mortar.Preconditioner`: `none` (default) or `interface`, the conjugate-mode preconditioner
 *   (see InterfacePreconditioner), developed in \cite Boon2023 for Darcy flow on both sides
 *   of every mortar.
 * - `Mortar.ResidualReduction` (default 1e-8), `Mortar.MaxIterations` (default 1000),
 *   `Mortar.GMResRestart` (default 100) and `Mortar.Verbosity` (default 1).
 */
template<typename M>
class InterfaceSolver
{
    using SolutionVector = typename M::SolutionVector;
    using Scalar = typename SolutionVector::field_type;
    using LinearOperator = Dune::LinearOperator<SolutionVector, SolutionVector>;
    using Preconditioner = Dune::Preconditioner<SolutionVector, SolutionVector>;
    using Constraints = CompatibilityConstraints<M>;

public:
    using Model = M;
    using Operator = InterfaceOperator<Model>;

    /*!
     * \brief Solver for the interface problem of the given model, in the coupling mode the
     *        model is in.
     * \param model The mortar model
     * \param paramGroup The group the `Mortar` parameters are read from
     */
    InterfaceSolver(std::shared_ptr<Model> model, const std::string& paramGroup = "")
    : model_(std::move(model))
    , operator_(model_)
    , solverName_(getParamFromGroup<std::string>(paramGroup, "Mortar.Solver", "cg"))
    , residualReduction_(getParamFromGroup<Scalar>(paramGroup, "Mortar.ResidualReduction", 1e-8))
    , maxIterations_(getParamFromGroup<int>(paramGroup, "Mortar.MaxIterations", 1000))
    , restart_(getParamFromGroup<int>(paramGroup, "Mortar.GMResRestart", 100))
    , verbosity_(getParamFromGroup<int>(paramGroup, "Mortar.Verbosity", 1))
    {
        if (solverName_ != "cg" && solverName_ != "gmres" && solverName_ != "bicgstab")
            DUNE_THROW(ParameterException, "Unknown interface solver '" << solverName_ << "'");

        const auto preconditionerName = getParamFromGroup<std::string>(paramGroup, "Mortar.Preconditioner", "none");
        if (preconditionerName == "interface")
            preconditioner_ = std::make_unique<InterfacePreconditioner<Model>>(model_, model_->couplingMode());
        else if (preconditionerName == "none")
            preconditioner_ = std::make_unique<NoPreconditioner<SolutionVector>>();
        else
            DUNE_THROW(ParameterException, "Unknown interface preconditioner '" << preconditionerName << "'");

        if (model_->couplingMode() == CouplingMode::natural)
        {
            constraints_ = std::make_unique<Constraints>(model_);
            if (!constraints_->empty())
            {
                constrainedOperator_ = std::make_unique<ConstrainedOperator<Operator, Constraints>>(operator_, *constraints_);
                constrainedPreconditioner_ = std::make_unique<ConstrainedPreconditioner<SolutionVector, Constraints>>(*preconditioner_, *constraints_);
            }
        }
    }

    /*!
     * \brief Solve the interface problem starting from the given mortar datum, which holds
     *        the solution on return; the subdomains are then solved with their data at it.
     */
    Dune::InverseOperatorResult solve(SolutionVector& x)
    {
        // the residual with the subdomains' own data at the initial datum is the right-hand
        // side of the homogeneous problem for the update
        SolutionVector rhs(x.size());
        model_->setHomogeneous(false);
        operator_.apply(x, rhs);
        rhs *= -1.0;
        model_->setHomogeneous(true);

        SolutionVector update(x.size());
        update = 0.0;
        Dune::InverseOperatorResult result;
        if (isConstrained())
        {
            constraints_->project(rhs);
            solve_(*constrainedOperator_, *constrainedPreconditioner_, update, rhs, result);
        }
        else
            solve_(operator_, *preconditioner_, update, rhs, result);
        if (!result.converged)
            DUNE_THROW(NumericalProblem, "The interface solver did not converge within " << result.iterations << " iterations");
        x += update;

        model_->setHomogeneous(false);
        residual_.resize(x.size());
        operator_.apply(x, residual_);
        return result;
    }

    //! The interface operator
    const Operator& linearOperator() const { return operator_; }

    //! The preconditioner, the identity if none was chosen
    Preconditioner& preconditioner() { return *preconditioner_; }

    //! Return true if the iterate is confined to the data compatible with floating subdomains
    bool isConstrained() const
    { return constraints_ && !constraints_->empty(); }

    //! The compatibility constraints, defined in natural mode
    const Constraints& constraints() const
    {
        if (!constraints_)
            DUNE_THROW(Dune::InvalidStateException, "Compatibility constraints exist in natural mode only");
        return *constraints_;
    }

    //! The residual with the subdomains' own data at the solution of the last solve
    const SolutionVector& residual() const { return residual_; }

    //! The mortar model
    Model& model() { return *model_; }
    //! The mortar model
    const Model& model() const { return *model_; }

private:
    void solve_(LinearOperator& op, Preconditioner& prec, SolutionVector& x, SolutionVector& b, Dune::InverseOperatorResult& result) const
    {
        // a single unconstrained unknown is solved directly: the operator applied to one is
        // its value, and a Krylov method exhausts its space after that one application
        if (x.size() == 1 && !isConstrained())
        {
            SolutionVector unit(1), image(1);
            unit = 1.0;
            op.apply(unit, image);
            x[0] = b[0];
            x[0] /= image[0][0];
            result.clear();
            result.converged = true;
            result.iterations = 1;
            return;
        }
        if (solverName_ == "cg")
            Dune::CGSolver<SolutionVector>(op, prec, residualReduction_, maxIterations_, verbosity_).apply(x, b, result);
        else if (solverName_ == "gmres")
            Dune::RestartedGMResSolver<SolutionVector>(op, prec, residualReduction_, restart_, maxIterations_, verbosity_).apply(x, b, result);
        else
            Dune::BiCGSTABSolver<SolutionVector>(op, prec, residualReduction_, maxIterations_, verbosity_).apply(x, b, result);
    }

    std::shared_ptr<Model> model_;
    Operator operator_;
    std::string solverName_;
    Scalar residualReduction_;
    int maxIterations_;
    int restart_;
    int verbosity_;
    std::unique_ptr<Preconditioner> preconditioner_;
    std::unique_ptr<Constraints> constraints_;
    std::unique_ptr<ConstrainedOperator<Operator, Constraints>> constrainedOperator_;
    std::unique_ptr<ConstrainedPreconditioner<SolutionVector, Constraints>> constrainedPreconditioner_;
    SolutionVector residual_;
};

} // end namespace Dumux::Mortar

#endif
