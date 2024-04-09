// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Reference implementation of a Newton solver.
 */
#ifndef DUMUX_NEWTON_SOLVER_DROPLET_HH
#define DUMUX_NEWTON_SOLVER_DROPLET_HH

#include <cmath>
#include <memory>
#include <iostream>
#include <type_traits>
#include <algorithm>
#include <numeric>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/pdesolver.hh>
#include <dumux/common/variablesbackend.hh>

#include <dumux/io/format.hh>

#include <dumux/linear/matrixconverter.hh>
#include <dumux/assembly/partialreassembler.hh>

#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief An implementation of a Newton solver
 * \tparam Assembler the assembler
 * \tparam LinearSolver the linear solver
 * \tparam Comm the communication object used to communicate with all processes
 * \note If you want to specialize only some methods but are happy with the
 *       defaults of the reference solver, derive your solver from
 *       this class and simply overload the required methods.
 */
template <class Assembler, class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator> >
class NewtonSolverDroplet : public NewtonSolver<Assembler, LinearSolver>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver>;

protected:
    using Backend = VariablesBackend<typename ParentType::Variables>;
    using SolutionVector = typename Backend::DofVector;
    using ResidualVector = typename Assembler::ResidualType;
    using LinearAlgebraNativeBackend = VariablesBackend<ResidualVector>;
private:
    using Scalar = typename Assembler::Scalar;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using ConvergenceWriter = ConvergenceWriterInterface<SolutionVector, ResidualVector>;
    using TimeLoop = TimeLoopBase<Scalar>;

    // enable models with primary variable switch
    // TODO: Always use ParentType::Variables once we require assemblers to export variables
    static constexpr bool assemblerExportsVariables = Detail::PDESolver::assemblerExportsVariables<Assembler>;
    using PriVarSwitchVariables
        = std::conditional_t<assemblerExportsVariables,
                             typename ParentType::Variables,
                             Detail::Newton::PriVarSwitchVariables<Assembler>>;
    using PrimaryVariableSwitchAdapter = Dumux::PrimaryVariableSwitchAdapter<PriVarSwitchVariables>;

public:
    using typename ParentType::Variables;
    using Communication = Comm;

    /*!
     * \brief The Constructor
     */
    NewtonSolverDroplet(std::shared_ptr<Assembler> assembler,
                 std::shared_ptr<LinearSolver> linearSolver,
                 const Communication& comm = Dune::MPIHelper::getCommunication(),
                 const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver, comm, paramGroup)
    {}

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const override
    {
        const auto dropSolver = this->assembler().problem().dropSolver();

        if (!(dropSolver->checkDropletsGeometry()))
            return false;

        return ParentType::newtonConverged();
    }

};

} // end namespace Dumux

#endif
