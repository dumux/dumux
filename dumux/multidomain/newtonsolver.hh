// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Nonlinear
 * \copydoc Dumux::MultiDomainNewtonSolver
 *
 */
#ifndef DUMUX_MULTIDOMAIN_NEWTON_SOLVER_HH
#define DUMUX_MULTIDOMAIN_NEWTON_SOLVER_HH

#include <memory>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief Newton sover for coupled problems
 */
template <class Assembler, class LinearSolver, class CouplingManager,
          class Reassembler = DefaultPartialReassembler,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class MultiDomainNewtonSolver: public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    using Scalar = typename Assembler::Scalar;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using SolutionVector = typename Assembler::ResidualType;
    using ConvergenceWriter = ConvergenceWriterInterface<SolutionVector>;

public:

    /*!
     * \brief Constructor for stationary problems
     */
    MultiDomainNewtonSolver(std::shared_ptr<Assembler> assembler,
                            std::shared_ptr<LinearSolver> linearSolver,
                            std::shared_ptr<CouplingManager> couplingManager,
                            const Comm& comm = Dune::MPIHelper::getCollectiveCommunication(),
                            const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver, comm, paramGroup)
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief Constructor for instationary problems
     */
    MultiDomainNewtonSolver(std::shared_ptr<Assembler> assembler,
                            std::shared_ptr<LinearSolver> linearSolver,
                            std::shared_ptr<CouplingManager> couplingManager,
                            std::shared_ptr<TimeLoop<Scalar>> timeLoop,
                            const Comm& comm = Dune::MPIHelper::getCollectiveCommunication(),
                            const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver, timeLoop, comm, paramGroup)
    , couplingManager_(couplingManager)
    {}


    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void newtonBeginStep(const SolutionVector& uCurrentIter) override
    {
        ParentType::newtonBeginStep(uCurrentIter);
        couplingManager_->updateSolution(uCurrentIter);
    }


    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param assembler The jacobian assembler
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void newtonEndStep(SolutionVector &uCurrentIter, const SolutionVector &uLastIter) override
    {
        ParentType::newtonEndStep(uCurrentIter, uLastIter);
        couplingManager_->updateSolution(uCurrentIter);
    }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
