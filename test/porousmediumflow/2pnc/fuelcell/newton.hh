// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_LOGARITHM_TRANSFORMATION_NEWTON_SOLVER_HH
#define DUMUX_LOGARITHM_TRANSFORMATION_NEWTON_SOLVER_HH

#include <algorithm>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

template <class Assembler, class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator> >
class LogarithmTransformationNewtonSolver : public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using Scalar = typename Assembler::Scalar;
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    using Indices = typename Assembler::GridVariables::VolumeVariables::Indices;

    using typename ParentType::Backend;
    using typename ParentType::SolutionVector;
    using typename ParentType::ResidualVector;

public:
    using ParentType::ParentType;
    using typename ParentType::Variables;

private:

    void newtonBeginStep(const Variables& vars) override
    {
        ParentType::newtonBeginStep(vars);
        uLastIter_ = Backend::dofs(vars);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     * Throws Dumux::NumericalProblem if the linear solver didn't converge.
     */
    bool solveLinearSystem_(ResidualVector& deltaU) override
    {
        static const bool useLogTranformation = getParam<bool>("Newton.UseLogTransformation", false);
        if (!useLogTranformation)
            return this->linearSolver().solve(
                this->assembler().jacobian(), deltaU, this->assembler().residual()
            );

        // Solve the linear system in a modified way
        auto JLogC = this->assembler().jacobian(); // copy
        auto deltaLogC = deltaU; // copy

        // modify the Jacobian by multiplying each molefraction derivative by the molefraction itself
        // JLogC = dR/dlnC = J * dC/dlnC
        // dx/dlnx = (dlnx/dx)^-1 = (1/x)^-1 = x
        for (auto rowIt = JLogC.begin(); rowIt != JLogC.end(); ++rowIt)
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                for (int i = 0; i < colIt->M(); ++i)
                    for (int j = 2; j < colIt->N(); ++j)
                        (*colIt)[i][j] *= uLastIter_[colIt.index()][j];

        bool converged = this->linearSolver().solve(
            JLogC,
            deltaLogC,
            this->assembler().residual()
        );

        if (!converged)
            DUNE_THROW(Dumux::NumericalProblem, "Linear solver did not converge.");

        // Transform the solution back to the original space
        for (int i = 0; i < deltaU.size(); ++i)
        {
            deltaU[i] = deltaLogC[i]; // copy back

            // overwrite the solution with the transformed solution for the mole fractions
            // lnCNew = lnCOld - deltaLnC
            // => CNew = COld * exp(-deltaLnC) = COld - deltaC
            // => deltaC = COld - COld * exp(-deltaLnC)
            for (int j = 2; j < deltaU[i].size(); ++j)
                deltaU[i][j] = uLastIter_[i][j] -  uLastIter_[i][j]*std::exp(-deltaLogC[i][j]);
        }

        return converged;
    }
private:
    SolutionVector uLastIter_;

};

} // end namespace Dumux

#endif
