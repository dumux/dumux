// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Constrained-pressure-residual AMG-based solver
 */
#ifndef DUMUX_LINEAR_CPR_AMG_SOLVER_HH
#define DUMUX_LINEAR_CPR_AMG_SOLVER_HH

#include <type_traits>
#include <memory>
#include <tuple>

#include <dune/common/parametertree.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/exceptions.hh>

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#include <dumux/common/math.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/parallelmatrixadapter.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/linear/symmetrize_constraints.hh>

namespace Dumux::Detail {

/*!
 * \ingroup Constrained-pressure-residual AMG preconditioner
 *
 * PR-AMG is a multiplicative preconditioner of the form
 * P^-1 = ILU(A)*(I - A*AMG_p(A)) + AMG_p(A)
 * where AMG_p only acts on the pressure unknowns
 *
 * \tparam M Type of the matrix.
 * \tparam X Type of the update.
 * \tparam Y Type of the defect.
 * \tparam l Preconditioner block level
 */
template<class M, class X, class Y, int l = 1>
class CPRAMGPreconditioner : public Dune::Preconditioner<X,Y>
{
    using PressureMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using PressureVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using PressureLOP = Dune::MatrixAdapter<PressureMatrix, PressureVector, PressureVector>;
public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = M;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief Scalar type underlying the field_type.
    using scalar_field_type = Dune::Simd::Scalar<field_type>;

    CPRAMGPreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& op,
        const Dune::ParameterTree& params
    )
    : op_(op)
    {
        pressureIndex_ = getParamFromGroup<int>(params.template get<std::string>("ParameterGroup"), "PressureIndex", 0);

        // copy sparsity pattern
        Dune::MatrixIndexSet sparsityPattern(op_->getmat().N(), op_->getmat().M());
        sparsityPattern.import(op_->getmat());
        sparsityPattern.exportIdx(pressureMatrix_());

        // initialize pressure matrix
        for (auto rowIt = op_->getmat().begin(); rowIt != op_->getmat().end(); ++rowIt)
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                pressureMatrix_()[rowIt.index()][colIt.index()][0][0] = (*colIt)[pressureIndex_][pressureIndex_];

        // build AMG_p(A) with SSOR smoother
        amg_ = std::make_shared<Dune::Amg::AMG<
            PressureLOP, PressureVector,
            Dune::SeqSOR<PressureMatrix, PressureVector, PressureVector>
        >>(std::make_shared<PressureLOP>(pressureMatrix_()), params);

        // build ILU(A)
        ilu_ = std::make_shared<Dune::SeqILU<M, X, Y>>(op_->getmat(), 1, 1.0);
    }

    /*!
     * \brief Prepare the preconditioner.
     */
    void pre(X& update, Y& currentDefect) override {}

    /*!
     * \brief Apply the preconditioner
     *
     * \param update The update to be computed.
     * \param currentDefect The current defect.
     *
     * The currentDefect has be be in a consistent representation,
     * Definition 2.3 Blatt and Bastian (2009) https://doi.org/10.1504/IJCSE.2008.021112
     * The update is initially zero. At exit the update has to be
     * in a consistent representation. This usually requires communication.
     */
    void apply(X& update, const Y& currentDefect) override
    {
        update = 0.0;

        PressureVector pDefect(currentDefect.size());
        for (std::size_t i = 0; i < currentDefect.size(); ++i)
            pDefect[i][0] = currentDefect[i][pressureIndex_];

        // v = AMG_p(pDefect)
        auto v = update;
        {
            PressureVector pUpdate(update.size());
            applyAMG_(pUpdate, pDefect);
            for (std::size_t i = 0; i < pUpdate.size(); ++i)
                v[i][pressureIndex_] = pUpdate[i][0];
        }

        // w = defect - Av;
        auto w = update;
        op_->apply(v, w);
        w *= -1;
        w += currentDefect;

        // u = v + ILU(w)
        applyILU_(update, w);
        update += v;
    }

    /*!
     * \brief Clean up.
     */
    void post(X& update) override {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    template<class Sol, class Rhs>
    void applyAMG_(Sol& sol, Rhs& rhs) const
    {
        amg_->pre(sol, rhs);
        amg_->apply(sol, rhs);
        amg_->post(sol);
    }

    template<class Sol, class Rhs>
    void applyILU_(Sol& sol, Rhs& rhs) const
    {
        ilu_->pre(sol, rhs);
        ilu_->apply(sol, rhs);
        ilu_->post(sol);
    }

    static PressureMatrix& pressureMatrix_() {
        static PressureMatrix p;
        return p;
    }

    std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>> op_;
    int pressureIndex_;
    //PressureMatrix pressureMatrix_;
    std::shared_ptr<Dune::Preconditioner<PressureVector, PressureVector>> amg_;
    std::shared_ptr<Dune::Preconditioner<X, Y>> ilu_;
};

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief An CPR-AMG preconditioned BiCGSTAB solver using dune-istl
 */
template<class LSTraits, class LATraits>
using CPRAMGBiCGSTABIstlSolver =
    Detail::IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::BiCGSTABSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlSolvers::IstlDefaultBlockLevelPreconditionerFactory<Detail::CPRAMGPreconditioner>,
        // the Dune::ILU preconditioners don't accept multi-type matrices
        /*convert multi-type istl types?*/ true
    >;

} // end namespace Dumux

#endif
