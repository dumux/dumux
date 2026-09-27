// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Preconditioners for mortar-coupling models.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_PRECONDITIONERS_HH
#define DUMUX_MULTIDOMAIN_MORTAR_PRECONDITIONERS_HH

#include <algorithm>
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dumux/common/exceptions.hh>

#include "couplingmode.hh"
#include "projectors.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief The identity as a preconditioner of the interface operator.
 */
template<typename X>
struct NoPreconditioner : public Dune::Preconditioner<X, X>
{
    void pre (X&, X&) override {}
    //! Return the defect x unchanged as the update r
    void apply (X& r, const X& x) override { r = x; }
    void post (X&) override {}

    //! The preconditioner acts on the global mortar vector without communication
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }
};

/*!
 * \ingroup MortarCoupling
 * \brief Preconditions the interface operator with subdomain solves in the conjugate coupling
 *        mode: an operator built from essential (Dirichlet) solves is preconditioned by natural
 *        (Neumann) solves and vice versa, so the preconditioner approximates the inverse of the
 *        sum of the local Steklov-Poincaré operators the interface operator is built from. For
 *        the flux-mortar variant this is the interface preconditioner of \cite Boon2023; for a
 *        value mortar it is its conjugate, of Neumann-Neumann type.
 *
 * One application computes \f$ N\,d = M_\Lambda^{-1}\big(\sum_i P_i \Theta_i P_i^{\mathsf{T}}\big)
 * M_\Lambda^{-1} d \f$, where \f$M_\Lambda\f$ is the mortar mass matrix, \f$P_i\f$ the coupling
 * matrices, and \f$\Theta_i\f$ is realized without assembly by one subdomain solve in the
 * conjugate mode. The defect handed to a preconditioner is a functional on the mortar space,
 * while imposed trace data and the returned update are functions in it; the two mass-matrix
 * solves are the Riesz maps that mediate between the two. They are solved with the sparse
 * mortar mass matrix by conjugate gradients to a residual reduction of \f$10^{-13}\f$.
 *
 * A subdomain whose entire boundary couples to mortars has a singular conjugate (Neumann)
 * problem, defined only for data with zero total flux and only up to a constant. For each such
 * subdomain a coarse function \f$\varphi_i = M_\Lambda^{-1}c_i\f$ is added, with \f$(c_i)_k\f$
 * the integral of mortar basis function \f$k\f$ over that subdomain's boundary, and the
 * preconditioner becomes the balanced form of \cite Mandel1993,
 * \f$ B = P_0 + (I - P_0 S)\,N\,(I - S P_0) \f$ with the coarse projection
 * \f$ P_0 = \Phi(\Phi^{\mathsf{T}}S\Phi)^{-1}\Phi^{\mathsf{T}} \f$: the pre-balancing makes
 * every singular subdomain's data compatible, the post-balancing removes the constant left
 * undetermined by its solve. \f$S\Phi\f$ is computed once with one operator sweep per coarse
 * function; each application then still costs a single subdomain sweep. The singular subdomain
 * solves themselves must be made unique by the problem, e.g. by an internal Dirichlet constraint
 * pinning one degree of freedom.
 *
 * \note Applications assume the model is in its homogeneous configuration, as it is inside a
 *       Krylov solver driven by the interface operator.
 */
template<typename M>
class InterfacePreconditioner
: public Dune::Preconditioner<typename M::SolutionVector, typename M::SolutionVector>
{
    using SolutionVector = typename M::SolutionVector;
    using Scalar = typename SolutionVector::field_type;
    using MassMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>>;
    using BlockVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;

    struct MortarBlock
    {
        std::size_t id;
        std::size_t offset;
        MassMatrix mass;
        std::vector<Scalar> basisIntegrals;
    };

    struct ModeGuard
    {
        M& model;
        CouplingMode restoreTo;

        // a destructor must not throw, so a failing restore during unwinding is swallowed
        ~ModeGuard()
        {
            try { model.setCouplingMode(restoreTo); }
            catch (...) {}
        }
    };

public:
    using Model = M;

    /*!
     * \brief Preconditioner for the interface operator of the given model.
     * \param model The mortar model the interface operator is built from
     * \param operatorMode The coupling mode of the interface operator; the subdomain solves
     *        of the preconditioner run in the conjugate mode
     */
    explicit InterfacePreconditioner(std::shared_ptr<Model> model,
                                     CouplingMode operatorMode = CouplingMode::essential)
    : model_{std::move(model)}
    , operatorMode_{operatorMode}
    {
        model_->decomposition().visitMortars([&] (const auto& mortarPtr) {
            blocks_.push_back(assembleMassData_(*mortarPtr));
        });
        buildCoarseSpace_();
    }

    void pre (SolutionVector&, SolutionVector&) override {}
    void post (SolutionVector&) override {}

    //! Apply the preconditioner to the defect d, the result is the update v
    void apply (SolutionVector& v, const SolutionVector& d) override
    {
        ensureCoarseSetup_();
        if (phi_.empty())
        {
            v = conjugateSweep_(d);
            return;
        }

        const auto n = phi_.size();
        Dune::DynamicVector<Scalar> rhs(n), beta0(n), beta1(n);
        for (std::size_t i = 0; i < n; ++i)
            rhs[i] = phi_[i]*d;
        gramInverse_.mv(rhs, beta0);

        auto balanced = d;
        for (std::size_t i = 0; i < n; ++i)
            balanced.axpy(-beta0[i], sPhi_[i]);

        const auto w = conjugateSweep_(balanced);

        for (std::size_t i = 0; i < n; ++i)
            rhs[i] = sPhi_[i]*w;
        gramInverse_.mv(rhs, beta1);

        v = w;
        for (std::size_t i = 0; i < n; ++i)
            v.axpy(beta0[i] - beta1[i], phi_[i]);
    }

    //! The preconditioner acts on the global mortar vector without communication
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    static CouplingMode conjugateOf_(CouplingMode mode)
    {
        return mode == CouplingMode::essential ? CouplingMode::natural
                                               : CouplingMode::essential;
    }

    //! One conjugate-mode sweep between the two Riesz maps
    SolutionVector conjugateSweep_(const SolutionVector& d)
    {
        const auto z = rieszRepresentative_(d);
        ModeGuard guard{*model_, operatorMode_};
        model_->setCouplingMode(conjugateOf_(operatorMode_));
        model_->setMortar(z);
        model_->solveSubDomains();
        SolutionVector w(d.size());
        w = 0.0;
        model_->assembleMortarResidual(w);
        return rieszRepresentative_(w);
    }

    //! One sweep of the interface operator itself, in the operator's own mode
    SolutionVector operatorSweep_(const SolutionVector& x)
    {
        model_->setMortar(x);
        model_->solveSubDomains();
        SolutionVector w(x.size());
        w = 0.0;
        model_->assembleMortarResidual(w);
        return w;
    }

    //! Solve with the mortar mass matrix, turning a functional into a function
    SolutionVector rieszRepresentative_(const SolutionVector& d) const
    {
        SolutionVector result(d.size());
        result = 0.0;
        for (const auto& block : blocks_)
        {
            const auto n = block.mass.N();
            BlockVector rhs(n), representative(n);
            for (std::size_t i = 0; i < n; ++i)
                rhs[i] = d[block.offset + i];
            representative = 0.0;

            Dune::MatrixAdapter<MassMatrix, BlockVector, BlockVector> massOperator(block.mass);
            Dune::SeqSSOR<MassMatrix, BlockVector, BlockVector> ssor(block.mass, 1, 1.0);
            Dune::CGSolver<BlockVector> cg(massOperator, ssor, massSolverReduction_, massSolverMaxIterations_, 0);
            Dune::InverseOperatorResult solveResult;
            cg.apply(representative, rhs, solveResult);
            if (!solveResult.converged)
                DUNE_THROW(NumericalProblem, "The mortar mass matrix solve did not converge within "
                           << solveResult.iterations << " iterations");

            for (std::size_t i = 0; i < n; ++i)
                result[block.offset + i] = representative[i];
        }
        return result;
    }

    template<typename MortarGridGeometry>
    MortarBlock assembleMassData_(const MortarGridGeometry& mortar) const
    {
        auto mass = Detail::assembleMortarMassMatrix<Scalar>(mortar);
        auto integrals = Detail::mortarBasisIntegrals<Scalar>(mortar);
        return {
            model_->decomposition().id(mortar),
            model_->mortarDofOffset(mortar),
            std::move(mass),
            std::move(integrals)
        };
    }

    /*!
     * \brief Builds the coarse space: per floating subdomain the Riesz representative of its
     *        boundary-indicator moments, carrying the signs its natural-mode data carries,
     *        and, when any subdomain floats, additionally the mortar basis functions at the
     *        boundary vertices of every mortar, the interface endpoints at cross points, where
     *        an equal-weight sum of subdomain solves represents the operator worst.
     */
    void buildCoarseSpace_()
    {
        std::unordered_map<std::size_t, SolutionVector> columns;
        model_->visitCouplings([&] (const auto& mortar, const auto& solver, std::size_t subDomainId) {
            if (!solver.isFloating())
                return;

            auto& column = columns[subDomainId];
            if (column.size() != model_->numMortarDofs())
            {
                column.resize(model_->numMortarDofs());
                column = 0.0;
            }

            const auto mortarId = model_->decomposition().id(mortar);
            const auto it = std::find_if(blocks_.begin(), blocks_.end(),
                                         [&] (const auto& b) { return b.id == mortarId; });
            const auto sign = static_cast<Scalar>(model_->orientation(subDomainId, mortarId));
            for (std::size_t k = 0; k < it->basisIntegrals.size(); ++k)
                column[it->offset + k] += sign*it->basisIntegrals[k];
        });

        if (columns.empty())
            return;

        for (auto& [id, column] : columns)
            phi_.push_back(rieszRepresentative_(column));

        model_->decomposition().visitMortars([&] (const auto& mortarPtr) {
            using MortarGridGeometry = std::remove_cvref_t<decltype(*mortarPtr)>;
            if constexpr (Detail::mortarSpaceOrder<MortarGridGeometry> == 1)
            {
                const auto offset = model_->mortarDofOffset(*mortarPtr);
                for (const auto dof : boundaryVertexDofs_(*mortarPtr))
                {
                    SolutionVector column(model_->numMortarDofs());
                    column = 0.0;
                    column[offset + dof] = 1.0;
                    phi_.push_back(std::move(column));
                }
            }
        });
    }

    template<typename MortarGridGeometry>
    std::vector<std::size_t> boundaryVertexDofs_(const MortarGridGeometry& mortar) const
    {
        static constexpr int dim = MortarGridGeometry::GridView::dimension;
        std::vector<std::size_t> dofs;
        for (const auto& element : elements(mortar.gridView()))
            for (const auto& is : intersections(mortar.gridView(), element))
            {
                if (!is.boundary())
                    continue;
                const auto refElement = Dune::referenceElement<Scalar, dim>(element.type());
                for (int i = 0; i < refElement.size(is.indexInInside(), 1, dim); ++i)
                    dofs.push_back(mortar.vertexMapper().subIndex(
                        element, refElement.subEntity(is.indexInInside(), 1, i, dim), dim
                    ));
            }
        std::sort(dofs.begin(), dofs.end());
        dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
        return dofs;
    }

    //! Computes the coarse operator images on first application, when the model is homogeneous
    void ensureCoarseSetup_()
    {
        if (coarseSetupDone_ || phi_.empty())
        {
            coarseSetupDone_ = true;
            return;
        }

        const auto n = phi_.size();
        sPhi_.reserve(n);
        for (const auto& phi : phi_)
            sPhi_.push_back(operatorSweep_(phi));

        gramInverse_.resize(n, n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                gramInverse_[i][j] = phi_[i]*sPhi_[j];
        gramInverse_.invert();
        coarseSetupDone_ = true;
    }

    static constexpr Scalar massSolverReduction_ = 1e-13;
    static constexpr int massSolverMaxIterations_ = 100;

    std::shared_ptr<Model> model_;
    CouplingMode operatorMode_;
    std::vector<MortarBlock> blocks_;

    std::vector<SolutionVector> phi_;
    std::vector<SolutionVector> sPhi_;
    Dune::DynamicMatrix<Scalar> gramInverse_;
    bool coarseSetupDone_ = false;
};

} // end namespace Dumux::Mortar

#endif
