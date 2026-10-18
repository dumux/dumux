// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Contact
 * \brief An assembler adding contact contributions to the system of another assembler
 */
#ifndef DUMUX_SOLIDMECHANICS_CONTACT_ASSEMBLER_HH
#define DUMUX_SOLIDMECHANICS_CONTACT_ASSEMBLER_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/jacobianpattern.hh>

namespace Dumux::Contact {

/*!
 * \ingroup Contact
 * \brief Wraps an assembler and adds the residual and Jacobian of a contact model
 *
 * The residual of the wrapped assembler, including its constraints, is extended by the contact
 * residual in all equations that are not constrained. With coupled pairs, the Jacobian is a
 * separate matrix whose pattern is the pattern of the grid discretization plus the entries
 * between each surface node and the corners of its candidate facets; the wrapped Jacobian is
 * copied into it after each assembly. Without coupled pairs, the contact Jacobian is added to
 * the wrapped Jacobian only where its pattern has entries, so the linearization is inexact while
 * the residual stays exact.
 *
 * The contact candidates are kept until updateCandidates or addMissingCandidates is called; both
 * rebuild the pattern.
 *
 * \tparam Assembler an assembler with element-local residual and Jacobian
 * \tparam ContactModel provides updateCandidates, addMissingCandidates, addResidual,
 *         addJacobian and jacobianPattern (e.g. PenaltyContact)
 */
template<class Assembler, class ContactModel>
class ContactAssembler
{
public:
    using Scalar = typename Assembler::Scalar;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using SolutionVector = typename Assembler::SolutionVector;
    using ResidualType = typename Assembler::ResidualType;
    using GridVariables = typename Assembler::GridVariables;
    using GridGeometry = typename Assembler::GridGeometry;
    using Problem = typename Assembler::Problem;

    ContactAssembler(std::shared_ptr<Assembler> assembler,
                     std::shared_ptr<ContactModel> contact,
                     bool couplePairs = true)
    : assembler_(assembler)
    , contact_(contact)
    , couplePairs_(couplePairs)
    {}

    void setLinearSystem()
    {
        assembler_->setLinearSystem();
        residual_ = std::make_shared<ResidualType>(assembler_->numDofs());
        if (couplePairs_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern_();
        }
    }

    //! Search the contact candidates for the solution u
    void updateCandidates(const SolutionVector& u)
    {
        contact_->updateCandidates(u);
        if (jacobian_)
            setJacobianPattern_();
    }

    //! Add the candidates in contact range at u that were missing; returns whether there were any
    bool addMissingCandidates(const SolutionVector& u)
    {
        const bool added = contact_->addMissingCandidates(u);
        if (added && jacobian_)
            setJacobianPattern_();
        return added;
    }

    //! Whether the candidates searched at uOld cover all contacts at u
    bool candidatesCover(const SolutionVector& uOld, const SolutionVector& u) const
    { return contact_->candidatesCover(uOld, u); }

    void assembleJacobianAndResidual(const SolutionVector& u)
    {
        assembler_->assembleJacobianAndResidual(u);
        if (couplePairs_)
            copyEntries_(assembler_->jacobian(), *jacobian_);
        *residual_ = assembler_->residual();

        const auto constrained = constrainedEquations_();
        addContactResidual_(*residual_, u, constrained);

        auto& A = jacobian();
        contact_->addJacobian(u, [&](std::size_t row, std::size_t col, const auto& block)
        {
            auto& matrixRow = A[row];
            const auto entry = matrixRow.find(col);
            if (entry == matrixRow.end())
            {
                if (couplePairs_)
                    DUNE_THROW(Dune::InvalidStateException, "Contact Jacobian entry outside the pattern");
                return;
            }
            for (int i = 0; i < block.N(); ++i)
                if (!constrained[row*blockSize + i])
                    for (int k = 0; k < block.M(); ++k)
                        (*entry)[i][k] += block[i][k];
        });
    }

    void assembleResidual(const SolutionVector& u)
    {
        *residual_ = 0.0;
        assembleResidual(*residual_, u);
    }

    void assembleResidual(ResidualType& r, const SolutionVector& u) const
    {
        assembler_->assembleResidual(r, u);
        addContactResidual_(r, u, constrainedEquations_());
    }

    JacobianMatrix& jacobian()
    { return couplePairs_ ? *jacobian_ : assembler_->jacobian(); }

    ResidualType& residual()
    { return *residual_; }

    std::size_t numDofs() const
    { return assembler_->numDofs(); }

    const Problem& problem() const
    { return assembler_->problem(); }

    const GridGeometry& gridGeometry() const
    { return assembler_->gridGeometry(); }

    const auto& gridView() const
    { return assembler_->gridView(); }

    GridVariables& gridVariables()
    { return assembler_->gridVariables(); }

    const GridVariables& gridVariables() const
    { return assembler_->gridVariables(); }

    const SolutionVector& prevSol() const
    { return assembler_->prevSol(); }

    bool isStationaryProblem() const
    { return assembler_->isStationaryProblem(); }

    void updateGridVariables(const SolutionVector& u)
    { assembler_->updateGridVariables(u); }

    void resetTimeStep(const SolutionVector& u)
    { assembler_->resetTimeStep(u); }

    const ContactModel& contact() const
    { return *contact_; }

private:
    static constexpr int blockSize = ResidualType::block_type::dimension;

    void setJacobianPattern_()
    {
        const auto numDofs = assembler_->numDofs();
        auto pattern = getJacobianPattern<true>(assembler_->gridGeometry());
        contact_->jacobianPattern([&](std::size_t row, std::size_t col){ pattern.add(row, col); });
        jacobian_->setSize(numDofs, numDofs);
        pattern.exportIdx(*jacobian_);
    }

    // the pattern of target contains the pattern of source, both with sorted column indices
    static void copyEntries_(const JacobianMatrix& source, JacobianMatrix& target)
    {
        for (std::size_t i = 0; i < source.N(); ++i)
        {
            auto targetEntry = target[i].begin();
            const auto targetEnd = target[i].end();
            for (auto entry = source[i].begin(); entry != source[i].end(); ++entry, ++targetEntry)
            {
                for (; targetEntry.index() < entry.index(); ++targetEntry)
                    *targetEntry = 0.0;
                *targetEntry = *entry;
            }
            for (; targetEntry != targetEnd; ++targetEntry)
                *targetEntry = 0.0;
        }
    }

    std::vector<bool> constrainedEquations_() const
    {
        std::vector<bool> constrained(assembler_->numDofs()*blockSize, false);
        if constexpr (Experimental::Detail::hasGlobalConstraints<Problem>())
        {
            for (const auto& constraint : assembler_->problem().constraints())
            {
                const auto& info = constraint.constraintInfo();
                for (int eqIdx = 0; eqIdx < info.size(); ++eqIdx)
                    if (info.isConstraintEquation(eqIdx))
                        constrained[constraint.dofIndex()*blockSize + eqIdx] = true;
            }
        }
        return constrained;
    }

    void addContactResidual_(ResidualType& r, const SolutionVector& u, const std::vector<bool>& constrained) const
    {
        ResidualType contactResidual(r.size());
        contactResidual = 0.0;
        contact_->addResidual(contactResidual, u);
        for (std::size_t dofIdx = 0; dofIdx < r.size(); ++dofIdx)
            for (int i = 0; i < blockSize; ++i)
                if (!constrained[dofIdx*blockSize + i])
                    r[dofIdx][i] += contactResidual[dofIdx][i];
    }

    std::shared_ptr<Assembler> assembler_;
    std::shared_ptr<ContactModel> contact_;
    bool couplePairs_;
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualType> residual_;
};

} // end namespace Dumux::Contact

#endif
