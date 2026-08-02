// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Assembles the residual and Jacobian rows belonging to one subdomain
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_ASSEMBLER_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_ASSEMBLER_HH

#include <memory>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

#include <dumux/nonlinear/preconditioning/matrixview.hh>
#include <dumux/nonlinear/preconditioning/partition.hh>
#include <dumux/nonlinear/preconditioning/residualview.hh>

namespace Dumux::NonlinearPreconditioning {

/*!
 * \ingroup Nonlinear
 * \brief Assembles into a compressed subdomain view by restricting the element loop
 *
 * The element loop range alone determines what is frozen: an element outside the range is never given a
 * turn, so its primary variables are never perturbed and no Jacobian column is ever opened for it. Its
 * value is still read whenever an in-range element's residual needs it, so the true two-point flux law
 * is used unmodified at the artificial cut. No boundary condition is imposed and the user's problem is
 * not wrapped.
 *
 * Two loop ranges are meaningful. Restricting to \f$ N_i \f$ yields the square block
 * \f$ A_i = R_i J R_i^T \f$ that the local nonlinear solve needs. Extending by the fringe to
 * \f$ N_i \cup \Gamma_i \f$ additionally opens the columns outside \f$ N_i \f$, giving the rectangular
 * \f$ R_i J \f$ with global columns that the exact RASPEN Jacobian needs.
 *
 * \note The user's problem, grid geometry and grid variables are shared, not copied. The assembler holds
 *       them by pointer and forwards to them so that the unmodified local assemblers see exactly the
 *       types and objects the monolithic assembler would.
 */
template<class TypeTag, DiffMethod diffMethod, bool isImplicit = true>
class SubdomainAssembler
{
    using ThisType = SubdomainAssembler<TypeTag, diffMethod, isImplicit>;
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using LocalAssembler = typename Dumux::Detail::LocalAssemblerChooser_t<TypeTag, ThisType, diffMethod, isImplicit>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

    static_assert(GG::discMethod == DiscretizationMethods::cctpfa,
                  "Nonlinear domain decomposition currently supports cell-centred TPFA only");
    static_assert(diffMethod == DiffMethod::numeric,
                  "Nonlinear domain decomposition currently supports numeric differentiation only");

public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GG;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ResidualType = SubdomainResidualView<PrimaryVariables>;
    using JacobianMatrix = SubdomainMatrixView<Dune::FieldMatrix<Scalar, numEq, numEq>>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Index = Partition::Index;

    /*!
     * \brief Assemble the rows of one subdomain of a stationary problem
     * \param includeFringe whether to extend the loop to \f$ \Gamma_i \f$, opening the global columns
     */
    using TimeLoop = TimeLoopBase<Scalar>;

    SubdomainAssembler(std::shared_ptr<const Problem> problem,
                       std::shared_ptr<const GridGeometry> gridGeometry,
                       std::shared_ptr<GridVariables> gridVariables,
                       std::shared_ptr<const Partition> partition,
                       const std::vector<std::vector<Index>>& influence,
                       SubdomainIndex subdomain,
                       bool includeFringe)
    : problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , partition_(partition)
    , subdomain_(subdomain)
    , elements_(includeFringe ? partition->dofsWithFringe(subdomain) : partition->dofs(subdomain))
    , jacobian_(*partition, influence, subdomain, includeFringe)
    , residual_(*partition, subdomain)
    {}

    //! transient problems additionally need the time level and the previous solution
    void setTimeLoop(std::shared_ptr<const TimeLoop> timeLoop) { timeLoop_ = timeLoop; }
    void setPreviousSolution(const SolutionVector& prevSol) { prevSol_ = &prevSol; }

    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        jacobian_.setToZero();
        residual_.setToZero();

        for (const auto eIdx : elements_)
        {
            const auto element = gridGeometry_->element(eIdx);
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(jacobian_, residual_, *gridVariables_,
                                                       static_cast<const DefaultPartialReassembler*>(nullptr));
        }
    }

    void assembleResidual(const SolutionVector& curSol)
    {
        residual_.setToZero();

        for (const auto eIdx : partition_->dofs(subdomain_))
        {
            const auto element = gridGeometry_->element(eIdx);
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleResidual(residual_);
        }
    }

    JacobianMatrix& jacobian() { return jacobian_; }
    const JacobianMatrix& jacobian() const { return jacobian_; }
    ResidualType& residual() { return residual_; }
    const ResidualType& residual() const { return residual_; }

    const Partition& partition() const { return *partition_; }
    SubdomainIndex subdomain() const { return subdomain_; }
    std::size_t numDofs() const { return partition_->dofs(subdomain_).size(); }

    //! the interface the unmodified local assemblers expect of an assembler
    const Problem& problem() const { return *problem_; }
    const GridGeometry& gridGeometry() const { return *gridGeometry_; }
    GridVariables& gridVariables() { return *gridVariables_; }
    const GridVariables& gridVariables() const { return *gridVariables_; }
    LocalResidual localResidual() const { return LocalResidual(problem_.get(), timeLoop_.get()); }
    bool isStationaryProblem() const { return !static_cast<bool>(timeLoop_); }
    const SolutionVector& prevSol() const { return *prevSol_; }

private:
    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<const GridGeometry> gridGeometry_;
    std::shared_ptr<GridVariables> gridVariables_;
    std::shared_ptr<const Partition> partition_;
    SubdomainIndex subdomain_;
    std::vector<Index> elements_;

    JacobianMatrix jacobian_;
    ResidualType residual_;

    std::shared_ptr<const TimeLoop> timeLoop_;
    const SolutionVector* prevSol_ = nullptr;
};

} // end namespace Dumux::NonlinearPreconditioning

#endif
