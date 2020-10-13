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
 * \ingroup Assembly
 * \brief Assembler class for residuals and jacobian matrices for
 *        grid-based numerical schemes and instationary problems.
 */
#ifndef DUMUX_INSTATIONARY_ASSEMBLER_HH
#define DUMUX_INSTATIONARY_ASSEMBLER_HH

#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/localview.hh>

#include <dumux/assembly/partialreassembler.hh>
#include <dumux/timestepping/multistagetimestepper.hh>
#include <dumux/timestepping/timelevel.hh>

// TODO: include box/cc local assemblers
#include "fv/instationarylocalassembler.hh"
#include "fem/instationarylocalassembler.hh"
#include "jacobianpattern.hh"

namespace Dumux {
namespace Impl {

    template<class Assembler, DiscretizationMethod dm, DiffMethod diffMethod> struct LocalAssemblerChooser;

    // TODO: Currently the fv local residuals require type tag etc...
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::tpfa> { using type = CCLocalAssembler<... };
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::mpfa> { using type = CCLocalAssembler<... };
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::staggered> { using type = ... };
    template<class Assembler, DiffMethod diffMethod>
    struct LocalAssemblerChooser<Assembler, DiscretizationMethod::box, diffMethod>
    { using type = FVInstationaryLocalAssembler<Assembler, diffMethod>; };
    template<class Assembler, DiffMethod diffMethod>
    struct LocalAssemblerChooser<Assembler, DiscretizationMethod::fem, diffMethod>
    { using type = FEInstationaryLocalAssembler<Assembler, diffMethod>; };

    template<class Assembler, DiscretizationMethod dm, DiffMethod diffMethod>
    using LocalAssemblerType = typename LocalAssemblerChooser<Assembler, dm, diffMethod>::type;

} // end namespace detail

//! Default types used for the linear system
template<class Scalar, int numEq>
struct DefaultLinearSystemTraits
{
private:
    using PrimaryVariables = Dune::FieldVector<Scalar, numEq>;
    using BlockType = Dune::FieldMatrix<Scalar, numEq, numEq>;

public:
    using ResidualVector = Dune::BlockVector<PrimaryVariables>;
    using JacobianMatrix = Dune::BCRSMatrix<BlockType>;
};

//! forward declaration, TODO: REMOVE
template<class P> struct ProblemTraits;

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for grid-based numerical schemes
 * \tparam LO The local operator (evaluation of the terms of the equation)
 * \tparam diffMethod The differentiation method to compute derivatives
 * \tparam LST The linear system traits (types used for jacobians and residuals)
 */
template< class LO, DiffMethod diffMethod,
          class LST = DefaultLinearSystemTraits<typename LO::GridVariables::Scalar,
                                                LO::GridVariables::PrimaryVariables::size()> >
class InstationaryAssembler
{
    using ThisType = InstationaryAssembler<LO, diffMethod, LST>;

    using GG = typename LO::GridVariables::GridGeometry;
    using GridView = typename GG::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using LocalAssembler = Impl::LocalAssemblerType<ThisType, GG::discMethod, diffMethod>;
    using ElementVariables = typename LO::GridVariables::LocalView;

public:
    //! export the types used for the linear system
    using Scalar = typename LO::GridVariables::Scalar;
    using JacobianMatrix = typename LST::JacobianMatrix;
    using ResidualVector = typename LST::ResidualVector;
    using ResidualType = ResidualVector; // Required by NewtonSolver etc. (TODO: Which name is better!?)

    //! export the local operator type
    using LocalOperator = LO;

    //! the local operator states the type of variables needed for evaluation
    using GridVariables = typename LO::GridVariables;

    //! export a grid-independent alias such that e.g. the newton solver works
    //! for grid-based schemes as well as algebraic non-linear systems of equations
    using Variables = GridVariables;

    //! The variables depend on a user-defined problem (boundary conditions/parameters)
    using Problem = typename GridVariables::Problem;

    //! This assembler is for grid-based schemes, so the variables live on a grid geometry
    using GridGeometry = GG;

    //! Export type used for solution vectors (the variables depend upon the solution)
    using SolutionVector = typename GridVariables::SolutionVector;

    //! export the type for parameters of a stage in time integration
    using StageParams = MultiStageParams<Scalar>;

    /*!
     * \brief The Constructor.
     * \param gridVariables A grid variables instance
     */
    InstationaryAssembler(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    {}

    /*!
     * \brief Assembles the Jacobian matrix and the residual around the evaluation point
     *        that is determined by the grid variables that have been used for construction.
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(const GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        resetJacobian_(partialReassembler);
        resetResidual_();

        auto feGeometry = localView(gridGeometry());
        std::vector<ElementVariables> elemVars;
        elemVars.reserve(stageParams_->size());

        for (unsigned int i = 0; i < stageParams_->size()-1; ++i)
            elemVars.emplace_back(prevStageVariables_[i]);
        elemVars.emplace_back(gridVariables);

        assemble_([&](const Element& element)
        {
            feGeometry.bind(element);
            for (auto& vars : elemVars)
                vars.bind(element, feGeometry);

            LocalAssembler localAssembler(element, feGeometry, elemVars, *stageParams_);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *residual_, partialReassembler);
        });

        // TODO: Put these into discretization-specific helpers
        enforceDirichletConstraints_(gridVariables, *jacobian_, *residual_);
        enforceInternalConstraints_(*jacobian_, *residual_);
        enforcePeriodicConstraints_(*jacobian_, *residual_);
    }

    /*!
     * \brief Assembles the Jacobian matrix around the evaluation point that is
     *        determined by the grid variables that have been used for construction.
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobian(const GridVariables& gridVariables,
                          const PartialReassembler* partialReassembler = nullptr)
    {
        resetJacobian_(partialReassembler);
        // TODO: WHat to do here?
    }

    /*!
     * \brief Assembles the residual for a given state represented by the provided
     *        grid variables object, using the internal residual vector to store the result.
     */
    void assembleResidual(const GridVariables& gridVars)
    {
        resetResidual_();
        assembleResidual(*residual_, gridVars);
    }

    /*!
     * \brief Assembles the residual for a given state represented by the provided
     *        grid variables object, using the provided residual vector to store the result.
     */
    void assembleResidual(ResidualVector& r, const GridVariables& gridVariables) const
    {
        auto feGeometry = localView(gridGeometry());
        std::vector<ElementVariables> elemVars;
        elemVars.reserve(stageParams_->size());

        for (unsigned int i = 0; i < stageParams_->size()-1; ++i)
            elemVars.emplace_back(prevStageVariables_[i]);
        elemVars.emplace_back(gridVariables);

        assemble_([&](const Element& element)
        {
            feGeometry.bind(element);
            for (auto& vars : elemVars)
                vars.bind(element, feGeometry);

            LocalAssembler localAssembler(element, feGeometry, elemVars, *stageParams_);
            localAssembler.assembleResidual(r);
        });
    }

    //! TODO: Do we want to remove this interface?
    //!       Should it be the assembler's job to compute the norm?
    Scalar residualNorm(const SolutionVector& x) const
    { DUNE_THROW(Dune::NotImplemented, "ResidualNorm"); }

    /*!
     * \brief Tells the assembler which jacobian and residual to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the jacobian matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<SolutionVector> r)
    {
        jacobian_ = A;
        residual_ = r;

        // check and/or set the BCRS matrix's build mode
        if (jacobian_->buildMode() == JacobianMatrix::BuildMode::unknown)
            jacobian_->setBuildMode(JacobianMatrix::random);
        else if (jacobian_->buildMode() != JacobianMatrix::BuildMode::random)
            DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        jacobian_->setBuildMode(JacobianMatrix::random);
        residual_ = std::make_shared<SolutionVector>();

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief Resizes the jacobian and sets the jacobian' sparsity pattern.
     */
    void setJacobianPattern()
    {
        // resize the jacobian and the residual
        const auto numDofs = this->numDofs();
        jacobian_->setSize(numDofs, numDofs);

        // create occupation pattern of the jacobian
        // TODO: HOW TO DETERMINE PATTERN DEPENDING ON TIME SCHEME??
        const auto occupationPattern = getJacobianPattern</*isImplicit*/true>(gridGeometry());

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    //! Resizes the residual
    void setResidualSize()
    { residual_->resize(numDofs()); }

    //! Returns the number of degrees of freedom
    std::size_t numDofs() const
    { return gridGeometry().numDofs(); }

    //! The global finite volume geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! The gridview
    const GridView& gridView() const
    { return gridGeometry().gridView(); }

    //! The jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! The residual vector (rhs)
    ResidualVector& residual()
    { return *residual_; }

    // prepare a stage within time integration step
    void prepareStage(GridVariables& gridVars,
                      std::shared_ptr<const StageParams> params)
    {
        stageParams_ = params;
        const auto curStage = params->size() - 1;

        prevStageVariables_.push_back(gridVars);

        // make sure the variables live on the correct time level
        const auto t = params->timeAtStage(curStage);
        const auto prevT = params->timeAtStage(0);
        const auto dtFraction = params->timeStepFraction(curStage);
        TimeLevel<Scalar> timeLevel(t, prevT, dtFraction);

        gridVars.updateTime(timeLevel);
    }

    // TODO: IMPLEMENT
    void clearStages()
    {
        prevStageVariables_.clear();
    }

protected:
    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if (!residual_)
        {
            residual_ = std::make_shared<ResidualVector>();
            setResidualSize();
        }

        (*residual_) = 0.0;
    }

    // reset the Jacobian matrix to 0.0
    template <class PartialReassembler = DefaultPartialReassembler>
    void resetJacobian_(const PartialReassembler *partialReassembler = nullptr)
    {
        if (!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern();
        }

        if (partialReassembler)
            partialReassembler->resetJacobian(*this);
        else
            *jacobian_ = 0.0;
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if something throwed during assembly
     */
    template<typename AssembleElementFunc>
    void assemble_(AssembleElementFunc&& assembleElement) const
    {
        // a state that will be checked on all processes
        bool succeeded = false;

        // try assembling using the local assembly function
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                assembleElement(element);

            // if we get here, everything worked well on this process
            succeeded = true;
        }
        // throw exception if a problem ocurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
        }

        // make sure everything worked well on all processes
        if (gridView().comm().size() > 1)
            succeeded = gridView().comm().min(succeeded);

        // if not succeeded rethrow the error on all processes
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    void enforceDirichletConstraints_(const GridVariables& gridVars, JacobianMatrix& jac, SolutionVector& res)
    {
        if constexpr (GridGeometry::discMethod == DiscretizationMethod::fem)
        {
            // TODO: This is the FEM implementation. We should outsource this into
            // discretization-specific helper classes
            const auto& basis = gridGeometry().feBasis();
            const auto& curSol = gridVars.dofs();

            // container to store which DOFs are Dirichlet and values
            using PrimaryVariables = typename ProblemTraits<Problem>::PrimaryVariables;
            static constexpr auto numEq = PrimaryVariables::size();

            std::vector< std::bitset<numEq> > isDirichlet(gridGeometry().numDofs());
            SolutionVector values; values.resize(gridGeometry().numDofs());

            // container to store local coordinates of an element
            using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
            std::vector<LocalCoordinate> localDofCoords;
            std::size_t curElementIndex = 0;
            bool firstCall = true;

            auto getDirichletValues = [&] (auto localIndex,
                                           const auto& localView,
                                           const auto& intersection)
            {
                const auto& element = localView.element();
                const auto& problem = gridVars.problem();
                const auto& timeLevel = gridVars.timeLevel();

                const auto bcTypes = problem.boundaryTypes(element, intersection, timeLevel);
                if (bcTypes.hasDirichlet())
                {
                    // get the local coordinates of the dofs of the element
                    const auto eIdx = gridGeometry().elementMapper().index(element);
                    if (eIdx != curElementIndex || firstCall)
                    {
                        std::vector<double> coords;
                        const auto& interp = localView.tree().finiteElement().localInterpolation();
                        interp.interpolate([&] (const LocalCoordinate& x) { return x[0]; }, coords);

                        localDofCoords.resize(coords.size());
                        for (unsigned int i = 0; i < coords.size(); ++i)
                            localDofCoords[i][0] = coords[i];

                        interp.interpolate([&] (const LocalCoordinate& x) { return x[1]; }, coords);
                        for (unsigned int i = 0; i < coords.size(); ++i)
                            localDofCoords[i][1] = coords[i];

                        firstCall = false;
                    }

                    curElementIndex = eIdx;
                    const auto index = localView.index(localIndex);
                    for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        if (bcTypes.isDirichlet(eqIdx))
                            isDirichlet[index][eqIdx] = true;
                    values[index] = problem.dirichlet(element, localDofCoords[localIndex], timeLevel);
                }
            };

            Dune::Functions::forEachBoundaryDOF(basis, getDirichletValues);

            // modify matrix accordingly
            for (auto rIt = jac.begin(); rIt != jac.end(); ++rIt)
            {
                const auto index = rIt.index();
                if (isDirichlet[index].any())
                {
                    for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        if (isDirichlet[index][eqIdx])
                            res[index][eqIdx] = curSol[index][eqIdx] - values[index][eqIdx];

                    for (auto cIt = rIt->begin(); cIt != rIt->end(); ++cIt)
                    {
                        for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        {
                            if (isDirichlet[index][eqIdx])
                            {
                                (*cIt)[eqIdx] = 0.0;
                                if (index == cIt.index())
                                    (*cIt)[eqIdx][eqIdx] = 1.0;
                            }
                        }
                    }
                }
                // SET 0 in off-diagonals?
                // else
                // {
                //     for (auto cIt = rIt->begin(); cIt != rIt->end(); ++cIt)
                //     {
                //         for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                //         {
                //             if (isDirichlet[cIt.index()][eqIdx])
                //             {
                //                 (*cIt)[eqIdx] = 0.0;
                //                 // if (index == cIt.index())
                //                 //     (*cIt)[eqIdx][eqIdx] = 1.0;
                //             }
                //         }
                //     }
                // }
            }
        }
    }

    void enforceInternalConstraints_(JacobianMatrix& jac, SolutionVector& res)
    { /*TODO: Implement*/ }

    void enforcePeriodicConstraints_(JacobianMatrix& jac, SolutionVector& res)
    { /*TODO: Implement*/ }

private:
    //! TODO DOc this
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualVector> residual_;

    //! parameters containing information on the current stage of time integration
    std::shared_ptr<const StageParams> stageParams_;

    //! TODO: Doc this
    std::vector<GridVariables> prevStageVariables_;
};

} // namespace Dumux

#endif
