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
 * \brief Assembler class for residuals and jacobian matrices
 *        for grid-based numerical schemes.
 */
#ifndef DUMUX_ASSEMBLER_HH
#define DUMUX_ASSEMBLER_HH

#include <cassert>

#include <dumux/discretization/method.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/timestepping/multistagetimestepper.hh>
#include <dumux/timestepping/timelevel.hh>

#include "jacobianpattern.hh"

// TODO: Remove when linear algebra traits are available
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux::Experimental {

//! Default types used for the linear system
template<class SolutionVector>
struct DefaultLinearSystemTraits;

//! Default linear system types for Dune::BlockVector
template<class PrimaryVariables>
struct DefaultLinearSystemTraits<Dune::BlockVector<PrimaryVariables>>
{
private:
    static constexpr int numEq = NumEqVectorTraits<PrimaryVariables>::numEq;
    using Scalar = typename PrimaryVariables::value_type;
    using MatrixBlockType = Dune::FieldMatrix<Scalar, numEq, numEq>;

public:
    using ResidualVector = Dune::BlockVector<PrimaryVariables>;
    using JacobianMatrix = Dune::BCRSMatrix<MatrixBlockType>;
};

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for grid-based numerical schemes
 * \tparam LA The local assembler (element-local assembly)
 * \tparam LST The linear system traits (types used for jacobians and residuals)
 */
template<class LA,
         class LST = DefaultLinearSystemTraits<typename LA::GridVariables::SolutionVector>>
class Assembler
{
    using GV = typename LA::GridVariables;
    using GG = typename GV::GridGeometry;

    using GridView = typename GG::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export the types used for the linear system
    using JacobianMatrix = typename LST::JacobianMatrix;
    using ResidualVector = typename LST::ResidualVector;
    using ResidualType = ResidualVector;

    //! export the types underlying discrete solutions
    using Scalar = typename GV::Scalar;
    using SolutionVector = typename GV::SolutionVector;

    //! export the types involved in element-local assembly
    using LocalAssembler = LA;
    using LocalOperator = typename LA::LocalOperator;

    //! the variables representing an evaluation point
    using GridVariables = GV;
    using Variables = GridVariables;

    //! export the underlying grid geometry type
    using GridGeometry = GG;

    //! export the type for parameters of a stage in time integration
    using StageParams = MultiStageParams<Scalar>;

    /*!
     * \brief The Constructor from a grid geometry.
     * \param gridGeometry A grid geometry instance
     * \param diffMethod The method used for computing partial derivatives
     * \note This assembler class is, after construction, defined for a specific equation
     *       (defined by the LocalOperator inside the LocalAssembler) and a specific grid
     *       geometry - which defines the connectivity of the degrees of freedom of the
     *       underlying discretization scheme on a particular grid. The evaluation point,
     *       consisting of a particular solution/variables/parameters may vary, and therefore,
     *       an instance of the grid variables is passed to the assembly functions.
     * \note This constructor (without time integration method) assumes hhat a stationary problem
     *       is assembled, and thus, an implicit jacobian pattern is used.
     */
    Assembler(std::shared_ptr<const GridGeometry> gridGeometry,
              DiffMethod diffMethod = DiffMethod::numeric)
    : gridGeometry_(gridGeometry)
    , isImplicit_(true)
    , diffMethod_(diffMethod)
    {}

    /*!
     * \brief The Constructor from a grid geometry.
     * \param gridGeometry A grid geometry instance
     * \param method the time integration method that will be used.
     */
    template<class TimeIntegrationMethod>
    Assembler(std::shared_ptr<const GridGeometry> gridGeometry,
              const TimeIntegrationMethod& method,
              DiffMethod diffMethod = DiffMethod::numeric)
    : gridGeometry_(gridGeometry)
    , isImplicit_(method.implicit())
    , diffMethod_(diffMethod)
    {}

    /*!
     * \brief Assembles the Jacobian matrix and the residual around the given evaluation point
     *        which is determined by the grid variables, containing all quantities required
     *        to evaluate the equations to be assembled.
     * \param gridVariables The variables corresponding to the given solution state
     * \note We assume the grid geometry on which the grid variables are defined
     *       to be the same as the one used to instantiate this class.
     * \todo TODO: The grid variables have to be non-const in order to compute
     *             the derivatives in the case of global caching. Could this be
     *             circumvented somehow?
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        assert(gridVariables.gridGeometry().numDofs() == gridGeometry().numDofs());

        resetJacobian_(partialReassembler);
        resetResidual_();

        assemble_([&](const Element& element)
        {
            auto localAssembler = this->makeLocalAssembler_(element, gridVariables);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *residual_, partialReassembler);
        });

        // TODO: Put these into discretization-specific helpers?
        enforceDirichletConstraints_(gridVariables, *jacobian_, *residual_);
        enforceInternalConstraints_(gridVariables, *jacobian_, *residual_);
        enforcePeriodicConstraints_(gridVariables, *jacobian_, *residual_);
    }

    /*!
     * \brief Assembles the Jacobian matrix of the discrete system of equations
     *        around a given state represented by the grid variables object.
     */
    void assembleJacobian(GridVariables& gridVariables)
    {
        assert(gridVariables.gridGeometry().numDofs() == gridGeometry().numDofs());

        resetJacobian_();

        assemble_([&](const Element& element)
        {
            auto localAssembler = this->makeLocalAssembler_(element, gridVariables);
            localAssembler.assembleJacobianAndResidual(*jacobian_);
        });

        // TODO: Put these into discretization-specific helpers?
        enforceDirichletConstraints_(gridVariables, *jacobian_, *residual_);
        enforceInternalConstraints_(gridVariables, *jacobian_, *residual_);
        enforcePeriodicConstraints_(gridVariables, *jacobian_, *residual_);
    }

    /*!
     * \brief Assembles the residual for a given state represented by the provided
     *        grid variables object, using the internal residual vector to store the result.
     */
    void assembleResidual(GridVariables& gridVariables)
    {
        resetResidual_();
        assembleResidual(*residual_, gridVariables);
    }

    /*!
     * \brief Assembles the residual for a given state represented by the provided
     *        grid variables object, using the provided residual vector to store the result.
     */
    void assembleResidual(ResidualVector& r, GridVariables& gridVariables) const
    {
        assemble_([&](const Element& element)
        { this->makeLocalAssembler_(element, gridVariables).assembleResidual(r); });
    }

    //! Will become obsolete once the new linear solvers are available
    Scalar residualNorm(GridVariables& gridVars) const
    {
        ResidualType residual(numDofs());
        assembleResidual(residual, gridVars);

        // for box communicate the residual with the neighboring processes
        if (GridGeometry::discMethod == DiscretizationMethod::box && gridView().comm().size() > 1)
        {
            using VertexMapper = typename GridGeometry::VertexMapper;
            VectorCommDataHandleSum<VertexMapper, SolutionVector, GridGeometry::GridView::dimension>
                sumResidualHandle(gridGeometry_->vertexMapper(), residual);
            gridView().communicate(sumResidualHandle,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);
        }

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView().comm().size() > 1)
            result2 = gridView().comm().sum(result2);

        using std::sqrt;
        return sqrt(result2);
    }


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
        // TODO: Does the bool need to be at compile time in getJacobianPattern?
        const auto occupationPattern = isImplicit_ ? getJacobianPattern<true>(gridGeometry())
                                                   : getJacobianPattern<false>(gridGeometry());

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    /*!
     * \brief Prepare for a new stage within a time integration step.
     *        This caches the given grid variables, which are then used as a
     *        representation of the previous stage. Moreover, the given grid
     *        variables are then updated to the time level of the upcoming stage.
     * \param gridVars the grid variables representing the previous stage
     * \param params the parameters with the weights to be used in the upcoming stage
     */
    void prepareStage(GridVariables& gridVars,
                      std::shared_ptr<const StageParams> params)
    {
        stageParams_ = params;
        const auto curStage = params->size() - 1;

        // we keep track of previous stages, they are needed for residual assembly
        prevStageVariables_.push_back(gridVars);

        // Now we update the time level of the given grid variables
        const auto t = params->timeAtStage(curStage);
        const auto prevT = params->timeAtStage(0);
        const auto dtFraction = params->timeStepFraction(curStage);
        gridVars.updateTime(TimeLevel<Scalar>{t, prevT, dtFraction});
    }

    /*!
     * \brief Remove traces from stages within a time integration step.
     */
    void clearStages()
    {
        prevStageVariables_.clear();
        stageParams_ = nullptr;
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

    //! The grid view
    const GridView& gridView() const
    { return gridGeometry().gridView(); }

    //! The jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! The residual vector (rhs)
    ResidualVector& residual()
    { return *residual_; }

protected:
    // make a local assembler instance
    LocalAssembler makeLocalAssembler_(const Element& element, GridVariables& gridVars) const
    {
        // instationary assembly
        if (!stageParams_)
            return LocalAssembler{element, gridVars, diffMethod_};

        // stationary assembly
        return LocalAssembler{element, prevStageVariables_, gridVars, stageParams_, diffMethod_};
    }

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
    { /*TODO: Implement / where to put this? Currently, local assemblers do this*/ }

    void enforceInternalConstraints_(const GridVariables& gridVars, JacobianMatrix& jac, SolutionVector& res)
    { /*TODO: Implement / where to put this? Currently, local assemblers do this*/ }

    void enforcePeriodicConstraints_(const GridVariables& gridVars, JacobianMatrix& jac, SolutionVector& res)
    { /*TODO: Implement / where to put this? Currently, local assemblers do this*/ }

private:
    //! the grid geometry on which it is assembled
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! states if an implicit of explicit scheme is used (affects jacobian pattern)
    bool isImplicit_;

    //! the method to compute the derivatives
    DiffMethod diffMethod_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualVector> residual_;

    //! parameters containing information on the current stage of time integration
     std::shared_ptr<const StageParams> stageParams_ = nullptr;

     //! keep track of the states of previous stages within a time integration step
     std::vector<GridVariables> prevStageVariables_;
};

} // end namespace Dumux::Experimental

#endif
