// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
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
 * \ingroup MultiDomainModel
 * \brief Base class for the coupled models of equal dimension
 */
#ifndef DUMUX_MULTIDOMAIN_MODEL_FOR_STAGGERED_HH
#define DUMUX_MULTIDOMAIN_MODEL_FOR_STAGGERED_HH

#include <dune/geometry/type.hh>
#include <dune/istl/bvector.hh>

#include <dumux/multidomain/properties.hh>
#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/subproblemlocaljacobian.hh>
#include <dumux/multidomain/newtoncontroller.hh>
#include <dumux/implicit/adaptive/gridadaptproperties.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{
/*!
 * \ingroup MultiDomainModel
 * \brief The base class for implicit models of equal dimension.
 */
template<class TypeTag>
class MultiDomainModelForStaggered
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, Model);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using JacobianAssembler = typename GET_PROP_TYPE(TypeTag, JacobianAssembler);
    using StokesLocalJacobian = typename GET_PROP_TYPE(TypeTag, StokesLocalJacobian);
    using DarcyLocalJacobian = typename GET_PROP_TYPE(TypeTag, DarcyLocalJacobian);
    using NewtonMethod = typename GET_PROP_TYPE(TypeTag, NewtonMethod);
    using NewtonController = typename GET_PROP_TYPE(TypeTag, NewtonController);
    using SubProblemBlockIndices = typename GET_PROP(TypeTag, SubProblemBlockIndices);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

    using StokesLocalResidual = typename GET_PROP_TYPE(StokesProblemTypeTag, LocalResidual);
    using DarcyLocalResidual = typename GET_PROP_TYPE(DarcyProblemTypeTag, LocalResidual);

    using StokesVertexMapper = typename GET_PROP_TYPE(StokesProblemTypeTag, VertexMapper);
    using DarcyVertexMapper = typename GET_PROP_TYPE(DarcyProblemTypeTag, VertexMapper);

    using StokesElementMapper = typename GET_PROP_TYPE(StokesProblemTypeTag, ElementMapper);
    using DarcyElementMapper = typename GET_PROP_TYPE(DarcyProblemTypeTag, ElementMapper);

    enum {
        dim = StokesGridView::dimension,
    };

    typename SubProblemBlockIndices::StokesIdx stokesIdx;
    typename SubProblemBlockIndices::DarcyIdx darcyIdx;

    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
     // copying a model is not a good idea
    MultiDomainModelForStaggered(const MultiDomainModelForStaggered&) = delete;

    /*!
     * \brief The constructor.
     */
    MultiDomainModelForStaggered()
    : problemPtr_(nullptr) {}

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     * \note at this point the sub problems are already initialized
     */
    void init(Problem &problem)
    {
        problemPtr_ = &problem;

        // resize the current and previous solution
        const unsigned int darcyNumDofs = asImp_().darcyNumDofs();
        uCur_[stokesIdx][cellCenterIdx].resize(problem_().stokesProblem().model().numCellCenterDofs());
        uCur_[stokesIdx][faceIdx].resize(problem_().stokesProblem().model().numFaceDofs());
        uCur_[darcyIdx].resize(darcyNumDofs);
        uPrev_[stokesIdx][cellCenterIdx].resize(problem_().stokesProblem().model().numCellCenterDofs());
        uPrev_[stokesIdx][faceIdx].resize(problem_().stokesProblem().model().numFaceDofs());
        uPrev_[darcyIdx].resize(darcyNumDofs);

        // initialize local jocabian and jacobian assembler
        // \note these local jacobians are not the localjacobian objects of the submodels
        // \todo generalize so that there is only one unique local jac object

        std::cout << "in init(Problem& problem) of multidomain/staggeredgrid/model " << std::endl;
        stokesLocalJacobian_.init(problem_());
        darcyLocalJacobian_.init(problem_());
        jacAsm_ = std::make_shared<JacobianAssembler>();
        jacAsm_->init(problem_());

        // set the model to the state of the initial solution
        // defined by the problem
        asImp_().copySubProblemSolutions_();

        // also set the solution of the "previous" time step to the
        // initial solution.
        uPrev_ = uCur_;
    }

    /*!
     * \brief Compute the global residual for an arbitrary solution
     *        vector.
     *
     * \param residual Stores the result
     * \param u The solution for which the residual ought to be calculated
     */
    Scalar globalResidual(SolutionVector &residual,
                          const SolutionVector &u)
    {
        SolutionVector tmp(curSol());
        curSol() = u;
        Scalar res = globalResidual(residual);
        curSol() = tmp;
        return res;
    }

    /*!
     * \brief Compute the global residual for the current solution
     *        vector.
     *
     * \param residual Stores the result
     */
    Scalar globalResidual(SolutionVector &residual)
    {
        residual = 0;

        for (const auto& element : elements(stokesGridView_()))
        {
            stokesLocalResidual().eval(element);
        }

        for (const auto& element : elements(darcyGridView_()))
        {
            darcyLocalResidual().eval(element);

            const unsigned int dofIdxGlobal = problem_().darcyProblem().model().dofMapper().index(element);
            residual[darcyIdx][dofIdxGlobal] = darcyLocalResidual().residual(0);
        }

        // calculate the square norm of the residual
        return residual.two_norm();
    }

    void adaptVariableSize()
    {
        uCur_[stokesIdx].resize(asImp_().stokesNumDofs());
        uCur_[darcyIdx].resize(asImp_().darcyNumDofs());
    }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    const SolutionVector &curSol() const
    { return uCur_; }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    SolutionVector &curSol()
    { return uCur_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    const SolutionVector &prevSol() const
    { return uPrev_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    SolutionVector &prevSol()
    { return uPrev_; }

    /*!
     * \brief Returns the operator assembler for the global jacobian of
     *        the problem.
     */
    JacobianAssembler &jacobianAssembler()
    { return *jacAsm_; }

    /*!
     * \copydoc jacobianAssembler()
     */
    const JacobianAssembler &jacobianAssembler() const
    { return *jacAsm_; }

    /*!
     * \brief Returns the local jacobian which calculates the local
     *        stiffness matrix for an arbitrary Stokes element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    StokesLocalJacobian &stokesLocalJacobian()
    { return stokesLocalJacobian_; }
    /*!
     * \copydoc stokesLocalJacobian()
     */
    const StokesLocalJacobian &stokesLocalJacobian() const
    { return stokesLocalJacobian_; }

    /*!
     * \brief Returns the local jacobian which calculates the local
     *        stiffness matrix for an arbitrary Darcy element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    DarcyLocalJacobian &darcyLocalJacobian()
    { return darcyLocalJacobian_; }
    /*!
     * \copydoc darcyLocalJacobian()
     */
    const DarcyLocalJacobian &darcyLocalJacobian() const
    { return darcyLocalJacobian_; }

    /*!
     * \brief Returns the stokes local residual function.
     */
    StokesLocalResidual& stokesLocalResidual()
    { return stokesLocalJacobian().localResidual(); }
    /*!
     * \copydoc stokesLocalResidual()
     */
    const StokesLocalResidual& stokesLocalResidual() const
    { return stokesLocalJacobian().localResidual(); }

    /*!
     * \brief Returns the low dim local residual function.
     */
    DarcyLocalResidual& darcyLocalResidual()
    { return darcyLocalJacobian().localResidual(); }
    /*!
     * \copydoc darcyLocalResidual()
     */
    const DarcyLocalResidual& darcyLocalResidual() const
    { return darcyLocalJacobian().localResidual(); }

    /*!
     * \brief Returns the maximum relative shift between two vectors of
     *        primary variables.
     *
     * \param priVars1 The first vector of primary variables
     * \param priVars2 The second vector of primary variables
     */
    template <class PrimaryVariableVector>
    Scalar relativeShiftAtDof(const PrimaryVariableVector &priVars1,
                              const PrimaryVariableVector &priVars2)
    {
        Scalar result = 0.0;
        for (int j = 0; j < priVars1.size(); ++j)
        {
            Scalar eqErr = std::abs(priVars1[j] - priVars2[j]);
            eqErr /= std::max<Scalar>(1.0, std::abs(priVars1[j] + priVars2[j])/2);

            result = std::max(result, eqErr);
        }
        return result;
    }

    /*!
     * \brief Try to progress the model to the next timestep.
     *
     * \param solver The non-linear solver
     * \param controller The controller which specifies the behaviour
     *                   of the non-linear solver
     */
    bool update(NewtonMethod &solver,
                NewtonController &controller)
    {
        asImp_().updateBegin();
        bool converged = solver.execute(controller);

        if (converged)
            asImp_().updateSuccessful();
        else
            asImp_().updateFailed();

        return converged;
    }

    /*!
     * \brief Check the plausibility of the current solution
     *
     *        This has to be done by the actual model, it knows
     *        best, what (ranges of) variables to check.
     *        This is primarily a hook
     *        which the actual model can overload.
     */
    void checkPlausibility() const
    { }

    /*!
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primarily a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    {
        problem_().stokesProblem().model().updateBegin();
        problem_().darcyProblem().model().updateBegin();

        if (problem_().stokesProblem().gridChanged() || problem_().darcyProblem().gridChanged())
        {
            // resize matrix adjust prev sol
            uPrev_ = uCur_;
            //! Recompute the assembly map in the localjacobian
            stokesLocalJacobian_.init(problem_());
            darcyLocalJacobian_.init(problem_());
            jacAsm_->init(problem_());
        }
    }


    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateSuccessful()
    {}

    void newtonEndStep()
    {
        problem_().stokesProblem().model().newtonEndStep();
        problem_().darcyProblem().model().newtonEndStep();
    }

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        uCur_ = uPrev_;

        // call the respective methods in the sub problems
        problem_().darcyProblem().model().updateFailed();
        problem_().stokesProblem().model().updateFailed();
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and
     *        the result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        // make the current solution the previous one.
        uPrev_ = uCur_;

        // call the respective methods in the sub problems
        problem_().darcyProblem().model().advanceTimeLevel();
        problem_().stokesProblem().model().advanceTimeLevel();
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs) of the stokes problem
     */
    std::size_t stokesNumDofs() const
    { return problem_().stokesProblem().model().numDofs(); }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs) of the low dim problem
     */
    std::size_t darcyNumDofs() const
    { return problem_().darcyProblem().model().numDofs(); }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    std::size_t numDofs() const
    { return asImp_().stokesNumDofs() + asImp_().darcyNumDofs(); }

    /*!
     * \brief Mapper for the entities of the subproblems
     */
    const StokesVertexMapper& stokesVertexMapper() const
    {
        return problem_().stokesProblem().vertexMapper();
    }

    const StokesElementMapper& stokesElementMapper() const
    {
        return problem_().stokesProblem().elementMapper();
    }

    const DarcyVertexMapper& darcyVertexMapper() const
    {
        return problem_().darcyProblem().vertexMapper();
    }

    const DarcyElementMapper& darcyElementMapper() const
    {
        return problem_().darcyProblem().elementMapper();
    }

    /*!
     * \brief Resets the Jacobian matrix assembler, so that the
     *        boundary types can be altered.
     */
    void resetJacobianAssembler ()
    {
        jacAsm_ = std::make_shared<JacobianAssembler>();
        jacAsm_->init(problem_());
    }

    /*!
     * \brief Copies the global solution vector to the sub problems
     */
    void copySolutionToSubProblems()
    {
        problem_().stokesProblem().model().curSol() = uCur_[stokesIdx];
        problem_().darcyProblem().model().curSol() = uCur_[darcyIdx];
    }

protected:
    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem_()
    { return *problemPtr_; }
    /*!
     * \copydoc problem_()
     */
    const Problem &problem_() const
    { return *problemPtr_; }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const StokesGridView &stokesGridView_() const
    { return problem_().stokesProblem().gridView(); }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const DarcyGridView &darcyGridView_() const
    { return problem_().darcyProblem().gridView(); }

    /*!
     * \brief Copies the solution vectors of the sub problems
     */
    void copySubProblemSolutions_()
    {
        uCur_[stokesIdx] = problem_().stokesProblem().model().curSol();
        uCur_[darcyIdx] = problem_().darcyProblem().model().curSol();
    }

    // the problem we want to solve. defines the constitutive
    // relations, matxerial laws, etc.
    Problem *problemPtr_;

    // calculates the local jacobian matrix for a given Stokes/Darcy element
    StokesLocalJacobian stokesLocalJacobian_;
    DarcyLocalJacobian darcyLocalJacobian_;

    // Linearizes the problem at the current time step using the
    // local jacobian
    std::shared_ptr<JacobianAssembler> jacAsm_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    SolutionVector uCur_;
    SolutionVector uPrev_;

private:

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

};

} // end namespace Dumux

#endif
