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
 * \ingroup MixedDimension
 * \brief Base class for the coupled models of mixed dimension
 */
#ifndef DUMUX_MIXEDDIMENSION_MODEL_HH
#define DUMUX_MIXEDDIMENSION_MODEL_HH

#include <dune/geometry/type.hh>
#include <dune/istl/bvector.hh>

#include <dumux/mixeddimension/properties.hh>
#include <dumux/mixeddimension/assembler.hh>
#include <dumux/mixeddimension/subproblemlocaljacobian.hh>
#include <dumux/mixeddimension/newtoncontroller.hh>
#include <dumux/implicit/adaptive/gridadaptproperties.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup MixedDimension
 * \brief The base class for implicit models of mixed dimension.
 */
template<class TypeTag>
class MixedDimensionModel
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, BulkLocalJacobian) BulkLocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, LowDimLocalJacobian) LowDimLocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;
    typedef typename GET_PROP(TypeTag, SubProblemBlockIndices) SubProblemBlockIndices;

    // obtain the type tags of the sub problems
    typedef typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag) BulkProblemTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag) LowDimProblemTypeTag;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, GridView) BulkGridView;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView) LowDimGridView;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, LocalResidual) BulkLocalResidual;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, LocalResidual) LowDimLocalResidual;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, VertexMapper) BulkVertexMapper;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, VertexMapper) LowDimVertexMapper;

    typedef typename GET_PROP_TYPE(BulkProblemTypeTag, ElementMapper) BulkElementMapper;
    typedef typename GET_PROP_TYPE(LowDimProblemTypeTag, ElementMapper) LowDimElementMapper;

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    enum {
        bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox),
        lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox),
        bulkDofCodim = bulkIsBox ? bulkDim : 0,
        lowDimDofCodim = lowDimIsBox ? lowDimDim : 0
    };

    typename SubProblemBlockIndices::BulkIdx bulkIdx;
    typename SubProblemBlockIndices::LowDimIdx lowDimIdx;

    typedef typename BulkGridView::template Codim<0>::Entity BulkElement;
    typedef typename LowDimGridView::template Codim<0>::Entity LowDimElement;

    typedef typename Dune::ReferenceElements<Scalar, bulkDim> BulkReferenceElements;
    typedef typename Dune::ReferenceElements<Scalar, lowDimDim> LowDimReferenceElements;

public:
     // copying a model is not a good idea
    MixedDimensionModel(const MixedDimensionModel&) = delete;

    /*!
     * \brief The constructor.
     */
    MixedDimensionModel()
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
        const unsigned int bulkNumDofs = asImp_().bulkNumDofs();
        const unsigned int lowDimNumDofs = asImp_().lowDimNumDofs();
        uCur_[bulkIdx].resize(bulkNumDofs);
        uCur_[lowDimIdx].resize(lowDimNumDofs);
        uPrev_[bulkIdx].resize(bulkNumDofs);
        uPrev_[lowDimIdx].resize(lowDimNumDofs);

        // initialize local jocabian and jacobian assembler
        // \note these local jacobians are not the localjacobian objects of the submodels
        // \todo generalize so that there is only one unique local jac object
        bulkLocalJacobian_.init(problem_());
        lowDimLocalJacobian_.init(problem_());
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

        for (const auto& element : elements(bulkGridView_()))
        {
            bulkLocalResidual().eval(element);

            if (bulkIsBox)
            {
                for (int i = 0; i < element.subEntities(bulkDim); ++i)
                {
                    const unsigned int dofIdxGlobal = problem_().bulkProblem().model().dofMapper().subIndex(element, i, bulkDim);
                    residual[bulkIdx][dofIdxGlobal] += bulkLocalResidual().residual(i);
                }
            }
            else
            {
                const unsigned int dofIdxGlobal = problem_().bulkProblem().model().dofMapper().index(element);
                residual[bulkIdx][dofIdxGlobal] = bulkLocalResidual().residual(0);
            }
        }

        for (const auto& element : elements(lowDimGridView_()))
        {
            lowDimLocalResidual().eval(element);

            if (lowDimIsBox)
            {
                for (int i = 0; i < element.subEntities(lowDimDim); ++i)
                {
                    const unsigned int dofIdxGlobal = problem_().lowDimProblem().model().dofMapper().subIndex(element, i, lowDimDim);
                    residual[lowDimIdx][dofIdxGlobal] += lowDimLocalResidual().residual(i);
                }
            }
            else
            {
                const unsigned int dofIdxGlobal = problem_().lowDimProblem().model().dofMapper().index(element);
                residual[lowDimIdx][dofIdxGlobal] = lowDimLocalResidual().residual(0);
            }
        }

        // calculate the square norm of the residual
        return residual.two_norm();
    }

    void adaptVariableSize()
    {
        uCur_[bulkIdx].resize(asImp_().bulkNumDofs());
        uCur_[lowDimIdx].resize(asImp_().lowDimNumDofs());
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
     *        stiffness matrix for an arbitrary bulk element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    BulkLocalJacobian &bulkLocalJacobian()
    { return bulkLocalJacobian_; }
    /*!
     * \copydoc bulkLocalJacobian()
     */
    const BulkLocalJacobian &bulkLocalJacobian() const
    { return bulkLocalJacobian_; }

    /*!
     * \brief Returns the local jacobian which calculates the local
     *        stiffness matrix for an arbitrary low dimensional element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    LowDimLocalJacobian &lowDimLocalJacobian()
    { return lowDimLocalJacobian_; }
    /*!
     * \copydoc lowDimLocalJacobian()
     */
    const LowDimLocalJacobian &lowDimLocalJacobian() const
    { return lowDimLocalJacobian_; }

    /*!
     * \brief Returns the bulk local residual function.
     */
    BulkLocalResidual& bulkLocalResidual()
    { return bulkLocalJacobian().localResidual(); }
    /*!
     * \copydoc bulkLocalResidual()
     */
    const BulkLocalResidual& bulkLocalResidual() const
    { return bulkLocalJacobian().localResidual(); }

    /*!
     * \brief Returns the low dim local residual function.
     */
    LowDimLocalResidual& lowDimLocalResidual()
    { return lowDimLocalJacobian().localResidual(); }
    /*!
     * \copydoc lowDimLocalResidual()
     */
    const LowDimLocalResidual& lowDimLocalResidual() const
    { return lowDimLocalJacobian().localResidual(); }

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
        problem_().bulkProblem().model().updateBegin();
        problem_().lowDimProblem().model().updateBegin();

        if (problem_().bulkProblem().gridChanged() || problem_().lowDimProblem().gridChanged())
        {
            // resize matrix adjust prev sol
            uPrev_ = uCur_;
            //! Recompute the assembly map in the localjacobian
            bulkLocalJacobian_.init(problem_());
            lowDimLocalJacobian_.init(problem_());
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
        problem_().bulkProblem().model().newtonEndStep();
        problem_().lowDimProblem().model().newtonEndStep();

        //this is needed in case one of the subproblems has a phaseswitch during a newton step
        asImp_().copySubProblemSolutions_();
    }

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
       // call the respective methods in the sub problems
        problem_().lowDimProblem().model().updateFailed();
        problem_().bulkProblem().model().updateFailed();

        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        uCur_ = uPrev_;
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
        // call the respective methods in the sub problems
        problem_().lowDimProblem().model().advanceTimeLevel();
        problem_().bulkProblem().model().advanceTimeLevel();

        // make the current solution the previous one.
        uPrev_ = uCur_;
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs) of the bulk problem
     */
    std::size_t bulkNumDofs() const
    { return problem_().bulkProblem().model().numDofs(); }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs) of the low dim problem
     */
    std::size_t lowDimNumDofs() const
    { return problem_().lowDimProblem().model().numDofs(); }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    std::size_t numDofs() const
    { return asImp_().bulkNumDofs() + asImp_().lowDimNumDofs(); }

    /*!
     * \brief Mapper for the entities of the subproblems
     */
    const BulkVertexMapper& bulkVertexMapper() const
    {
        return problem_().bulkProblem().vertexMapper();
    }

    const BulkElementMapper& bulkElementMapper() const
    {
        return problem_().bulkProblem().elementMapper();
    }

    const LowDimVertexMapper& lowDimVertexMapper() const
    {
        return problem_().lowDimProblem().vertexMapper();
    }

    const LowDimElementMapper& lowDimElementMapper() const
    {
        return problem_().lowDimProblem().elementMapper();
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
        problem_().bulkProblem().model().curSol() = uCur_[bulkIdx];
        problem_().lowDimProblem().model().curSol() = uCur_[lowDimIdx];
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
    const BulkGridView &bulkGridView_() const
    { return problem_().bulkProblem().gridView(); }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const LowDimGridView &lowDimGridView_() const
    { return problem_().lowDimProblem().gridView(); }

    /*!
     * \brief Copies the solution vectors of the sub problems
     */
    void copySubProblemSolutions_()
    {
        uCur_[bulkIdx] = problem_().bulkProblem().model().curSol();
        uCur_[lowDimIdx] = problem_().lowDimProblem().model().curSol();
    }

    // the problem we want to solve. defines the constitutive
    // relations, matxerial laws, etc.
    Problem *problemPtr_;

    // calculates the local jacobian matrix for a given bulk/lowdim element
    BulkLocalJacobian bulkLocalJacobian_;
    LowDimLocalJacobian lowDimLocalJacobian_;

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
