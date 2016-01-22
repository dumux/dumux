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
 * \brief The base class of models which consist of two arbitrary
 *        sub-models which are coupled
 */

#ifndef DUMUX_MULTIDOMAIN_MODEL_HH
#define DUMUX_MULTIDOMAIN_MODEL_HH

#include <dune/common/deprecated.hh>

#include "properties.hh"
#include "propertydefaults.hh"
#include "problem.hh"
#include "convergencewriter.hh"
#include "newtoncontroller.hh"

namespace Dumux
{
/*!
 * \ingroup MultidomainModel
 * \brief The base class of models which consist of two arbitrary
 *        sub-models which are coupled
 */
template<class TypeTag>
class MultiDomainModel
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubDomain1TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubDomain2TypeTag;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, Problem) SubDomainProblem1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, Problem) SubDomainProblem2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, Model) SubDomainModel1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, Model) SubDomainModel2;

    typedef typename GET_PROP_TYPE(TypeTag, SplitAndMerge) SplitAndMerge;

    enum {
        numEq1 = GET_PROP_VALUE(TypeTag, NumEq1),
        numEq2 = GET_PROP_VALUE(TypeTag, NumEq2)
    };

public:
    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The problem
     */
    void init(Problem &problem)
    {
        problem_ = &problem;

        // the two sub models have already been initialized by the
        // sub-problems!
        jacAsm_ = std::make_shared<JacobianAssembler>();
        jacAsm_->init(problem);

        uCur_ = std::make_shared<SolutionVector>(jacAsm_->gridFunctionSpace());
        uPrev_ = std::make_shared<SolutionVector>(jacAsm_->gridFunctionSpace());

        *uCur_= 0;
        *uPrev_= 0;

        SplitAndMerge::mergeSolVectors(sdModel1().curSol(),
                                       sdModel2().curSol(),
                                       *uCur_);
        SplitAndMerge::mergeSolVectors(sdModel1().prevSol(),
                                       sdModel2().prevSol(),
                                       *uPrev_);
    }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    SolutionVector &curSol()
    { return *uCur_; }

    //! \brief \copydoc curSol()
    const SolutionVector &curSol() const
    { return *uCur_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    SolutionVector &prevSol()
    { return *uPrev_; }

    //! \brief \copydoc prevSol()
    const SolutionVector &prevSol() const
    { return *uPrev_; }

    /*!
     * \brief Returns the operator assembler for the global jacobian of
     *        the problem.
     */
    JacobianAssembler &jacobianAssembler()
    { return *jacAsm_; }

    //! \brief \copydoc jacobianAssembler()
    const JacobianAssembler &jacobianAssembler() const
    { return *jacAsm_; }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem()
    { return *problem_; }

    //! \brief \copydoc problem()
    const Problem &problem() const
    { return *problem_; }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    SubDomainProblem1 &sdProblem1()
    { return problem().sdProblem1(); }

    //! \brief \copydoc sdProblem1()
    const SubDomainProblem1 &sdProblem1() const
    { return problem().sdProblem1(); }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    SubDomainProblem2 &sdProblem2()
    { return problem().sdProblem2(); }

    //! \brief \copydoc sdProblem2()
    const SubDomainProblem2 &sdProblem2() const
    { return problem().sdProblem2(); }

    /*!
     * \brief A reference to the first sub-problem's model.
     */
    SubDomainModel1 &sdModel1()
    { return sdProblem1().model(); }

    //! \brief \copydoc sdModel1()
    const SubDomainModel1 &sdModel1() const
    { return sdProblem1().model(); }

    /*!
     * \brief A reference to the second sub-problem's model.
     */
    SubDomainModel2 &sdModel2()
    { return sdProblem2().model(); }

    //! \brief \copydoc sdModel2()
    const SubDomainModel2 &sdModel2() const
    { return sdProblem2().model(); }

    //! \copydoc Dumux::ImplicitModel::update()
    bool update(NewtonMethod &solver,
                NewtonController &controller)
    {
#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().base().size(); ++i)
            Valgrind::CheckDefined(curSol().base()[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        bool converged = solver.execute(controller);
        if (!converged)
            asImp_().updateFailed();
        else
            asImp_().updateSuccessful();

#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().base().size(); ++i) {
            Valgrind::CheckDefined(curSol().base()[i]);
        }
#endif // HAVE_VALGRIND

        return converged;
    }


    //! \copydoc Dumux::ImplicitModel::checkPlausibility()
    void checkPlausibility() const
    { }

    //! \copydoc Dumux::ImplicitModel::updateBegin()
    void updateBegin()
    {
        sdModel1().updateBegin();
        sdModel2().updateBegin();

        SplitAndMerge::mergeSolVectors(sdModel1().curSol(), sdModel2().curSol(), *uCur_);
    }

    //! \copydoc Dumux::ImplicitModel::updateSuccessful()
    void updateSuccessful()
    {
        sdModel1().updateSuccessful();
        sdModel2().updateSuccessful();
    }

    /*!
     * \brief Called by the problem if a timeintegration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should perpare the model for the next time integration.
     * Note, that the advanceTimeLevel() methods of the sub-models
     * have already been called by the individual sub problems...
     */
    void advanceTimeLevel()
    {
        // merge the two sub-vectors together
        SplitAndMerge::mergeSolVectors(sdModel1().curSol(), sdModel2().curSol(), *uCur_);
        SplitAndMerge::mergeSolVectors(sdModel1().prevSol(), sdModel2().prevSol(), *uPrev_);
    }

    //! \copydoc Dumux::ImplicitModel::updateFailed()
    void updateFailed()
    {
        sdModel1().updateFailed();
        sdModel2().updateFailed();

        // merge the two sub-vectors together
        SplitAndMerge::mergeSolVectors(sdModel1().curSol(), sdModel2().curSol(), *uCur_);
    }

     /*!
      * \brief Called by the update() method if a try was
      *         unsuccessful. This is primary a hook which the
      *         actual model can overload.
      */
    void updateFailedTry()
    {
        sdModel1().updateFailedTry();
        sdModel2().updateFailedTry();

        // merge the two sub-vectors together
        SplitAndMerge::mergeSolVectors(sdModel1().curSol(), sdModel2().curSol(), *uCur_);
    }

    /*!
     * \brief Calculate the global residual.
     *
     * \param globResidual The global residual
     *
     * The global deflection of the balance equation from zero.
     */
    void evalGlobalResidual(SolutionVector &globResidual)
    {
        DUNE_THROW(Dune::NotImplemented, "");
    }

    //! \copydoc Dumux::ImplicitModel::serialize()
    template <class Restarter>
    void serialize(Restarter &res)
    {
        sdProblem1().serialize(res);
        sdProblem2().serialize(res);
    }

    //! \copydoc Dumux::ImplicitModel::deserialize()
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        sdProblem1().deserialize(res);
        sdProblem2().deserialize(res);
        wasRestarted_ = true;
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    DUNE_DEPRECATED_MSG("numDofs() is deprecated.")
    size_t numDofs() const
    {
        return sdModel1().numDofs()*numEq1 + sdModel2().numDofs()*numEq2;
    }

    //! \copydoc Dumux::ImplicitModel::resetJacobianAssembler()
    void resetJacobianAssembler()
    {
        jacAsm_.template reset<JacobianAssembler>(0);
        jacAsm_ = std::make_shared<JacobianAssembler>(asImp_(), problem());
    }


protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // the problem we want to solve. defines the constitutive
    // relations, material laws, etc.
    Problem *problem_;

    // the jacobian assembler
    std::shared_ptr<JacobianAssembler> jacAsm_;

    // cur is the current solution, prev the solution of the previous
    // time step
    std::shared_ptr<SolutionVector> uCur_;
    std::shared_ptr<SolutionVector> uPrev_;

    bool wasRestarted_;
};
} // namespace Dumux

#endif // DUMUX_MULTIDOMAIN_MODEL_HH
