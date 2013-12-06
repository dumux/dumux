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

/*
 * \file
* \brief docme
*/

#ifndef DUMUX_MULTIDOMAIN_MODEL_HH
#define DUMUX_MULTIDOMAIN_MODEL_HH

#include "multidomainproperties.hh"
#include "multidomainpropertydefaults.hh"
#include "multidomainproblem.hh"
#include "multidomainnewtoncontroller.hh"
//#include "coupledjacobianassembler.hh"

/*
* \brief docme
*/

namespace Dumux
{

/*!
 * \defgroup ModelCoupling Coupled implicit models
 */


/*!
 * \ingroup ModelCoupling
 *
 * \brief The base class of models which consist of two arbitrary
 *        sub-models which are coupled
 */

/*
* \brief docme
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

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubTypeTag2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, Problem) SubProblem1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Problem) SubProblem2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, Model) SubModel1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Model) SubModel2;

    enum {
        numEq1 = GET_PROP_VALUE(TypeTag, NumEq1),
        numEq2 = GET_PROP_VALUE(TypeTag, NumEq2)
    };

public:
    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem docme
     *
     */
    void init(Problem &problem)
    {
        problem_ = &problem;

        // the two sub models have already been initialized by the
        // sub-problems!
        jacAsm_ = new JacobianAssembler();
        jacAsm_->init(problem);

//        uCur_.resize(asImp_().numDofs());
//        uPrev_.resize(asImp_().numDofs());

        uCur_.resize(jacAsm_->residual().size());
        uPrev_.resize(jacAsm_->residual().size());

        uCur_= 0;
        uPrev_= 0;

        typedef Dumux::SplitAndMerge<TypeTag> Common;
        Common::mergeSolVectors(subModel1().curSol(),
                                 subModel2().curSol(),
                                 uCur_);
        Common::mergeSolVectors(subModel1().prevSol(),
                                 subModel2().prevSol(),
                                 uPrev_);
    }

    /*
    * \brief docme
    * \param u docme
    * \param tmp docme
    */

    Scalar globalResidual(const SolutionVector &u, SolutionVector &tmp)
    {
        DUNE_THROW(Dune::NotImplemented, "");
#if 0
          SolutionVector tmpU(asImp_(), 0.0);
        tmpU = uCur_;
        uCur_ = u;
        localJacobian_.evalGlobalResidual(tmp);

        Scalar result = tmp.two_norm();
        /*
        Scalar result = 0;
        for (int i = 0; i < (*tmp).size(); ++i) {
            for (int j = 0; j < numEq; ++j)
                result += std::abs((*tmp)[i][j]);
        }
        */
        uCur_ = tmpU;
        return result;
#endif
    };

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
    const JacobianAssembler &jacobianAssembler() const
    { return *jacAsm_; }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem()
    { return *problem_; }
    /*!
     * \copydoc problem()
     */
    const Problem &problem() const
    { return *problem_; }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    SubProblem1 &subProblem1()
    { return problem().subProblem1(); }
    /*!
     * \copydoc subProblem1()
     */
    const SubProblem1 &subProblem1() const
    { return problem().subProblem1(); }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    SubProblem2 &subProblem2()
    { return problem().subProblem2(); }
    /*!
     * \copydoc subProblem2()
     */
    const SubProblem2 &subProblem2() const
    { return problem().subProblem2(); }

    /*!
     * \brief A reference to the first sub-problem's model.
     */
    SubModel1 &subModel1()
    { return subProblem1().model(); }
    /*!
     * \copydoc subModel1()
     */
    const SubModel1 &subModel1() const
    { return subProblem1().model(); }

    /*!
     * \brief A reference to the second sub-problem's model.
     */
    SubModel2 &subModel2()
    { return subProblem2().model(); }
    /*!
     * \copydoc subModel2()
     */
    const SubModel2 &subModel2() const
    { return subProblem2().model(); }

    /*!
     * \brief Try to progress the model to the next timestep.
     *
     * \param solver docme
     * \param controller docme
     *
     */
    bool update(NewtonMethod &solver,
                NewtonController &controller)
    {
#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().size(); ++i)
            Valgrind::CheckDefined(curSol()[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        bool converged = solver.execute(controller);
        if (!converged)
            asImp_().updateFailed();
        else
            asImp_().updateSuccessful();

#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().size(); ++i) {
            Valgrind::CheckDefined(curSol()[i]);
        }
#endif // HAVE_VALGRIND

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
     *        apply the newton method. This is primary a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    {
        subModel1().updateBegin();
        subModel2().updateBegin();

        typedef Dumux::SplitAndMerge<TypeTag> Common;
        Common::mergeSolVectors(subModel1().curSol(), subModel2().curSol(), uCur_);
    }


    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateSuccessful()
    {
        subModel1().updateSuccessful();
        subModel2().updateSuccessful();
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
        typedef Dumux::SplitAndMerge<TypeTag> Common;
        // merge the two sub-vectors together
        Common::mergeSolVectors(subModel1().curSol(), subModel2().curSol(), uCur_);
        Common::mergeSolVectors(subModel1().prevSol(), subModel2().prevSol(), uPrev_);
    };

    /*!
     * \brief Called by the update() method if a try was ultimately
     *        unsuccessful. This is primary a hook which the
     *        actual model can overload.
     */
    void updateFailed()
    {
        subModel1().updateFailed();
        subModel2().updateFailed();

        typedef Dumux::SplitAndMerge<TypeTag> Common;
        // merge the two sub-vectors together
        Common::mergeSolVectors(subModel1().curSol(), subModel2().curSol(), uCur_);
    };

    /*!
     * \brief Called by the update() method if a try was
     *         unsuccessful. This is primary a hook which the
     *         actual model can overload.
     */
    void updateFailedTry()
    {
        subModel1().updateFailedTry();
        subModel2().updateFailedTry();

        typedef Dumux::SplitAndMerge<TypeTag> Common;
        // merge the two sub-vectors together
        Common::mergeSolVectors(subModel1().curSol(), subModel2().curSol(), uCur_);
    };

    /*!
     * \brief Calculate the global residual.
     *
     * \param globResidual docme
     *
     * The global deflection of the mass balance from zero.
     */
    void evalGlobalResidual(SolutionVector &globResidual)
    {
        DUNE_THROW(Dune::NotImplemented, "");
    }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \param res docme
     *
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        subProblem1().serialize(res);
        subProblem2().serialize(res);
    }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \param res docme
     *
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        subProblem1().deserialize(res);
        subProblem2().deserialize(res);
        wasRestarted_ = true;
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDofs() const
    {
		return subModel1().numDofs()*numEq1 + subModel2().numDofs()*numEq2;
    }

    void resetJacobianAssembler()
    {
        delete jacAsm_;
        jacAsm_ = new JacobianAssembler(asImp_(), problem());
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
    JacobianAssembler *jacAsm_;

    // cur is the current solution, prev the solution of the previous
    // time step
    SolutionVector uCur_;
    SolutionVector uPrev_;

    bool wasRestarted_;
};
}

#endif
