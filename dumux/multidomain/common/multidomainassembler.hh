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
 * \brief An assembler for the global Jacobian matrix for multidomain models.
 */

#ifndef DUMUX_MULTIDOMAIN_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_ASSEMBLER_HH

#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/multidomain/constraints.hh>

#include "multidomainproperties.hh"
#include "multidomainpropertydefaults.hh"

namespace Dumux {

/*!
 * \ingroup MultidomainModel
 * \brief An assembler for the global Jacobian matrix for multidomain models.
 */
template<class TypeTag>
class MultiDomainAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubDomain1TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubDomain2TypeTag;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, Problem) SubDomainProblem1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, Problem) SubDomainProblem2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, LocalFEMSpace) FEM1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, LocalFEMSpace) FEM2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, ScalarGridFunctionSpace) ScalarGridFunctionSpace1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, ScalarGridFunctionSpace) ScalarGridFunctionSpace2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, GridFunctionSpace) GridFunctionSpace1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, GridFunctionSpace) GridFunctionSpace2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, LocalOperator) LocalOperator1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, LocalOperator) LocalOperator2;

    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGridFunctionSpace) MultiDomainGridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainCondition) MultiDomainCondition;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainSubProblem1) MultiDomainSubProblem1;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainSubProblem2) MultiDomainSubProblem2;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainCouplingLocalOperator) MultiDomainCouplingLocalOperator;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainCoupling) MultiDomainCoupling;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainConstraintsTrafo) MultiDomainConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGridOperator) MultiDomainGridOperator;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;

    // copying the jacobian assembler is not a good idea
    MultiDomainAssembler(const MultiDomainAssembler &);

public:
    //! \brief The constructor
    MultiDomainAssembler()
    {
        globalProblem_ = 0;
        sdProblem1_= 0;
        sdProblem2_= 0;
    }

    //! \brief The destructor
    ~MultiDomainAssembler()
    { }

    //! \copydoc ImplicitAssembler::init()
    void init(Problem& problem)
    {
        globalProblem_ = &problem;
        sdProblem1_ = &globalProblem_->sdProblem1();
        sdProblem2_ = &globalProblem_->sdProblem2();

        fem1_ = std::make_shared<FEM1>(globalProblem_->sdGridView1());
        fem2_ = std::make_shared<FEM2>(globalProblem_->sdGridView2());

        scalarGridFunctionSpace1_ = std::make_shared<ScalarGridFunctionSpace1>(globalProblem_->sdGridView1(),
                                                                                *fem1_);
        scalarGridFunctionSpace2_ = std::make_shared<ScalarGridFunctionSpace2>(globalProblem_->sdGridView2(),
                                                                                *fem2_);

        gridFunctionSpace1_ = std::make_shared<GridFunctionSpace1>(*scalarGridFunctionSpace1_);
        gridFunctionSpace2_ = std::make_shared<GridFunctionSpace2>(*scalarGridFunctionSpace2_);

        mdGridFunctionSpace_ = std::make_shared<MultiDomainGridFunctionSpace>(globalProblem_->mdGrid(),
                                                                               *gridFunctionSpace1_,
                                                                               *gridFunctionSpace2_);

        localOperator1_ = std::make_shared<LocalOperator1>(sdProblem1_->model());
        localOperator2_ = std::make_shared<LocalOperator2>(sdProblem2_->model());

        condition1_ = std::make_shared<MultiDomainCondition>(0);
        condition2_ = std::make_shared<MultiDomainCondition>(1);

        mdSubProblem1_ = std::make_shared<MultiDomainSubProblem1>(*localOperator1_, *condition1_);
        mdSubProblem2_ = std::make_shared<MultiDomainSubProblem2>(*localOperator2_, *condition2_);

        couplingLocalOperator_ = std::make_shared<MultiDomainCouplingLocalOperator>(*globalProblem_);
        mdCoupling_ = std::make_shared<MultiDomainCoupling>(*mdSubProblem1_, *mdSubProblem2_, *couplingLocalOperator_);

        constraintsTrafo_ = std::make_shared<MultiDomainConstraintsTrafo>();

        mdGridOperator_ = std::make_shared<MultiDomainGridOperator>(*mdGridFunctionSpace_, *mdGridFunctionSpace_,
                                                                     *constraintsTrafo_, *constraintsTrafo_,
                                                                     *mdSubProblem1_, *mdSubProblem2_, *mdCoupling_);

        matrix_ = std::make_shared<JacobianMatrix>(*mdGridOperator_);
        *matrix_ = 0;

        residual_ = std::make_shared<SolutionVector>(*mdGridFunctionSpace_);
    }

    //! \copydoc ImplicitAssembler::assemble()
    void assemble()
    {
        // assemble the matrix
        *matrix_ = 0;
        mdGridOperator_->jacobian(globalProblem_->model().curSol(), *matrix_);

        // calculate the global residual
        *residual_ = 0;
        mdGridOperator_->residual(globalProblem_->model().curSol(), *residual_);
    }

    //! \copydoc ImplicitAssembler::reassembleAll()
    void reassembleAll()
    { }

    //! \copydoc ImplicitAssembler::matrix()
    const JacobianMatrix &matrix() const
    { return *matrix_; }
    JacobianMatrix &matrix()
    { return *matrix_; }

    //! \copydoc ImplicitAssembler::residual()
    const SolutionVector &residual() const
    { return *residual_; }
    SolutionVector &residual()
    { return *residual_; }

    /*!
     * \brief Return constant reference to the multidomain gridfunctionspace
     */
    MultiDomainGridFunctionSpace &gridFunctionSpace() const
    { return *mdGridFunctionSpace_; }

    /*!
     * \brief Return constant reference to the multidomain gridfunctionspace
     */
    MultiDomainGridFunctionSpace &mdGridFunctionSpace() const
    { return *mdGridFunctionSpace_; }

    /*!
     * \brief Return the multidomain constraints transformation
     */
    MultiDomainConstraintsTrafo &constraintsTrafo() const
    { return *constraintsTrafo_; }

private:
    Problem *globalProblem_;
    SubDomainProblem1 *sdProblem1_;
    SubDomainProblem2 *sdProblem2_;

    std::shared_ptr<FEM1> fem1_;
    std::shared_ptr<FEM2> fem2_;

    std::shared_ptr<ScalarGridFunctionSpace1> scalarGridFunctionSpace1_;
    std::shared_ptr<ScalarGridFunctionSpace2> scalarGridFunctionSpace2_;

    std::shared_ptr<GridFunctionSpace1> gridFunctionSpace1_;
    std::shared_ptr<GridFunctionSpace2> gridFunctionSpace2_;
    std::shared_ptr<MultiDomainGridFunctionSpace> mdGridFunctionSpace_;

    std::shared_ptr<LocalOperator1> localOperator1_;
    std::shared_ptr<LocalOperator2> localOperator2_;

    std::shared_ptr<MultiDomainCondition> condition1_;
    std::shared_ptr<MultiDomainCondition> condition2_;

    std::shared_ptr<MultiDomainSubProblem1> mdSubProblem1_;
    std::shared_ptr<MultiDomainSubProblem2> mdSubProblem2_;

    std::shared_ptr<MultiDomainCouplingLocalOperator> couplingLocalOperator_;
    std::shared_ptr<MultiDomainCoupling> mdCoupling_;

    std::shared_ptr<MultiDomainConstraintsTrafo> constraintsTrafo_;
    std::shared_ptr<MultiDomainGridOperator> mdGridOperator_;

    std::shared_ptr<JacobianMatrix> matrix_;

    std::shared_ptr<SolutionVector> residual_;
};

} // namespace Dumux

#endif // DUMUX_MULTIDOMAIN_ASSEMBLER_HH
