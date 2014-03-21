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

#include "multidomainproperties.hh"
#include "multidomainpropertydefaults.hh"
#include <dune/pdelab/constraints/constraintsparameters.hh>
#include <dune/pdelab/multidomain/constraints.hh>


namespace Dumux {

/*!
 * \brief Prevents the setting of a Dirichlet constraint anywhere
 */
struct NoDirichletConstraints :
  public Dune::PDELab::DirichletConstraintsParameters
{
  template<typename I>
  bool isDirichlet(const I& intersection, const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord) const
  {
    return false;
  }
};

/*!
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

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, Constraints) Constraints1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, Constraints) Constraints2;

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

    //! \copydoc Dumux::ImplicitModel::init()
    void init(Problem& problem)
    {
        globalProblem_ = &problem;
        sdProblem1_ = &globalProblem_->sdProblem1();
        sdProblem2_ = &globalProblem_->sdProblem2();

        fem1_ = Dune::make_shared<FEM1>();
        fem2_ = Dune::make_shared<FEM2>();

        constraints1_ = Dune::make_shared<Constraints1>();
        constraints2_ = Dune::make_shared<Constraints2>();

        scalarGridFunctionSpace1_ = Dune::make_shared<ScalarGridFunctionSpace1>(globalProblem_->sdGridView1(),
                                                                                *fem1_, *constraints1_);
        scalarGridFunctionSpace2_ = Dune::make_shared<ScalarGridFunctionSpace2>(globalProblem_->sdGridView2(),
                                                                                *fem2_, *constraints2_);
        // constraints store indices of ghost dofs
        constraints1_->compute_ghosts(*scalarGridFunctionSpace1_);
        constraints2_->compute_ghosts(*scalarGridFunctionSpace2_);
        Valgrind::CheckDefined(*constraints1_);
        Valgrind::CheckDefined(*constraints2_);
//        std::cerr  << __FILE__ << ":" << __LINE__ << "\n";

        gridFunctionSpace1_ = Dune::make_shared<GridFunctionSpace1>(*scalarGridFunctionSpace1_);
        gridFunctionSpace2_ = Dune::make_shared<GridFunctionSpace2>(*scalarGridFunctionSpace2_);

        mdGridFunctionSpace_ = Dune::make_shared<MultiDomainGridFunctionSpace>(globalProblem_->mdGrid(),
                                                                               *gridFunctionSpace1_,
                                                                               *gridFunctionSpace2_);

        localOperator1_ = Dune::make_shared<LocalOperator1>(sdProblem1_->model());
        localOperator2_ = Dune::make_shared<LocalOperator2>(sdProblem2_->model());

        condition1_ = Dune::make_shared<MultiDomainCondition>(0);
        condition2_ = Dune::make_shared<MultiDomainCondition>(1);

        mdSubProblem1_ = Dune::make_shared<MultiDomainSubProblem1>(*localOperator1_, *condition1_);
        mdSubProblem2_ = Dune::make_shared<MultiDomainSubProblem2>(*localOperator2_, *condition2_);

        couplingLocalOperator_ = Dune::make_shared<MultiDomainCouplingLocalOperator>(*globalProblem_);
        mdCoupling_ = Dune::make_shared<MultiDomainCoupling>(*mdSubProblem1_, *mdSubProblem2_, *couplingLocalOperator_);

        // TODO proper constraints stuff
        constraintsTrafo_ = Dune::make_shared<MultiDomainConstraintsTrafo>();

        NoDirichletConstraints dirichletVal;
        auto constraints = Dune::PDELab::MultiDomain::constraints<Scalar>(*mdGridFunctionSpace_,
                                                                          Dune::PDELab::MultiDomain::constrainSubProblem(*mdSubProblem1_,
                                                                                                                         dirichletVal),
                                                                          Dune::PDELab::MultiDomain::constrainSubProblem(*mdSubProblem2_,
                                                                                                                         dirichletVal));
        constraints.assemble(*constraintsTrafo_);

        mdGridOperator_ = Dune::make_shared<MultiDomainGridOperator>(*mdGridFunctionSpace_, *mdGridFunctionSpace_,
                                                                     *constraintsTrafo_, *constraintsTrafo_,
                                                                     *mdSubProblem1_, *mdSubProblem2_, *mdCoupling_);

        matrix_ = Dune::make_shared<JacobianMatrix>(*mdGridOperator_);
        *matrix_ = 0;

        residual_.resize(matrix_->N());
    }

    //! \copydoc Dumux::ImplicitModel::assemble()
    void assemble()
    {
//        std::cerr  << __FILE__ << ":" << __LINE__ << "\n";

    	// assemble the matrix
    	*matrix_ = 0;

        residual_ = 0;
    	mdGridOperator_->jacobian(globalProblem_->model().curSol(), *matrix_);
//    	printmatrix(std::cout, matrix_->base(), "global stiffness matrix", "row", 11, 3);

    	// calculate the global residual
    	residual_ = 0;
    	mdGridOperator_->residual(globalProblem_->model().curSol(), residual_);
//    	printvector(std::cout, residual_, "residual", "row", 200, 1, 3);
    }

    //! \copydoc Dumux::ImplicitModel::reassembleAll()
    void reassembleAll()
    { }

    //! \copydoc Dumux::ImplicitModel::matrix()
    const JacobianMatrix &matrix() const
    { return *matrix_; }

    /*!
     * \brief Return constant reference to global Jacobian matrix.
     *        This is not very nice, but required for the AMG solver
     *
     * \return A const reference to matrix
     */
    JacobianMatrix &matrix()
    { return *matrix_; }

    //! \copydoc Dumux::ImplicitModel::residual()
    const SolutionVector &residual() const
    { return residual_; }
    SolutionVector &residual()
    { return residual_; }

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

    Dune::shared_ptr<FEM1> fem1_;
    Dune::shared_ptr<FEM2> fem2_;

    Dune::shared_ptr<Constraints1> constraints1_;
    Dune::shared_ptr<Constraints2> constraints2_;

    Dune::shared_ptr<ScalarGridFunctionSpace1> scalarGridFunctionSpace1_;
    Dune::shared_ptr<ScalarGridFunctionSpace2> scalarGridFunctionSpace2_;

    Dune::shared_ptr<GridFunctionSpace1> gridFunctionSpace1_;
    Dune::shared_ptr<GridFunctionSpace2> gridFunctionSpace2_;
    Dune::shared_ptr<MultiDomainGridFunctionSpace> mdGridFunctionSpace_;

    Dune::shared_ptr<LocalOperator1> localOperator1_;
    Dune::shared_ptr<LocalOperator2> localOperator2_;

    Dune::shared_ptr<MultiDomainCondition> condition1_;
    Dune::shared_ptr<MultiDomainCondition> condition2_;

    Dune::shared_ptr<MultiDomainSubProblem1> mdSubProblem1_;
    Dune::shared_ptr<MultiDomainSubProblem2> mdSubProblem2_;

    Dune::shared_ptr<MultiDomainCouplingLocalOperator> couplingLocalOperator_;
    Dune::shared_ptr<MultiDomainCoupling> mdCoupling_;

    Dune::shared_ptr<MultiDomainConstraintsTrafo> constraintsTrafo_;
    Dune::shared_ptr<MultiDomainGridOperator> mdGridOperator_;

    Dune::shared_ptr<JacobianMatrix> matrix_;

    SolutionVector residual_;
};

} // namespace Dumux

#endif // DUMUX_MULTIDOMAIN_ASSEMBLER_HH
