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
#ifndef DUMUX_MULTIDOMAIN_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_ASSEMBLER_HH

#include "multidomainproperties.hh"
#include "multidomainpropertydefaults.hh"
#include <dune/pdelab/constraints/constraintsparameters.hh>
#include <dune/pdelab/multidomain/constraints.hh>

namespace Dumux {

//! Prevents the setting of a dirichlet constraint anywhere
struct NoDirichletConstraints :
  public Dune::PDELab::DirichletConstraintsParameters
{
  template<typename I>
  bool isDirichlet(const I& intersection, const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord) const
  {
    return false;
  }
};

template<class TypeTag>
class MultiDomainAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, MDGrid) MDGrid;

    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem2TypeTag) SubTypeTag2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, Problem) SubProblem1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Problem) SubProblem2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, LocalFEMSpace) FEM1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, LocalFEMSpace) FEM2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, ScalarGridFunctionSpace) ScalarGridFunctionSpace1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, ScalarGridFunctionSpace) ScalarGridFunctionSpace2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, GridFunctionSpace) GridFunctionSpace1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, GridFunctionSpace) GridFunctionSpace2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, LocalOperator) LocalOperator1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, LocalOperator) LocalOperator2;

    typedef typename GET_PROP_TYPE(TypeTag, MDGridFunctionSpace) MDGridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, MDCondition) MDCondition;
    typedef typename GET_PROP_TYPE(TypeTag, MDSubProblem1) MDSubProblem1;
    typedef typename GET_PROP_TYPE(TypeTag, MDSubProblem2) MDSubProblem2;
    typedef typename GET_PROP_TYPE(TypeTag, MDCouplingLocalOperator) MDCouplingLocalOperator;
    typedef typename GET_PROP_TYPE(TypeTag, MDCoupling) MDCoupling;
    typedef typename GET_PROP_TYPE(TypeTag, MDConstraintsTrafo) MDConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, MDGridOperator) MDGridOperator;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;

    typedef typename GET_PROP_TYPE(SubTypeTag1, Constraints) Constraints1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Constraints) Constraints2;

    // copying the jacobian assembler is not a good idea
    MultiDomainAssembler(const MultiDomainAssembler &);

public:
    MultiDomainAssembler()
    {
        globalProblem_ = 0;
        problem1_= 0;
        problem2_= 0;
    }

    ~MultiDomainAssembler()
    { }

    /*!
     * \brief docme
     *
     * \param problem docme
     *
     */
    void init(Problem& problem)
    {
        globalProblem_ = &problem;
        problem1_ = &globalProblem_->subProblem1();
        problem2_ = &globalProblem_->subProblem2();

        fem1_ = Dune::make_shared<FEM1>();
        fem2_ = Dune::make_shared<FEM2>();

        constraints1_ = Dune::make_shared<Constraints1>();
        constraints2_ = Dune::make_shared<Constraints2>();

        scalarGridFunctionSpace1_ = Dune::make_shared<ScalarGridFunctionSpace1>(problem1_->gridView(),
        									*fem1_,
        									*constraints1_);
        scalarGridFunctionSpace2_ = Dune::make_shared<ScalarGridFunctionSpace2>(problem2_->gridView(),
        									*fem2_,
        									*constraints2_);
        // constraints store indices of ghost dofs
        constraints1_->compute_ghosts(*scalarGridFunctionSpace1_);
        constraints2_->compute_ghosts(*scalarGridFunctionSpace2_);
        Valgrind::CheckDefined(*constraints1_);
        Valgrind::CheckDefined(*constraints2_);
//        std::cerr  << __FILE__ << ":" << __LINE__ << "\n";

        gridFunctionSpace1_ = Dune::make_shared<GridFunctionSpace1>(*scalarGridFunctionSpace1_);
        gridFunctionSpace2_ = Dune::make_shared<GridFunctionSpace2>(*scalarGridFunctionSpace2_);

        mdGridFunctionSpace_ = Dune::make_shared<MDGridFunctionSpace>(globalProblem_->mdGrid(),
        											   *gridFunctionSpace1_,
        											   *gridFunctionSpace2_);

        localOperator1_ = Dune::make_shared<LocalOperator1>(problem1_->model());
        localOperator2_ = Dune::make_shared<LocalOperator2>(problem2_->model());

        condition1_ = Dune::make_shared<MDCondition>(0);
        condition2_ = Dune::make_shared<MDCondition>(1);

        mdSubProblem1_ = Dune::make_shared<MDSubProblem1>(*localOperator1_, *condition1_);
        mdSubProblem2_ = Dune::make_shared<MDSubProblem2>(*localOperator2_, *condition2_);

        couplingLocalOperator_ = Dune::make_shared<MDCouplingLocalOperator>(*globalProblem_);
        mdCoupling_ = Dune::make_shared<MDCoupling>(*mdSubProblem1_, *mdSubProblem2_, *couplingLocalOperator_);

        // TODO proper constraints stuff
        constraintsTrafo_ = Dune::make_shared<MDConstraintsTrafo>();

        NoDirichletConstraints dirichletVal;
        auto constraints = Dune::PDELab::MultiDomain::constraints<Scalar>(*mdGridFunctionSpace_,
                                                                          Dune::PDELab::MultiDomain::constrainSubProblem(*mdSubProblem1_,dirichletVal),
                                                                          Dune::PDELab::MultiDomain::constrainSubProblem(*mdSubProblem2_,dirichletVal));
        constraints.assemble(*constraintsTrafo_);

        mdGridOperator_ = Dune::make_shared<MDGridOperator>(*mdGridFunctionSpace_, *mdGridFunctionSpace_,
        									 *constraintsTrafo_, *constraintsTrafo_,
        									 *mdSubProblem1_, *mdSubProblem2_, *mdCoupling_);

        matrix_ = Dune::make_shared<JacobianMatrix>(*mdGridOperator_);
        *matrix_ = 0;

        residual_.resize(matrix_->N());
    }

    /*!
     * \brief docme
     */
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

    //! return const reference to matrix
    const JacobianMatrix &matrix() const
    { return *matrix_; }
    //! return reference to matrix
    // This is not very nice, but required for the AMG solver
    JacobianMatrix &matrix()
    { return *matrix_; }

    //! return const reference to residual
    const SolutionVector &residual() const
    { return residual_; }
    SolutionVector &residual()
    { return residual_; }

	//! return the multidomain gridfunctionspace
    MDGridFunctionSpace &gridFunctionSpace() const
    { return *mdGridFunctionSpace_; }
    MDGridFunctionSpace &mdGridFunctionSpace() const
    { return *mdGridFunctionSpace_; }
    
	//! return the multidomain constraints trafo
    MDConstraintsTrafo &constraintsTrafo() const
    { return *constraintsTrafo_; }

private:
    Problem *globalProblem_;
    SubProblem1 *problem1_;
    SubProblem2 *problem2_;

    Dune::shared_ptr<FEM1> fem1_;
    Dune::shared_ptr<FEM2> fem2_;

    Dune::shared_ptr<Constraints1> constraints1_;
    Dune::shared_ptr<Constraints2> constraints2_;

    Dune::shared_ptr<ScalarGridFunctionSpace1> scalarGridFunctionSpace1_;
    Dune::shared_ptr<ScalarGridFunctionSpace2> scalarGridFunctionSpace2_;

    Dune::shared_ptr<GridFunctionSpace1> gridFunctionSpace1_;
    Dune::shared_ptr<GridFunctionSpace2> gridFunctionSpace2_;
    Dune::shared_ptr<MDGridFunctionSpace> mdGridFunctionSpace_;

    Dune::shared_ptr<LocalOperator1> localOperator1_;
    Dune::shared_ptr<LocalOperator2> localOperator2_;

    Dune::shared_ptr<MDCondition> condition1_;
    Dune::shared_ptr<MDCondition> condition2_;

    Dune::shared_ptr<MDSubProblem1> mdSubProblem1_;
    Dune::shared_ptr<MDSubProblem2> mdSubProblem2_;

    Dune::shared_ptr<MDCouplingLocalOperator> couplingLocalOperator_;
    Dune::shared_ptr<MDCoupling> mdCoupling_;

    Dune::shared_ptr<MDConstraintsTrafo> constraintsTrafo_;
    Dune::shared_ptr<MDGridOperator> mdGridOperator_;

    Dune::shared_ptr<JacobianMatrix> matrix_;

    SolutionVector residual_;
};

} // namespace Dumux

#endif

