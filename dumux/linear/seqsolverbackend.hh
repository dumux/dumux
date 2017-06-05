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
 * \brief Dumux solver backend
 */
#ifndef DUMUX_SOLVER_BACKEND_HH
#define DUMUX_SOLVER_BACKEND_HH

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/basicproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dumux/linear/amgproperties.hh>

namespace Dumux
{

/*!
 * \ingroup Linear
 * \brief A general solver backend allowing arbitrary preconditioners and solvers.
 *
 * This class is used as a base class for specific solver-preconditioner
 * combinations. Several parameters from the group LinearSolver are read to
 * customize the solver and preconditioner:
 *
 * Verbosity: determines how verbose the linear solver should print output.
 *
 * MaxIterations: the maximum number of iterations for the linear solver.
 *
 * ResidualReduction: the threshold for declaration of convergence.
 *
 * PreconditionerRelaxation: relaxation parameter for the preconditioner.
 *
 * PreconditionerIterations: usually specifies the number of times the
 * preconditioner is applied. In case of ILU(n), it specifies the order of the
 * applied ILU.
 */
template <class TypeTag>
class IterativePrecondSolverBackend
{
public:

  template<class Preconditioner, class Solver, class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
    const int maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
    const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);

    Vector bTmp(b);

    const double relaxation = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, PreconditionerRelaxation);
    const int precondIter = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, PreconditionerIterations);

    Preconditioner precond(A, precondIter, relaxation);

    typedef Dune::MatrixAdapter<Matrix, Vector, Vector> MatrixAdapter;
    MatrixAdapter operatorA(A);

    Solver solver(operatorA, precond, residReduction, maxIter, verbosity);

    solver.apply(x, bTmp, result_);

    return result_.converged;
  }

  // solve with RestartedGMRes (needs restartGMRes as additional argument)
  template<class Preconditioner, class Solver, class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b, const int restartGMRes)
  {
    int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
    const int maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
    const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);

    Vector bTmp(b);

    const double relaxation = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, PreconditionerRelaxation);
    const int precondIter = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, PreconditionerIterations);

    Preconditioner precond(A, precondIter, relaxation);

    typedef Dune::MatrixAdapter<Matrix, Vector, Vector> MatrixAdapter;
    MatrixAdapter operatorA(A);

    Solver solver(operatorA, precond, residReduction, restartGMRes, maxIter, verbosity);

    solver.apply(x, bTmp, result_);

    return result_.converged;
  }

  const Dune::InverseOperatorResult& result() const
  {
    return result_;
  }

private:
  Dune::InverseOperatorResult result_;
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(n)-preconditioned BiCSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. It can be applied to
 * nonsymmetric matrices.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: ILU(n) incomplete LU factorization. The order n can be
 * provided by the parameter LinearSolver.PreconditionerIterations and controls
 * the fill-in. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class ILUnBiCGSTABBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  ILUnBiCGSTABBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqILUn<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::BiCGSTABSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential SOR-preconditioned BiCSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. It can be applied to
 * nonsymmetric matrices.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: SOR successive overrelaxation method. The relaxation is
 * controlled by the parameter LinearSolver.PreconditionerRelaxation. In each
 * preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class SORBiCGSTABBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  SORBiCGSTABBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqSOR<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::BiCGSTABSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential SSOR-preconditioned BiCGSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: SSOR symmetric successive overrelaxation method. The
 * relaxation is controlled by the parameter LinearSolver.PreconditionerRelaxation.
 * In each preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class SSORBiCGSTABBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  SSORBiCGSTABBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqSSOR<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::BiCGSTABSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential GS-preconditioned BiCGSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: GS Gauss-Seidel method. It can be damped by the relaxation
 * parameter LinearSolver.PreconditionerRelaxation. In each preconditioning step,
 * it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class GSBiCGSTABBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  GSBiCGSTABBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqGS<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::BiCGSTABSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential Jacobi-preconditioned BiCSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: Jacobi method. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation. In each preconditioning step, it is
 * applied as often as given by the parameter LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class JacBiCGSTABBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  JacBiCGSTABBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqJac<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::BiCGSTABSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(n)-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: ILU(n) incomplete LU factorization. The order n can be
 * provided by the parameter LinearSolver.PreconditionerIterations and controls
 * the fill-in. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class ILUnCGBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  ILUnCGBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqILUn<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::CGSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential SOR-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: SOR successive overrelaxation method. The relaxation is
 * controlled by the parameter LinearSolver.PreconditionerRelaxation. In each
 * preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class SORCGBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  SORCGBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqSOR<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::CGSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential SSOR-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: SSOR symmetric successive overrelaxation method. The
 * relaxation is controlled by the parameter LinearSolver.PreconditionerRelaxation.
 * In each preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class SSORCGBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  SSORCGBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqSSOR<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::CGSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential GS-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: GS Gauss-Seidel method. It can be damped by the relaxation
 * parameter LinearSolver.PreconditionerRelaxation. In each preconditioning step,
 * it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class GSCGBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  GSCGBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqGS<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::CGSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential Jacobi-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: Jacobi method. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation. In each preconditioning step, it is
 * applied as often as given by the parameter LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class JacCGBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  JacCGBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqJac<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::CGSolver<Vector> Solver;

    return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential SSOR-preconditioned GMRes solver.
 *
 * Solver: The GMRes (generalized minimal residual) method is an iterative
 * method for the numerical solution of a nonsymmetric system of linear
 * equations.\n
 * See: Saad, Y., Schultz, M. H. (1986). "GMRES: A generalized minimal residual
 * algorithm for solving nonsymmetric linear systems." SIAM J. Sci. and Stat.
 * Comput. 7: 856–869.
 *
 * Preconditioner: SSOR symmetric successive overrelaxation method. The
 * relaxation is controlled by the parameter LinearSolver.PreconditionerRelaxation.
 * In each preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class SSORRestartedGMResBackend: public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  SSORRestartedGMResBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    typedef Dune::SeqSSOR<Matrix, Vector, Vector, blockLevel> Preconditioner;
    typedef Dune::RestartedGMResSolver<Vector> Solver;
    const int restart = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, GMResRestart);

    return ParentType::template solve<Preconditioner, Solver>(A, x, b, restart);
  }
};

/*!
 * \ingroup Linear
 * \brief Base class for backend combinations of linear solvers and a ILU0 preconditioner
 *
 * This class is used as a base class for combinations of a specific linear
 * solver with the ILU(0) preconditioner. Several parameters from the group
 * LinearSolver are read to customize the solver and preconditioner:
 *
 * Verbosity: determines how verbose the linear solver should print output.
 *
 * MaxIterations: the maximum number of iterations for the linear solver.
 *
 * ResidualReduction: the threshold for declaration of convergence.
 *
 * PreconditionerRelaxation: relaxation parameter for the preconditioner.
 */
template <class TypeTag>
class ILU0SolverBackend
{
public:

  template<class Preconditioner, class Solver, class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
    const int maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
    const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);

    Vector bTmp(b);

    const double relaxation = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, PreconditionerRelaxation);

    Preconditioner precond(A, relaxation);

    typedef Dune::MatrixAdapter<Matrix, Vector, Vector> MatrixAdapter;
    MatrixAdapter operatorA(A);

    Solver solver(operatorA, precond, residReduction, maxIter, verbosity);

    solver.apply(x, bTmp, result_);

    return result_.converged;
  }

  // solve with RestartedGMRes (needs restartGMRes as additional argument)
  template<class Preconditioner, class Solver, class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b, const int restartGMRes)
  {
    int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
    const int maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
    const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);

    Vector bTmp(b);

    const double relaxation = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, PreconditionerRelaxation);

    Preconditioner precond(A, relaxation);

    typedef Dune::MatrixAdapter<Matrix, Vector, Vector> MatrixAdapter;
    MatrixAdapter operatorA(A);

    Solver solver(operatorA, precond, residReduction, restartGMRes, maxIter, verbosity);

    solver.apply(x, bTmp, result_);

    return result_.converged;
  }

  const Dune::InverseOperatorResult& result() const
  {
    return result_;
  }

private:
  Dune::InverseOperatorResult result_;
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(0)-preconditioned BiCGSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. It can be applied to
 * nonsymmetric matrices.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: ILU(0) incomplete LU factorization. The order 0 indicates
 * that no fill-in is allowed. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class ILU0BiCGSTABBackend : public ILU0SolverBackend<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef ILU0SolverBackend<TypeTag> ParentType;
    enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
  public:

  ILU0BiCGSTABBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
      typedef Dune::SeqILU0<Matrix, Vector, Vector, blockLevel> Preconditioner;
      typedef Dune::BiCGSTABSolver<Vector> Solver;

      return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};



/*!
  \brief Schur complement inverse operator.
*/
template<class AType, class BType, class CType, class DType, class InvOpAType, class X, class Y>
class SchurComplement : public Dune::LinearOperator<X,Y>
{
public:
    // export types
    typedef DType matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

    //! constructor: just store a reference to a matrix
    explicit SchurComplement (const AType& A, const BType& B,
                              const CType& C, const DType& D,
                              const std::shared_ptr<InvOpAType>& invA)
    : A_(A), B_(B), C_(C), D_(D), invA_(invA) {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
        // apply B, note x and aTmp1 have different size (Bx)
        X aTmp1(A_.N());
        aTmp1 = 0.0;
        B_.mv(x, aTmp1);

        // apply A^-1 (A^-1Bx)
        auto aTmp2 = aTmp1;

        // if we use the exact Schur complement the eigenvalues are 1 and GMRes converges in 2 steps
        // ... a single V-cycle of AMG instead
        Dune::InverseOperatorResult result;
        invA_->apply(aTmp2, aTmp1, result);
        // invA_.pre(aTmp2, aTmp1);
        // invA_.apply(aTmp2, aTmp1);
        // invA_.post(aTmp2);

        // apply C (CA^-1Bx)
        C_.mv(aTmp2, y);
        // switch signs
        y *= -1.0;

        // add Dx (Dx - CA^-1Bx) -> the full Schur complement
        D_.umv(x, y);
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
        auto tmp = y;
        this->apply(x, tmp);
        y.axpy(alpha, tmp);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
      return Dune::SolverCategory::sequential;
    }

  private:
    const AType& A_;
    const BType& B_;
    const CType& C_;
    const DType& D_;
    std::shared_ptr<InvOpAType> invA_;
};

/*!
  \brief Schur complement inverse operator.
*/
template<class AType, class BType, class CType, class DType, class X, class Y>
class SchurApproximate : public Dune::LinearOperator<X,Y>
{
public:
    // export types
    typedef DType matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

    //! constructor: just store a reference to a matrix
    explicit SchurApproximate (const AType& A, const BType& B,
                               const CType& C, const DType& D)
    : A_(A), B_(B), C_(C), D_(D) {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
        auto viscosity = 1.0;
        auto theta = 1.0;
        auto rho = 1.0;

        // theta*rho*C*B*x
        X tmp(A_.N());
        tmp = 0.0;
        B_.mv(x, tmp);
        C_.mv(tmp, y);
        y *= theta*rho;

        // + D*x
        D_.umv(x, y);

        // + I*mu*x
        y.axpy(viscosity, x);
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
        auto tmp = y;
        this->apply(x, tmp);
        y.axpy(alpha, tmp);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
      return Dune::SolverCategory::sequential;
    }

  private:
    const AType& A_;
    const BType& B_;
    const CType& C_;
    const DType& D_;
};

/*! \brief The schur complement preconditioner

   \tparam M The matrix type to operate on
   \tparam X Type of the update
   \tparam Y Type of the defect
   \tparam l The block level to invert. Default is 1
*/
template<class TypeTag, class AInvOpType, class M, class X, class Y, int l=1>
class SchurComplementPreconditioner : public Dune::Preconditioner<X,Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \note this is an upper triangular preconditioner of the form
            | A  B |^-1
            | 0  S |
     */
    SchurComplementPreconditioner (const M& matrix,
                                   const std::shared_ptr<AInvOpType>& aInvOp)
      : M_(matrix), aInvOp_(aInvOp)
    {}

    virtual void pre (X& x, Y& b)
    {
      DUNE_UNUSED_PARAMETER(x);
      DUNE_UNUSED_PARAMETER(b);
    }

    /*!
       \brief Apply the preconditioner.

       \copydoc Preconditioner::apply(X&,const Y&)
     */
    virtual void apply (X& v, const Y& b)
    {
        using namespace Dune::Indices;
        using VVector = Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables)>;
        using PVector = Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables)>;

        using BType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockCCToFace;
        using AType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockFaceToFace;
        using CType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockFaceToCC;
        using DType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockCCToCC;

        // Currently Navier-Stokes system is assembled the other way around
        // Consider flipping indices
        const AType& A = M_[_1][_1];
        const BType& B = M_[_1][_0];
        const CType& C = M_[_0][_1];
        const DType& D = M_[_0][_0];

        const VVector& f = b[_1];
        const PVector& g = b[_0];

        // first solve pressure by applying the Schur inverse operator using a GMRes solver
        static const int schurMaxIter = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, SchurIterations);
        static const int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, PreconditionerVerbosity);
        static const double schurResred = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, SchurResidualReduction);
        static const double velResred = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, VelocityResidualReduction);
        using Identity = Dune::Richardson<PVector, PVector>;
        Identity identity(1.0);
        // SchurComplement<AType, BType, CType, DType, AInvOpType, PVector, PVector> schurOperator(A, B, C, D, aInvOp_);
        SchurApproximate<AType, BType, CType, DType, PVector, PVector> schurOperator(A, B, C, D);
        Dune::RestartedGMResSolver<PVector> invOpS(schurOperator, identity, schurResred*g.two_norm(), 100, schurMaxIter, verbosity);
        // Dune::BiCGSTABSolver<PVector> invOpS(schurOperator, identity, schurResred*g.two_norm(), schurMaxIter, verbosity);

        // apply S^-1 to solve pressure block
        PVector gTmp(g);
        Dune::InverseOperatorResult result;
        invOpS.apply(v[_0], gTmp, result);

        // prepare rhs for application of A^-1 (f - Bp)
        VVector vrhs(f);
        B.mmv(v[_0], vrhs);

        // then solve for velocity block v using the action of A^-1
        // approximate the inverse operator of A or use direct solver
        // here e.g. GMRes with AMG preconditioner
        Dune::InverseOperatorResult resultTmp;
        aInvOp_->apply(v[_1], vrhs, velResred*vrhs.two_norm(), resultTmp);
    }

    virtual void post (X& x)
    { DUNE_UNUSED_PARAMETER(x); }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    { return Dune::SolverCategory::sequential; }

  private:
    //! \brief The matrix we operate on.
    const M& M_;
    std::shared_ptr<AInvOpType> aInvOp_;
  };

template <class TypeTag>
class SchurComplementSolver
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

public:
    SchurComplementSolver() = default;
    SchurComplementSolver(const Problem& problem) {}

    // expects a system as a multi-type matrix
    // | A  B |
    // | C  D |

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        const int maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
        const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        const int restartGMRes = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, GMResRestart);

        Vector bTmp(b);

        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter operatorA(A);

        // LU decompose vel sub matrix and hand it to the preconditioner
        using namespace Dune::Indices;
        using VelMatrixType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockFaceToFace;
        using VVector = Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables)>;
        // using InnerSolver = Dune::UMFPack<VelMatrixType>;
        // auto invOpVel = std::make_shared<InnerSolver>(A[_1][_1], false);
        // using InnerPreconditioner = Dune::SeqJac<VelMatrixType, VVector, VVector, 1>;
        // auto innerPrec = std::make_shared<InnerPreconditioner>(A[_1][_1], 1, 1.0);

        // or use an AMG approximation instead instead
        using Comm = Dune::Amg::SequentialInformation;
        using LinearOperator = Dune::MatrixAdapter<VelMatrixType, VVector, VVector>;
        using ScalarProduct = Dune::SeqScalarProduct<VVector>;
        using Smoother = Dune::SeqSSOR<VelMatrixType, VVector, VVector>;
        auto comm = std::make_shared<Comm>();
        auto fop = std::make_shared<LinearOperator>(A[_1][_1]);
        auto sp = std::make_shared<ScalarProduct>();

        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<VelMatrixType, Dune::Amg::FirstDiagonal> >;
        Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(GET_PROP_TYPE(TypeTag, GridView)::Traits::Grid::dimension);
        params.setDebugLevel(verbosity);
        params.setGamma(1); // one V-cycle
        Criterion criterion(params);
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1.0;

        using InnerPreconditioner = Dune::Amg::AMG<LinearOperator, VVector, Smoother, Comm>;
        auto innerPrec = std::make_shared<InnerPreconditioner>(*fop, criterion, smootherArgs, *comm);
        static const int precVerbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, PreconditionerVerbosity);
        static const int precMaxIter = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, VelocityIterations);
        using InnerSolver = Dune::RestartedGMResSolver<VVector>;
        auto invOpVel = std::make_shared<InnerSolver>(*fop, *sp, *innerPrec, 1.0, 100, precMaxIter, precVerbosity);

        using Preconditioner = SchurComplementPreconditioner<TypeTag, InnerSolver, Matrix, Vector, Vector, 1>;
        Preconditioner schurPrecond(A, invOpVel);

        // using Solver = Dune::MINRESSolver<Vector>;
        using Solver = Dune::RestartedfGMResSolver<Vector>;
        // using Solver = Dune::BiCGSTABSolver<Vector>;
        Solver solver(operatorA, schurPrecond, residReduction, restartGMRes, maxIter, verbosity);
        // Solver solver(operatorA, schurPrecond, residReduction, maxIter, verbosity);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

private:
    Dune::InverseOperatorResult result_;
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(0)-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: ILU(0) incomplete LU factorization. The order 0 indicates
 * that no fill-in is allowed. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class ILU0CGBackend : public ILU0SolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef ILU0SolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  ILU0CGBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
      typedef Dune::SeqILU0<Matrix, Vector, Vector, blockLevel> Preconditioner;
      typedef Dune::CGSolver<Vector> Solver;

      return ParentType::template solve<Preconditioner, Solver>(A, x, b);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU0-preconditioned GMRes solver.
 *
 * Solver: The GMRes (generalized minimal residual) method is an iterative
 * method for the numerical solution of a nonsymmetric system of linear
 * equations.\n
 * See: Saad, Y., Schultz, M. H. (1986). "GMRES: A generalized minimal residual
 * algorithm for solving nonsymmetric linear systems." SIAM J. Sci. and Stat.
 * Comput. 7: 856–869.
 *
 * Preconditioner: ILU(0) incomplete LU factorization. The order 0 indicates
 * that no fill-in is allowed. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class ILU0RestartedGMResBackend : public ILU0SolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef ILU0SolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  ILU0RestartedGMResBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
      typedef Dune::SeqILU0<Matrix, Vector, Vector, blockLevel> Preconditioner;
      typedef Dune::RestartedGMResSolver<Vector> Solver;
      const int restart = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, GMResRestart);

      return ParentType::template solve<Preconditioner, Solver>(A, x, b, restart);
  }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(n)-preconditioned GMRes solver.
 *
 * Solver: The GMRes (generalized minimal residual) method is an iterative
 * method for the numerical solution of a nonsymmetric system of linear
 * equations.\n
 * See: Saad, Y., Schultz, M. H. (1986). "GMRES: A generalized minimal residual
 * algorithm for solving nonsymmetric linear systems." SIAM J. Sci. and Stat.
 * Comput. 7: 856–869.
 *
 * Preconditioner: ILU(n) incomplete LU factorization. The order n can be
 * provided by the parameter LinearSolver.PreconditionerIterations and controls
 * the fill-in. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template <class TypeTag>
class ILUnRestartedGMResBackend : public IterativePrecondSolverBackend<TypeTag>
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
  typedef IterativePrecondSolverBackend<TypeTag> ParentType;
  enum { blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel) };
public:

  ILUnRestartedGMResBackend(const Problem& problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
      typedef Dune::SeqILUn<Matrix, Vector, Vector, blockLevel> Preconditioner;
      typedef Dune::RestartedGMResSolver<Vector> Solver;
      const int restart = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, GMResRestart);

      return ParentType::template solve<Preconditioner, Solver>(A, x, b, restart);
  }
};

#if HAVE_SUPERLU
/*!
 * \ingroup Linear
 * \brief Direct linear solver using the SuperLU library.
 *
 * See: Li, X. S. (2005). "An overview of SuperLU: Algorithms, implementation,
 * and user interface." ACM Transactions on Mathematical Software (TOMS) 31(3): 302-325.
 * http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
 */
template <class TypeTag>
class SuperLUBackend
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:

  SuperLUBackend(const Problem& problem)
  : problem_(problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    Vector bTmp(b);

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum {blockSize = GET_PROP_VALUE(TypeTag, LinearSolverBlockSize)};
    typedef typename Dune::FieldMatrix<Scalar, blockSize, blockSize> MatrixBlock;
    typedef typename Dune::BCRSMatrix<MatrixBlock> ISTLMatrix;

    int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
    Dune::SuperLU<ISTLMatrix> solver(A, verbosity > 0);

    solver.apply(x, bTmp, result_);

    int size = x.size();
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < blockSize; j++)
        {
            using std::isnan;
            using std::isinf;
            if (isnan(x[i][j]) || isinf(x[i][j]))
            {
                result_.converged = false;
                break;
            }
        }
    }

    return result_.converged;
  }

  const Dune::InverseOperatorResult& result() const
  {
    return result_;
  }

private:
  Dune::InverseOperatorResult result_;
  const Problem& problem_;
};
#endif // HAVE_SUPERLU



#if HAVE_UMFPACK
/*!
 * \ingroup Linear
 * \brief Direct linear solver using the UMFPack library.
 *
 * See: Davis, Timothy A. (2004). "Algorithm 832". ACM Transactions on
 * Mathematical Software 30 (2): 196–199. doi:10.1145/992200.992206.
 * http://faculty.cse.tamu.edu/davis/suitesparse.html
 */
template <class TypeTag>
class UMFPackBackend
{
  typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:

  UMFPackBackend(const Problem& problem)
  : problem_(problem)
  {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    Vector bTmp(b);

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum {blockSize = GET_PROP_VALUE(TypeTag, LinearSolverBlockSize)};
    typedef typename Dune::FieldMatrix<Scalar, blockSize, blockSize> MatrixBlock;
    typedef typename Dune::BCRSMatrix<MatrixBlock> ISTLMatrix;

    int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
    Dune::UMFPack<ISTLMatrix> solver(A, verbosity > 0);

    solver.apply(x, bTmp, result_);

    int size = x.size();
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < blockSize; j++)
        {
            using std::isnan;
            using std::isinf;
            if (isnan(x[i][j]) || isinf(x[i][j]))
            {
                result_.converged = false;
                break;
            }
        }
    }

    return result_.converged;
  }

  const Dune::InverseOperatorResult& result() const
  {
    return result_;
  }

private:
  Dune::InverseOperatorResult result_;
  const Problem& problem_;
};
#endif // HAVE_UMFPACK

}
#endif
