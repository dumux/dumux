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

  ILUnBiCGSTABBackend()
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

  SORBiCGSTABBackend()
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

  SSORBiCGSTABBackend()
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

  GSBiCGSTABBackend()
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

  JacBiCGSTABBackend()
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

  ILUnCGBackend()
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

  SORCGBackend()
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

  SSORCGBackend()
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

  GSCGBackend()
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

  JacCGBackend()
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

  SSORRestartedGMResBackend()
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

  ILU0BiCGSTABBackend()
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

  ILU0CGBackend()
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

  ILU0RestartedGMResBackend()
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
  {}

  SuperLUBackend()
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
  {}

  UMFPackBackend()
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
};
#endif // HAVE_UMFPACK

}
#endif
