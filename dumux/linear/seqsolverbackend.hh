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
 * \brief Dumux sequential solver backends
 */
#ifndef DUMUX_SEQ_SOLVER_BACKEND_HH
#define DUMUX_SEQ_SOLVER_BACKEND_HH

#include <type_traits>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/linear/solver.hh>

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
class IterativePreconditionedSolverImpl
{
public:

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    static bool solve(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                      const std::string& modelParamGroup = "")
    {
        Preconditioner precond(A, s.precondIter(), s.relaxation());

        // make a linear operator from a matrix
        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter linearOperator(A);

        Solver solver(linearOperator, precond, s.residReduction(), s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    static bool solveWithGMRes(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                               const std::string& modelParamGroup = "")
    {
        // get the restart threshold
        const int restartGMRes = getParamFromGroup<double>(modelParamGroup, "LinearSolver.GMResRestart");

        Preconditioner precond(A, s.precondIter(), s.relaxation());

        // make a linear operator from a matrix
        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter linearOperator(A);

        Solver solver(linearOperator, precond, s.residReduction(), restartGMRes, s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    static bool solveWithILU0Prec(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                                  const std::string& modelParamGroup = "")
    {
        Preconditioner precond(A, s.relaxation());

        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter operatorA(A);

        Solver solver(operatorA, precond, s.residReduction(), s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    // solve with RestartedGMRes (needs restartGMRes as additional argument)
    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    static bool solveWithILU0PrecGMRes(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                                       const std::string& modelParamGroup = "")
    {
        // get the restart threshold
        const int restartGMRes = getParamFromGroup<int>(modelParamGroup, "LinearSolver.GMResRestart");

        Preconditioner precond(A, s.relaxation());

        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter operatorA(A);

        Solver solver(operatorA, precond, s.residReduction(), restartGMRes, s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }
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
class ILUnBiCGSTABBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqILUn<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "ILUn preconditioned BiCGSTAB solver";
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
class SORBiCGSTABBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqSOR<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "SOR preconditioned BiCGSTAB solver";
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
class SSORBiCGSTABBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqSSOR<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "SSOR preconditioned BiCGSTAB solver";
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
class GSBiCGSTABBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqGS<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "SSOR preconditioned BiCGSTAB solver";
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
class JacBiCGSTABBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqJac<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "Jac preconditioned BiCGSTAB solver";
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
class ILUnCGBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqILUn<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "ILUn preconditioned CG solver";
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
class SORCGBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqSOR<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "SOR preconditioned CG solver";
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
class SSORCGBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqSSOR<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "SSOR preconditioned CG solver";
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
class GSCGBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqGS<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "GS preconditioned CG solver";
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
class JacCGBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqJac<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "GS preconditioned CG solver";
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
class SSORRestartedGMResBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqSSOR<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::RestartedGMResSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithGMRes<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "SSOR preconditioned GMRes solver";
    }
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
class ILU0BiCGSTABBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqILU0<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithILU0Prec<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "ILU0 preconditioned BiCGSTAB solver";
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
class ILU0CGBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqILU0<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithILU0Prec<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "ILU0 preconditioned BiCGSTAB solver";
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
class ILU0RestartedGMResBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqILU0<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::RestartedGMResSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithILU0PrecGMRes<Preconditioner, Solver>(*this, A, x, b, paramGroup);
    }

    std::string name() const
    {
        return "ILU0 preconditioned BiCGSTAB solver";
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
class ILUnRestartedGMResBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockLevel = GET_PROP_VALUE(TypeTag, LinearSolverPreconditionerBlockLevel);
        using Preconditioner = Dune::SeqILUn<Matrix, Vector, Vector, blockLevel>;
        using Solver = Dune::RestartedGMResSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithGMRes<Preconditioner, Solver>(*this, A, x, b, paramGroup);
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
class SuperLUBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockSize = GET_PROP_VALUE(TypeTag, LinearSolverBlockSize);
        using MatrixBlock = typename Dune::FieldMatrix<double, blockSize, blockSize>;
        using ISTLMatrix = typename Dune::BCRSMatrix<MatrixBlock>;
        static_assert(std::is_same<Matrix, ISTLMatrix>::value, "SuperLU only works with BCRS matrices!");

        Dune::SuperLU<ISTLMatrix> solver(A, this->verbosity() > 0);

        Vector bTmp(b);
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

    std::string name() const
    {
        return "SuperLU solver";
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
class UMFPackBackend : public LinearSolver<TypeTag>
{
public:

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static const std::string paramGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);
        constexpr auto blockSize = GET_PROP_VALUE(TypeTag, LinearSolverBlockSize);
        using MatrixBlock = typename Dune::FieldMatrix<double, blockSize, blockSize>;
        using ISTLMatrix = typename Dune::BCRSMatrix<MatrixBlock>;
        static_assert(std::is_same<Matrix, ISTLMatrix>::value, "UMFPack only works with BCRS matrices!");

        Dune::UMFPack<ISTLMatrix> solver(A, this->verbosity() > 0);

        Vector bTmp(b);
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

    std::string name() const
    {
        return "UMFPack solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    Dune::InverseOperatorResult result_;
};
#endif // HAVE_UMFPACK

} // end namespace Dumux

#endif
