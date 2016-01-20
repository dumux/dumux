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
 * \brief Dumux solver backend the Pardiso direct solver.
 */
#ifndef DUMUX_PARDISO_BACKEND_HH
#define DUMUX_PARDISO_BACKEND_HH

#ifdef HAVE_PARDISO

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/implicit/box/properties.hh>
#include <dumux/decoupled/common/decoupledproperties.hh>



/* Change this, if your Fortran compiler does not append underscores. */
/* e.g. the AIX compiler:  #define F77_FUN(func) func                */
#ifdef AIX
#define F77_FUN(func) func
#else
#define F77_FUN(func) func ## _
#endif

/* PARDISO prototype. */
extern "C" int F77_FUN(pardisoinit)
    (void *pt, int *mtype, int *solver, int *iparm, double *dparm, int *error);

extern "C" int F77_FUN(pardiso)
    (void *pt, int *maxfct, int *mnum, int *mtype, int *phase, int *n,
     double *a, int *ia, int *ja, int *perm, int *nrhs, int *iparm,
     int *msglvl, double *b, double *x, int *error, double *dparm);



namespace Dumux
{
// forward declaration
template<class Matrix, class Vector> class Pardiso;


/*!
 * \ingroup Linear
 * \brief A solver backend for the Pardiso direct solver.
 *
 * See http://www.pardiso-project.org.
 */
template <class TypeTag>
class PardisoBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

public:

    PardisoBackend(const Problem& problem)
    {}

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        const int maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
        const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);

        Vector bTmp(b);
        Matrix ATemp(A);

        int numProcs = 1;
        if (ParameterTree::tree().hasKey("Pardiso.NumProcessors"))
            numProcs = std::max(numProcs, GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Pardiso, NumProcessors));

        Pardiso<Matrix, Vector> precond(ATemp, verbosity > 0, numProcs);

        typedef Dune::MatrixAdapter<Matrix, Vector, Vector> MatrixAdapter;
        MatrixAdapter operatorA(ATemp);

        Dune::BiCGSTABSolver<Vector> solver(operatorA, precond, residReduction, maxIter, verbosity);
        //  Dune::LoopSolver<Vector> solver(operatorA, precond, residReduction, maxIter, verbosity);

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
 * \brief A Pardiso preconditioner.
 *
 *     Put the Pardiso direct solver into the preconditioner framework.
 */
template<class Matrix, class Vector>
class Pardiso : public Dune::Preconditioner<Vector, Vector>
{
public:
    typedef typename Matrix::RowIterator RowIterator;
    typedef typename Matrix::ColIterator ColIterator;

    // define the category
    enum {
        //! \brief The category the preconditioner is part of
        category=Dune::SolverCategory::sequential
    };

    /*! \brief Constructor.

      Constructor gets all parameters to operate the preconditioner.
      \param A The matrix to operate on
      \param verbose Write messages if true
    */
    Pardiso (Matrix& A, bool verbose = false, int numProcs = 1)
        : verbose_(verbose), num_procs_(numProcs)
    {
        initAndFactor(A);
    }

    ~Pardiso()
    {
        if (phase_ != -1)
        {
            phase_ = -1;                 // Release internal memory.
            int idum;
            double ddum;

            F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                              &n_, &ddum, ia_, ja_, &idum, &nrhs_,
                              iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);
            delete a_;
            delete ia_;
            delete ja_;
        }
    }

    /*!
     * \brief Initialize and factorize the matrix
     *
     * \param A the matrix
     */
    void initAndFactor(Matrix& A)
    {
        mtype_ = 11;
        nrhs_ = 1;
        maxfct_ = 1;
        mnum_ = 1;
        msglvl_ = 0;
        error_ = 0;
        solver_ = 0; // solver_ = 0, choose sparse direct solver, = 1 multi-recursive iterative solver

        RowIterator i0 = A.begin();
        ColIterator j0 = (*i0).begin();

        systemsize_ = (*j0).N();
        n_ = A.N()*systemsize_;

        int nnz = 0;
        RowIterator endi = A.end();
        int rows = 0;
        for (RowIterator i = A.begin(); i != endi; ++i)
        {
            rows++;
            ColIterator endj = (*i).end();
            for (ColIterator j = (*i).begin(); j != endj; ++j) {
                nnz += systemsize_*systemsize_;
            }
        }
        if (verbose_)
            std::cout << "Pardiso: dimension = " << n_ << ", number of nonzeros = " << nnz << std::endl;

        a_ = new double[nnz];
        ia_ = new int[n_+1];
        ja_ = new int[nnz];

        int count = 0;
        for (RowIterator i = A.begin(); i != endi; ++i)
        {
            for (int iComp = 0; iComp < systemsize_; iComp++) {
                ia_[i.index()*systemsize_ + iComp] = count+1;
                ColIterator endj = (*i).end();
                for (ColIterator j = (*i).begin(); j != endj; ++j) {
                    for (int jComp = 0; jComp < systemsize_; jComp++) {
                        a_[count] = (*j)[iComp][jComp];
                        ja_[count] = j.index()*systemsize_ + jComp + 1;

                        count++;
                    }
                }
            }
        }
        ia_[n_] = count+1;

        F77_FUN(pardisoinit) (pt_,  &mtype_, &solver_, iparm_, dparm_, &error_);

        if (error_)
        {
            switch(error_)
            {
            case -10:
                DUNE_THROW(Dune::MathError, "Constructor Pardiso: pardisoinit failed. No license file found. Error code " << error_);
                break;
            case -11:
                DUNE_THROW(Dune::MathError, "Constructor Pardiso: pardisoinit failed. License has expired. Error code " << error_);
                break;
            case -12:
                DUNE_THROW(Dune::MathError, "Constructor Pardiso: pardisoinit failed. Wrong username or hostname. Error code " << error_);
                break;
            default:
                DUNE_THROW(Dune::MathError, "Constructor Pardiso: pardisoinit failed. Error code " << error_);
                break;
            }
        }

        phase_ = 11;
        int idum;
        double ddum;
        iparm_[2]  = num_procs_;

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                          &n_, a_, ia_, ja_, &idum, &nrhs_,
                          iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);

        if (error_ != 0)
            DUNE_THROW(Dune::MathError, "Constructor Pardiso: Reordering failed. Error code " << error_);

        if (verbose_)
            std::cout << "  Reordering completed. Number of nonzeros in factors = " << iparm_[17] << std::endl;

        phase_ = 22;

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                          &n_, a_, ia_, ja_, &idum, &nrhs_,
                          iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);

        if (error_ != 0)
            DUNE_THROW(Dune::MathError, "Constructor Pardiso: Factorization failed. Error code " << error_);

        if (verbose_)
            std::cout << "  Factorization completed." << std::endl;
    }

    /*!
     * \brief Factorize the matrix.
     *
     * \param A The matrix
     */
    void factorize (Matrix& A)
    {
        RowIterator i0 = A.begin();
        ColIterator j0 = (*i0).begin();

        systemsize_ = (*j0).N();
        n_ = A.N()*systemsize_;
        int nnz = 0;
        RowIterator endi = A.end();
        int rows = 0;
        for (RowIterator i = A.begin(); i != endi; ++i)
        {
            rows++;
            ColIterator endj = (*i).end();
            for (ColIterator j = (*i).begin(); j != endj; ++j) {
                nnz += systemsize_*systemsize_;
            }
        }
        if (verbose_)
            std::cout << "Pardiso: dimension = " << n_ << ", number of nonzeros = " << nnz << std::endl;

        a_ = new double[nnz];
        ia_ = new int[n_+1];
        ja_ = new int[nnz];

        int count = 0;
        for (RowIterator i = A.begin(); i != endi; ++i)
        {
            for (int iComp = 0; iComp < systemsize_; iComp++) {
                ia_[i.index()*systemsize_ + iComp] = count+1;
                ColIterator endj = (*i).end();
                for (ColIterator j = (*i).begin(); j != endj; ++j) {
                    for (int jComp = 0; jComp < systemsize_; jComp++) {
                        a_[count] = (*j)[iComp][jComp];
                        ja_[count] = j.index()*systemsize_ + jComp + 1;

                        count++;
                    }
                }
            }
        }
        ia_[n_] = count+1;

        F77_FUN(pardisoinit) (pt_,  &mtype_, &solver_, iparm_, dparm_, &error_);

        if (error_)
        {
            switch(error_)
            {
            case -10:
                DUNE_THROW(Dune::MathError, "Constructor Pardiso: pardisoinit failed. No license file found. Error code " << error_);
                break;
            case -11:
                DUNE_THROW(Dune::MathError, "Constructor Pardiso: pardisoinit failed. License has expired. Error code " << error_);
                break;
            case -12:
                DUNE_THROW(Dune::MathError, "Constructor Pardiso: pardisoinit failed. Wrong username or hostname. Error code " << error_);
                break;
            default:
                DUNE_THROW(Dune::MathError, "Constructor Pardiso: pardisoinit failed. Error code " << error_);
                break;
            }
        }


        phase_ = 11;
        int idum;
        double ddum;
        iparm_[2]  = num_procs_;

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                          &n_, a_, ia_, ja_, &idum, &nrhs_,
                          iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);

        if (error_ != 0)
            DUNE_THROW(Dune::MathError, "Constructor Pardiso: Reordering failed. Error code " << error_);

        if (verbose_)
            std::cout << "  Reordering completed. Number of nonzeros in factors = " << iparm_[17] << std::endl;

        phase_ = 22;

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                          &n_, a_, ia_, ja_, idum, &nrhs_,
                          iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);

        if (error_ != 0)
            DUNE_THROW(Dune::MathError, "Constructor Pardiso: Factorization failed. Error code " << error_);

        if (verbose_)
            std::cout << "  Factorization completed." << std::endl;
    }

    /*!
     * \brief Prepare the preconditioner.
     *
     * A solver solves a linear operator equation A(x)=b by applying
     * one or several steps of the preconditioner. The method pre()
     * is called before the first apply operation.
     * b and x are right hand side and solution vector of the linear
     * system respectively. It may. e.g., scale the system, allocate memory or
     * compute a (I)LU decomposition.
     * Note: The ILU decomposition could also be computed in the constructor
     * or with a separate method of the derived method if several
     * linear systems with the same matrix are to be solved.
     *
     * \param x The left hand side of the equation.
     * \param b The right hand side of the equation.
     */
    virtual void pre (Vector& x, Vector& b)
    {}

    /*! \brief Apply one step of the preconditioner to the system \f$ A(v)=d \f$.
     *
     * On entry \f$ v=0 \f$ and \f$ d=b-A(x) \f$ (although this might not be
     * computed in that way. On exit v contains the update, i.e
     * one step computes \f$ v = M^{-1} d \f$ where \f$ M \f$ is the
     * approximate inverse of the operator \f$ A \f$ characterizing
     * the preconditioner.
     * \param[out] v The update to be computed
     * \param d The current defect.
     */
    virtual void apply (Vector& v, const Vector& d)
    {
        phase_ = 33;

        iparm_[7] = 1;       /* Max numbers of iterative refinement steps. */
        int idum;
        double x[n_];
        double b[n_];
        for (typename Vector::size_type i = 0; i < v.size(); i++) {
            for (int comp = 0; comp < systemsize_; comp++) {
                x[i*systemsize_ + comp] = v[i][comp];
                b[i*systemsize_ + comp] = d[i][comp];
            }
        }

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                          &n_, a_, ia_, ja_, &idum, &nrhs_,
                          //&n_, &ddum, &idum, &idum, &idum, &nrhs_,
                          iparm_, &msglvl_, b, x, &error_, dparm_);

        if (error_ != 0)
            DUNE_THROW(Dune::MathError, "Pardiso.apply: Backsolve failed. Error code " << error_);

        for (typename Vector::size_type i = 0; i < v.size(); i++)
            for (int comp = 0; comp < systemsize_; comp++)
                v[i][comp] = x[i*systemsize_ + comp];
    }

    /*! \brief Clean up.
     *
     * This method is called after the last apply call for the
     * linear system to be solved. Memory may be deallocated safely
     * here. x is the solution of the linear equation.
     *
     * \param b The right hand side of the equation.
    */
    virtual void post (Vector& b)
    {
        phase_ = -1;                 // Release internal memory.
        int idum;
        double ddum;

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                          &n_, &ddum, ia_, ja_, &idum, &nrhs_,
                          iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);
        delete a_;
        delete ia_;
        delete ja_;
    }

private:
    int n_; //!< dimension of the system
    double *a_; //!< matrix values
    int *ia_; //!< indices to rows
    int *ja_; //!< column indices
    int mtype_; //!< matrix type, currently only 11 (real unsymmetric matrix) is supported
    int nrhs_; //!< number of right hand sides
    void *pt_[64]; //!< internal solver memory pointer
    int iparm_[64]; //!< Pardiso control parameters.
    int maxfct_;    //!< Maximum number of numerical factorizations.
    int mnum_;  //!<        Which factorization to use.
    int msglvl_;    //!< flag to print statistical information
    int error_;      //!< error flag
    bool verbose_;
    int num_procs_; //!< number of processors.
    int systemsize_;
    int phase_;
    double dparm_[64];
    int solver_;
};

}
#else
#warning "no Pardiso library available, reconfigure with correct --with-pardiso options"
#endif // HAVE_PARDISO

#endif // PARDISO_BACKEND
