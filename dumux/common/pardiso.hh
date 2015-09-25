// $Id: pardiso.hh 3826 2010-07-14 07:03:41Z bernd $
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Sequential preconditioner using the Pardiso direct solver.
 */
#ifndef DUMUX_PARDISO_HH
#define DUMUX_PARDISO_HH

#include <dune/istl/preconditioners.hh>

/* Change this, if your Fortran compiler does not append underscores. */
/* e.g. the AIX compiler:  #define F77_FUN(func) func                */


#ifdef HAVE_PARDISO

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

#endif /* HAVE_PARDISO */

namespace Dune {


/*! \brief The sequential Pardiso preconditioner.

  Put the Pardiso direct solver into the preconditioner framework.
*/
template<class M, class X, class Y>
class SeqPardiso : public Dune::Preconditioner<X,Y> {
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

    typedef typename M::RowIterator RowIterator;
    typedef typename M::ColIterator ColIterator;

    // define the category
    enum {
        //! \brief The category the preconditioner is part of
        category=SolverCategory::sequential
    };


    SeqPardiso(bool verbose = false, int n = 0, double w = 0.0)
    : verbose_(verbose)
    {

#ifdef HAVE_PARDISO
        mtype_ = 11;
        nrhs_ = 1;
        num_procs_ = 1;
        maxfct_ = 1;
        mnum_ = 1;
        msglvl_ = 0;
        error_ = 0;

        solver_ = 0; // solver_ = 0, choose sparse direct solver, = 1 multi-recursive iterative solver


        //F77_FUN(pardisoinit) (pt_, &mtype_, iparm_);
#else
        DUNE_THROW(Dune::NotImplemented, "no Pardiso library available, reconfigure with correct --with-pardiso options");
#endif
    }

    void factorize (M& A)
    {
#ifdef HAVE_PARDISO

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
            //if (A.rowdim(i.index()) != 1)
            //    DUNE_THROW(Dune::NotImplemented, "SeqPardiso: row blocksize != 1.");
            ColIterator endj = (*i).end();
            for (ColIterator j = (*i).begin(); j != endj; ++j) {
                //if (A.coldim(j.index()) != 1)
                //    DUNE_THROW(Dune::NotImplemented, "SeqPardiso: column blocksize != 1.");
                nnz += systemsize_*systemsize_;
            }
        }
        //std::cout << "rows = " << rows;
        if (verbose_)
            std::cout << "SeqPardiso: dimension = " << n_ << ", number of nonzeros = " << nnz << std::endl;

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

        /*std::cout << "systemsize_ =" << systemsize_ << ", n_ = " << n_ << ", nnz_ = " << nnz << std::endl;
          for (int i = 0; i <= n_; i++)
          std::cout << ia_[i] << std::endl;
        */
        /*std::cout << "ja_:" << std::endl;
          for (int i = 0; i <= nnz; i++)
          std::cout << ja_[i] << std::endl;
          std::cout << "a_:" << std::endl;
          for (int i = 0; i <= nnz; i++)
          std::cout << a_[i] << std::endl;
        */

        F77_FUN(pardisoinit) (pt_,  &mtype_, &solver_, iparm_, dparm_, &error_);

        if (error_)
	{
		switch(error_)
		{
		    case -10:
			DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: pardisoinit failed. No license file found. Error code " << error_);
			break;
		    case -11:
			DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: pardisoinit failed. License has expired. Error code " << error_);
			break;
		    case -12:
			DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: pardisoinit failed. Wrong username or hostname. Error code " << error_);
			break;
		    default:
			DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: pardisoinit failed. Error code " << error_);
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
            DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: Reordering failed. Error code " << error_);

        if (verbose_)
            std::cout << "  Reordering completed. Number of nonzeros in factors = " << iparm_[17] << std::endl;

        phase_ = 22;

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                         &n_, a_, ia_, ja_, idum, &nrhs_,
                         iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);

        if (error_ != 0)
            DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: Factorization failed. Error code " << error_);

        if (verbose_)
            std::cout << "  Factorization completed." << std::endl;
#endif
    }

    void initAndFactor(M& A)
    {
#ifdef HAVE_PARDISO

        mtype_ = 11;
        nrhs_ = 1;
        num_procs_ = 1;
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
            //if (A.rowdim(i.index()) != 1)
            //    DUNE_THROW(Dune::NotImplemented, "SeqPardiso: row blocksize != 1.");
            ColIterator endj = (*i).end();
            for (ColIterator j = (*i).begin(); j != endj; ++j) {
                //if (A.coldim(j.index()) != 1)
                //    DUNE_THROW(Dune::NotImplemented, "SeqPardiso: column blocksize != 1.");
                nnz += systemsize_*systemsize_;
            }
        }
        //std::cout << "rows = " << rows;
        if (verbose_)
            std::cout << "SeqPardiso: dimension = " << n_ << ", number of nonzeros = " << nnz << std::endl;

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

        /*std::cout << "systemsize_ =" << systemsize_ << ", n_ = " << n_ << ", nnz_ = " << nnz << std::endl;
          for (int i = 0; i <= n_; i++)
          std::cout << ia_[i] << std::endl;
        */
        /*std::cout << "ja_:" << std::endl;
          for (int i = 0; i <= nnz; i++)
          std::cout << ja_[i] << std::endl;
          std::cout << "a_:" << std::endl;
          for (int i = 0; i <= nnz; i++)
          std::cout << a_[i] << std::endl;
        */

        F77_FUN(pardisoinit) (pt_,  &mtype_, &solver_, iparm_, dparm_, &error_);

        if (error_)
	{
		switch(error_)
		{
		    case -10:
			DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: pardisoinit failed. No license file found. Error code " << error_);
			break;
		    case -11:
			DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: pardisoinit failed. License has expired. Error code " << error_);
			break;
		    case -12:
			DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: pardisoinit failed. Wrong username or hostname. Error code " << error_);
			break;
		    default:
			DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: pardisoinit failed. Error code " << error_);
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
            DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: Reordering failed. Error code " << error_);

        if (verbose_)
            std::cout << "  Reordering completed. Number of nonzeros in factors = " << iparm_[17] << std::endl;

        phase_ = 22;

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                           &n_, a_, ia_, ja_, &idum, &nrhs_,
                           iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);

        if (error_ != 0)
            DUNE_THROW(Dune::MathError, "Constructor SeqPardiso: Factorization failed. Error code " << error_);

        if (verbose_)
            std::cout << "  Factorization completed." << std::endl;

#else
        DUNE_THROW(Dune::NotImplemented, "no Pardiso library available, reconfigure with correct --with-pardiso options");
#endif
    }

    /*! \brief Constructor.

      Constructor gets all parameters to operate the preconditioner.
      \param A The matrix to operate on
      \param verbose Write messages if true
    */
    SeqPardiso (M& A, bool verbose = false)
    : verbose_(verbose)
    {
        initAndFactor(A);
    }

    SeqPardiso (M& A, int iter, double relaxation)
    : verbose_(false)
    {
        initAndFactor(A);
    }

	/*! \brief Prepare the preconditioner.

	A solver solves a linear operator equation A(x)=b by applying
    one or several steps of the preconditioner. The method pre()
    is called before the first apply operation.
    b and x are right hand side and solution vector of the linear
    system respectively. It may. e.g., scale the system, allocate memory or
    compute a (I)LU decomposition.
	  Note: The ILU decomposition could also be computed in the constructor
    or with a separate method of the derived method if several
    linear systems with the same matrix are to be solved.

    \param x The left hand side of the equation.
    \param b The right hand side of the equation.
	*/
    virtual void pre (X& x, Y& b) {}

	/*! \brief Apply one step of the preconditioner to the system \f$ A(v)=d \f$.

    On entry \f$ v=0 \f$ and \f$ d=b-A(x) \f$ (although this might not be
    computed in that way. On exit v contains the update, i.e
    one step computes \f$ v = M^{-1} d \f$ where \f$ M \f$ is the
    approximate inverse of the operator \f$ A \f$ characterizing
    the preconditioner.
    \param[out] v The update to be computed
    \param d The current defect.
	*/
    virtual void apply (X& v, const Y& d)
    {
#ifdef HAVE_PARDISO
        phase_ = 33;

        iparm_[7] = 1;       /* Max numbers of iterative refinement steps. */
        int idum;
        double x[2*n_];
        double b[2*n_];
        for (typename X::size_type i = 0; i < v.size(); i++) {
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
            DUNE_THROW(Dune::MathError, "SeqPardiso.apply: Backsolve failed. Error code " << error_);

        for (typename X::size_type i = 0; i < v.size(); i++)
            for (int comp = 0; comp < systemsize_; comp++)
                v[i][comp] = x[i*systemsize_ + comp];


        //phase_ = -1;                 // Release internal memory.
        //int idum;
        //double ddum;

        //F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
        //           &n_, &ddum, ia_, ja_, &idum, &nrhs_,
        //           iparm_, &msglvl_, &ddum, &ddum, &error_);
        //delete a_;
        //delete ia_;
        //delete ja_;
        //std::cout << "SeqPardiso: Backsolve completed." << std::endl;
#endif
    }

	/*! \brief Clean up.

	This method is called after the last apply call for the
    linear system to be solved. Memory may be deallocated safely
    here. x is the solution of the linear equation.

    \param x The right hand side of the equation.
	*/
    virtual void post (X& x)
    {
#ifdef HAVE_PARDISO
        phase_ = -1;                 // Release internal memory.
        int idum;
        double ddum;

        F77_FUN(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                           &n_, &ddum, ia_, ja_, &idum, &nrhs_,
                           iparm_, &msglvl_, &ddum, &ddum, &error_, dparm_);
        delete a_;
        delete ia_;
        delete ja_;
#endif
    }

    ~SeqPardiso()
    {
#ifdef HAVE_PARDISO
        if (phase_ != -1) {
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
#endif
    }

private:
    //M A_; //!< The matrix we operate on.
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
    int num_procs_; //!< number of processors.
    int systemsize_;
    int phase_;
    bool verbose_;
    double dparm_[64];
    int solver_;
//    int perm_[64];  // not yet used here
};

}






#endif

