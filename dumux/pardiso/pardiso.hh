#ifndef DUNE_PARDISO_HH
#define DUNE_PARDISO_HH

#include <dune/istl/preconditioners.hh>

/* Change this, if your Fortran compiler does not append underscores. */
/* e.g. the AIX compiler:  #define F77_FUNC(func) func                */

#ifdef AIX
#define F77_FUNC(func)  func 
#else
#define F77_FUNC(func)  func ## _
#endif


#ifdef HAVE_PARDISO 
/* PARDISO prototype. */
extern "C" int F77_FUNC(pardisoinit)
    (void *, int *, int *);

extern "C" int F77_FUNC(pardiso)
    (void *, int *, int *, int *, int *, int *, 
     double *, int *, int *, int *, int *, int *, 
     int *, double *, double *, int *);
#endif 

namespace Dune {


  /*! \brief The sequential Pardiso preconditioner.
    
     Put the Pardiso direct solver into the preconditioner framework.
   */
  template<class M, class X, class Y>
  class SeqPardiso : public Preconditioner<X,Y> {
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

    /*! \brief Constructor.

    Constructor gets all parameters to operate the prec.
    \param A The matrix to operate on.
    \param n The number of iterations to perform.
    \param w The relaxation factor.
    */
    SeqPardiso (const M& A)
      : A_(A)
    {
#ifdef HAVE_PARDISO
    	
    	mtype_ = 11;
    	nrhs_ = 1;
    	num_procs_ = 1;
        maxfct_ = 1;	
        mnum_   = 1;         
        msglvl_ = 0;        
        error_  = 0;        

    	n_ = A_.rowdim();
    	int nnz = 0;
    	RowIterator endi = A_.end();
    	for (RowIterator i = A_.begin(); i != endi; ++i)
    	{
    		if (A_.rowdim(i.index()) != 1)
    			DUNE_THROW(NotImplemented, "SeqPardiso: row blocksize != 1.");
    		ColIterator endj = (*i).end();
    		for (ColIterator j = (*i).begin(); j != endj; ++j) {
    			if (A_.coldim(j.index()) != 1)
    				DUNE_THROW(NotImplemented, "SeqPardiso: column blocksize != 1.");
    			nnz++;
    		}
    	}
		  
    	std::cout << "dimension = " << n_ << ", number of nonzeros = " << nnz << std::endl;
    	
    	a_ = new double[nnz];
    	ia_ = new int[n_+1];
    	ja_ = new int[nnz];
    	
    	int count = 0;
    	for (RowIterator i = A_.begin(); i != endi; ++i)
    	{
    		ia_[i.index()] = count+1;
    		ColIterator endj = (*i).end();
    		for (ColIterator j = (*i).begin(); j != endj; ++j) {
    			a_[count] = *j;
    			ja_[count] = j.index()+1;
    			
    			count++;
    		}
    	}
    	ia_[n_] = count+1;
    	
        F77_FUNC(pardisoinit) (pt_,  &mtype_, iparm_); 

        int phase = 11;
        int idum;
        double ddum;
        iparm_[2]  = num_procs_;
        
        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase,
    		       &n_, a_, ia_, ja_, &idum, &nrhs_,
    		       iparm_, &msglvl_, &ddum, &ddum, &error_);
      
        if (error_ != 0) 
        	DUNE_THROW(MathError, "Constructor SeqPardiso: Factorization failed. Error code " << error_);
        
        std::cout << "Constructor SeqPardiso: Factorization completed." << std::endl;
        
#else
        DUNE_THROW(NotImplemented, "no Pardiso library available, reconfigure with correct --with-pardiso options");
#endif
    }

    /*!
      \brief Prepare the preconditioner.
      
      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (X& x, Y& b) {}

    /*!
      \brief Apply the preconditioner.
      
      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (X& v, const Y& d)
    {
#ifdef HAVE_PARDISO
    	int phase = 33;

        iparm_[7] = 1;       /* Max numbers of iterative refinement steps. */
        int idum;

        double x[n_];
        double b[n_];
        for (int i = 0; i < n_; i++) {
        	x[i] = v[i];
        	b[i] = d[i];
        }
        
        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase,
    		       &n_, a_, ia_, ja_, &idum, &nrhs_,
    		       iparm_, &msglvl_, b, x, &error_);
        
        if (error_ != 0) 
        	DUNE_THROW(MathError, "SeqPardiso.apply: Backsolve failed. Error code " << error_);
        
        for (int i = 0; i < n_; i++) 
        	v[i] = x[i];
        
        std::cout << "SeqPardiso: Backsolve completed." << std::endl;
#endif
    }

    /*!
      \brief Clean up.
      
      \copydoc Preconditioner::post(X&)
    */
    virtual void post (X& x) {}

    ~SeqPardiso() 
    {
#ifdef HAVE_PARDISO
        int phase = -1;                 // Release internal memory. 
        int idum;
        double ddum;

        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase,
    		       &n_, &ddum, ia_, ja_, &idum, &nrhs_,
    		       iparm_, &msglvl_, &ddum, &ddum, &error_);
#endif
    }
  
  private:
    M A_; //!< The matrix we operate on.
    int n_; //!< dimension of the system
    double *a_; //!< matrix values 
    int *ia_; //!< indices to rows 
    int *ja_; //!< column indices
    int mtype_; //!< matrix type, currently only 11 (real unsymmetric matrix) is supported
    int nrhs_; //!< number of right hand sides
    void *pt_[64]; //!< internal solver memory pointer
    int iparm_[64]; //!< Pardiso control parameters.
    int maxfct_;	//!< Maximum number of numerical factorizations. 
    int mnum_;  //!<        Which factorization to use. 
    int msglvl_;    //!< flag to print statistical information 
    int error_;      //!< error flag 
    int num_procs_; //!< number of processors.
  };

}






#endif

