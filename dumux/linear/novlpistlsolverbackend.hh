// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
/*****************************************************************************
 *   This file is based on the file of DUNE-PDELab with the same name,       *
 *   http://www.dune-project.org/pdelab/.                                    *
 *   Some changes are made such that the backends work for the vectors and   *
 *   matrices used in Dumux.                                                 *
 *****************************************************************************/
#ifndef DUMUX_NOVLPISTLSOLVERBACKEND_HH
#define DUMUX_NOVLPISTLSOLVERBACKEND_HH

#include <cstddef>

#include <dune/common/deprecated.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/backend/parallelistlhelper.hh>

#include "istlvectorbackend.hh"
#include "seqistlsolverbackend.hh"

namespace Dumux {

    //! \addtogroup Backend
    //! \ingroup PDELab
    //! \{

    //========================================================
    // Generic support for nonoverlapping grids
    //========================================================

    //! Operator for the non-overlapping parallel case
    /**
     * Calculate \f$y:=Ax\f$.
     *
     * \tparam GFS The GridFunctionSpace the vectors apply to.
     * \tparam M   Type of the matrix.  Should be one of the ISTL matrix types.
     * \tparam X   Type of the vectors the matrix is applied to.
     * \tparam Y   Type of the result vectors.
     */
    template<class GFS, class M, class X, class Y>
    class NonoverlappingOperator : public Dune::AssembledLinearOperator<M,X,Y>
    {
    public:
      //! export type of matrix
      typedef M matrix_type;
      //! export type of vectors the matrix is applied to
      typedef X domain_type;
      //! export type of result vectors
      typedef Y range_type;
      //! export type of the entries for x
      typedef typename X::field_type field_type;

      //redefine the category, that is the only difference
      enum {category=Dune::SolverCategory::nonoverlapping};

      //! Construct a non-overlapping operator
      /**
       * \param gfs_    GridFunctionsSpace for the vectors.
       * \param A       Matrix for this operator.  This should be the locally
       *                assembled matrix.
       * \param helper_ Helper for parallel communication (not used).
       *
       * \note The constructed object stores references to all the objects
       *       given as parameters here.  They should be valid for as long as
       *       the constructed object is used.  They are not needed to
       *       destruct the constructed object.
       *
       * \deprecated The helper_ parameter is unused.  Use the constructor
       *             without the helper_ parameter instead.
       */
      NonoverlappingOperator (const GFS& gfs_, const M& A,
                              const Dune::PDELab::ParallelISTLHelper<GFS>& helper_)
        DUNE_DEPRECATED
        : gfs(gfs_), _A_(A)
      {
      }

      //! Construct a non-overlapping operator
      /**
       * \param gfs_ GridFunctionsSpace for the vectors.
       * \param A    Matrix for this operator.  This should be the locally
       *             assembled matrix.
       *
       * \note The constructed object stores references to all the objects
       *       given as parameters here.  They should be valid for as long as
       *       the constructed object is used.  They are not needed to
       *       destruct the constructed object.
       */
      NonoverlappingOperator (const GFS& gfs_, const M& A)
        : gfs(gfs_), _A_(A)
      { }

      //! apply operator
      /**
       * Compute \f$y:=A(x)\f$ on this process, then make y consistent (sum up
       * corresponding entries of y on the different processes and store the
       * result back in y on each process).
       */
      virtual void apply (const X& x, Y& y) const
      {
        // apply local operator; now we have sum y_p = sequential y
        _A_.mv(x,y);

        // accumulate y on border
        Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
        if (gfs.gridView().comm().size()>1)
          gfs.gridView().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
      /**
       * Compute \f$y:=\alpha A(x)\f$ on this process, then make y consistent
       * (sum up corresponding entries of y on the different processes and
       * store the result back in y on each process).
       */
      virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
      {
        // apply local operator; now we have sum y_p = sequential y
        _A_.usmv(alpha,x,y);

        // accumulate y on border
        Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
        if (gfs.gridView().comm().size()>1)
          gfs.gridView().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

      //! extract the matrix
      virtual const M& getmat () const
      {
        return _A_;
      }

    private:
      const GFS& gfs;
      const M& _A_;
    };

    // parallel scalar product assuming no overlap
    template<class GFS, class X>
    class NonoverlappingScalarProduct : public Dune::ScalarProduct<X>
    {
    public:
      //! export types
      typedef X domain_type;
      typedef typename X::ElementType field_type;

      //! define the category
      enum {category=Dune::SolverCategory::nonoverlapping};

      /*! \brief Constructor needs to know the grid function space
       */
      NonoverlappingScalarProduct (const GFS& gfs_, const Dune::PDELab::ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {}

      /*! \brief Dot product of two vectors.
        It is assumed that the vectors are consistent on the interior+border
        partition.
      */
      virtual field_type dot (const X& x, const X& y)
      {
        // do local scalar product on unique partition
        field_type sum = 0;
        for (typename X::size_type i=0; i<x.base().N(); ++i)
          for (typename X::size_type j=0; j<x[i].N(); ++j)
            sum += (x[i][j]*y[i][j])*helper.mask(i,j);

        // do global communication
        return gfs.gridView().comm().sum(sum);
      }

      /*! \brief Norm of a right-hand side vector.
        The vector must be consistent on the interior+border partition
      */
      virtual double norm (const X& x)
      {
        return sqrt(static_cast<double>(this->dot(x,x)));
      }

      /*! \brief make additive vector consistent
       */
      void make_consistent (X& x) const
      {
        Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,x);
        if (gfs.gridView().comm().size()>1)
          gfs.gridView().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
      }

    private:
      const GFS& gfs;
      const Dune::PDELab::ParallelISTLHelper<GFS>& helper;
    };

    // parallel Richardson preconditioner
    template<class GFS, class X, class Y>
    class NonoverlappingRichardson : public Dune::Preconditioner<X,Y>
    {
    public:
      //! \brief The domain type of the preconditioner.
      typedef X domain_type;
      //! \brief The range type of the preconditioner.
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;

      // define the category
      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::nonoverlapping
      };

      //! \brief Constructor.
      NonoverlappingRichardson (const GFS& gfs_, const Dune::PDELab::ParallelISTLHelper<GFS>& helper_)
        : gfs(gfs_), helper(helper_)
      {
      }

      /*!
        \brief Prepare the preconditioner.
      */
      virtual void pre (X& x, Y& b) {}

      /*!
        \brief Apply the precondioner.
      */
      virtual void apply (X& v, const Y& d)
      {
        v = d;
      }

      /*!
        \brief Clean up.
      */
      virtual void post (X& x) {}

    private:
      const GFS& gfs;
      const Dune::PDELab::ParallelISTLHelper<GFS>& helper;
    };

    //! parallel non-overlapping Jacobi preconditioner
    /**
     * \tparam Diagonal Vector type used to store the diagonal of the matrix
     * \tparam X        Vector type used to store the result of applying the
     *                  preconditioner.
     * \tparam Y        Vector type used to store the defect.
     *
     * The Jacobi preconditioner approximates the inverse of a matrix M by
     * taking the diagonal diag(M) and inverting that.  In the parallel case
     * the matrix M is assumed to be inconsistent, so diagonal entries for
     * dofs on the border are summed up over all relevant processes by this
     * precoditioner before the inverse is computed.
     */
    template<class Diagonal, class X, class Y>
    class NonoverlappingJacobi : public Dune::Preconditioner<X,Y>
    {
      typedef typename Diagonal::Backend DBackend;

      std::size_t gfsSize;
      Diagonal diagonal;

    public:
      //! The domain type of the operator.
      /**
       * The preconditioner is an inverse operator, so this is the output type
       * of the preconditioner.
       */
      typedef X domain_type;
      //! \brief The range type of the operator.
      /**
       * The preconditioner is an inverse operator, so this is the input type
       * of the preconditioner.
       */
      typedef Y range_type;
      //! \brief The field type of the preconditioner.
      typedef typename X::ElementType field_type;

      enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::nonoverlapping
      };

      //! \brief Constructor.
      /**
       * \param gfs The GridFunctionSpace the matrix and the vectors live on.
       * \param m   The matrix whose inverse the preconditioner should
       *            estimate.  m is assumed to be inconsistent (i.e. rows for
       *            dofs on the border only contain the contribution of the
       *            local process).
       *
       * The preconditioner does not store any reference to the gfs or the
       * matrix m.  The diagonal of m is copied, since it has to be made
       * consistent.
       */
      template<class GFS, class Matrix>
      NonoverlappingJacobi(const GFS& gfs, const Matrix &m) :
        gfsSize(gfs.size()), diagonal(gfs)
      {
        typedef typename Matrix::Backend MBackend;

        for(std::size_t i = 0; i < gfsSize; ++i)
          DBackend::access(diagonal, i) = MBackend::access(m, i, i);

        Dune::PDELab::AddDataHandle<GFS, Diagonal> addDH(gfs, diagonal);
        gfs.gridView().communicate(addDH,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);
      }

      //! Prepare the preconditioner.
      virtual void pre (X& x, Y& b) {}

      //! Apply the precondioner.
      /*
       * For this preconditioner, this method works with both consistent and
       * inconsistent vectors: if d is consistent, v will be consistent, if d
       * is inconsistent, v will be inconsistent.
       */
      virtual void apply (X& v, const Y& d)
      {
        typedef typename X::Backend XBackend;
        typedef typename Y::Backend YBackend;

        for(std::size_t i = 0; i < gfsSize; ++i)
          XBackend::access(v, i) =
            YBackend::access(d, i) / DBackend::access(diagonal, i);
      }

      //! Clean up.
      virtual void post (X& x) {}
    };

    //! \addtogroup PDELab_novlpsolvers Nonoverlapping Solvers
    //! \{

    //! \brief Nonoverlapping parallel CG solver without preconditioner
    template<class GFS>
    class ISTLBackend_NOVLP_CG_NOPREC
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_NOVLP_CG_NOPREC (const GFS& gfs_,
                                            unsigned maxiter_=5000,
                                            int verbose_=1)
        : gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        typedef NonoverlappingRichardson<GFS,V,W> PRICH;
        PRICH prich(gfs,phelper);
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Dune::CGSolver<V> solver(pop,psp,prich,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int verbose;
    };

    //! \brief Nonoverlapping parallel CG solver with Jacobi preconditioner
    template<class GFS>
    class ISTLBackend_NOVLP_CG_Jacobi
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int verbose;

    public:
      //! make a linear solver object
      /**
       * \param gfs_     A grid function space
       * \param maxiter_ Maximum number of iterations to do.
       * \param verbose_ Verbosity level, directly handed to the CGSolver.
       */
      explicit ISTLBackend_NOVLP_CG_Jacobi(const GFS& gfs_,
                                           unsigned maxiter_ = 5000,
                                           int verbose_ = 1) :
        gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
      {}

      //! compute global norm of a vector
      /**
       * \param v The vector to compute the norm of.  Should be an
       *          inconsistent vector (i.e. the entries corresponding a DoF on
       *          the border should only contain the summand of this process).
       */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      //! solve the given linear system
      /**
       * \param A         The matrix to solve.  Should be a matrix from one of
       *                  PDELabs ISTL backends (only ISTLBCRSMatrixBackend at
       *                  the moment).
       * \param z         The solution vector to be computed
       * \param r         Right hand side
       * \param reduction to be achieved
       *
       * Solve the linear system A*z=r such that
       * norm(A*z0-r)/norm(A*z-r) < reduction where z0 is the initial value of
       * z.
       */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);

        typedef typename M::ElementType MField;
        typedef typename Dune::PDELab::BackendVectorSelector<GFS,MField>::Type Diagonal;
        typedef NonoverlappingJacobi<Diagonal,V,W> PPre;
        PPre ppre(gfs,A);

        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Dune::CGSolver<V> solver(pop,psp,ppre,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      //! Return access to result data
      const Dune::PDELab::LinearSolverResult<double>& result() const
      { return res; }
    };

    //! \brief Nonoverlapping parallel BiCGStab solver without preconditioner
    template<class GFS>
    class ISTLBackend_NOVLP_BCGS_NOPREC
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_NOVLP_BCGS_NOPREC (const GFS& gfs_, unsigned maxiter_=5000, int verbose_=1)
        : gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_), verbose(verbose_)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      {
        typedef NonoverlappingOperator<GFS,M,V,W> POP;
        POP pop(gfs,A);
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        typedef NonoverlappingRichardson<GFS,V,W> PRICH;
        PRICH prich(gfs,phelper);
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Dune::BiCGSTABSolver<V> solver(pop,psp,prich,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      int verbose;
    };

    //! Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
    template<typename GFS>
    class ISTLBackend_NOVLP_ExplicitDiagonal
    {
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;

      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;

    public:
      /*! \brief make a linear solver object

        \param[in] gfs_ GridFunctionSpace, used to identify DoFs for parallel
        communication
      */
      explicit ISTLBackend_NOVLP_ExplicitDiagonal(const GFS& gfs_)
        : gfs(gfs_), phelper(gfs)
      {}

      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      template<class V>
      typename V::ElementType norm (const V& v) const
      {
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        V x(v); // make a copy because it has to be made consistent
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief solve the given linear system

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename W::ElementType reduction)
      {
        Dune::SeqJac<M,V,W> jac(A,1,1.0);
        jac.pre(z,r);
        jac.apply(z,r);
        jac.post(z);
        if (gfs.gridView().comm().size()>1)
        {
          Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,z);
          gfs.gridView().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
        }
        res.converged  = true;
        res.iterations = 1;
        res.elapsed    = 0.0;
        res.reduction  = reduction;
        res.conv_rate  = reduction; // pow(reduction,1.0/1)
      }

      /*! \brief Return access to result data */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }
    };
    //! \} Nonoverlapping Solvers

    /**
    * @brief Helper class for adding up matrix entries on border.
    * @tparam GridOperator The grid operator to work on.
    * @tparam MatrixType The MatrixType.
    */
    template<class GridOperator, class MatrixType>
    class VertexExchanger
    {
      typedef MatrixType Matrix;
      typedef typename GridOperator::Traits GridOperatorTraits;
      typedef typename GridOperatorTraits::JacobianField Scalar;
      typedef typename GridOperatorTraits::TrialGridFunctionSpace GFS;
      typedef typename GFS::Traits::GridViewType GridView;
      enum {dim = GridView::dimension};
      typedef typename GridView::Traits::Grid Grid;
      typedef typename Matrix::block_type BlockType;
      typedef typename GridView::template Codim<dim>::Iterator  VertexIterator;
      typedef typename Grid::Traits::GlobalIdSet IDS;
      typedef typename IDS::IdType IdType;
      typedef typename Matrix::RowIterator RowIterator;
      typedef typename Matrix::ColIterator ColIterator;

    public:
      /*! \brief Constructor. Sets up the local to global relations.

      \param[in] gridView The grid view to operate on.
      */
      VertexExchanger(const GridView& gridView)
        : gridView_(gridView)
      {
        gid2Index_.clear();
        index2GID_.clear();


        VertexIterator vertexEndIt = gridView_.template end<dim>();
        for (VertexIterator vertexIt = gridView_.template begin<dim>(); vertexIt != vertexEndIt; ++vertexIt)
        {
          if (vertexIt->partitionType() == Dune::BorderEntity)
          {
            int localIdx = gridView_.indexSet().index(*vertexIt);
            IdType globalIdx = gridView_.grid().globalIdSet().id(*vertexIt);

            std::pair<IdType,int> g2iPair(globalIdx, localIdx);
            gid2Index_.insert(g2iPair);

            std::pair<int,IdType> i2gPair(localIdx, globalIdx);
            index2GID_.insert(i2gPair);

          }
        }
      }

      //! A DataHandle class to exchange matrix sparsity patterns
      class MatPatternExchange
        : public Dune::CommDataHandleIF<MatPatternExchange,IdType> {
        typedef typename Matrix::RowIterator RowIterator;
        typedef typename Matrix::ColIterator ColIterator;
      public:
        //! Export type of data for message buffer
        typedef IdType DataType;

        /** @brief Returns true if data for given valid codim should be communicated   
         */
        bool contains (int dim, int codim) const
        {
          return (codim==dim);
        }

        /** @brief Returns true if size of data per entity of given dim and codim is a constant
        */
        bool fixedsize (int dim, int codim) const
        {
          return false;
        }

        /** @brief How many objects of type DataType have to be sent for a given entity
        */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
          int i = gridView_.indexSet().index(e);
          int n = 0;
          for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
              typename std::map<int,IdType>::const_iterator it = index2GID_.find(j.index());
              if (it != index2GID_.end())
                n++;
            }
          
          return n;
        }

        /** @brief Pack data from user to message buffer
        */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
          int i = gridView_.indexSet().index(e);
          for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
              typename std::map<int,IdType>::const_iterator it=index2GID_.find(j.index());
              if (it != index2GID_.end())
                buff.write(it->second);
            }
        
        }

        /** @brief Unpack data from message buffer to user
        */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
          int i = gridView_.indexSet().index(e);
          for (size_t k = 0; k < n; k++)
            {
              IdType id;
              buff.read(id);
              // only add entries corresponding to border entities
              typename std::map<IdType,int>::const_iterator it = gid2Index_.find(id);
              if (it != gid2Index_.end() 
                  && sparsity_[i].find(it->second) == sparsity_[i].end()
                  && helper_.ghost(it->second,0) != 1<<24)
                sparsity_[i].insert(it->second);
            }
        }

        /**
         * @brief Get the communicated sparsity pattern
         * @return the vector with the sparsity pattern
         */
        std::vector<std::set<int> >& sparsity ()
        {
          return sparsity_;
        }

        /** @brief Constructor
            @param[in] gridView Grid view.
            @param[in] g2i Global to local index map.
            @param[in] i2g Local to global index map.
            @param[in] A Matrix to operate on.
            @param[in] helper parallel istl helper.
        */
        MatPatternExchange (const GridView& gridView,
                            const std::map<IdType,int>& g2i,
                            const std::map<int,IdType>& i2g, Matrix& A, 
                            const Dune::PDELab::ParallelISTLHelper<GFS>& helper)
          : gridView_(gridView), gid2Index_(g2i), index2GID_(i2g),
            sparsity_(A.N()), A_(A), helper_(helper)
        {}

      private:
        const GridView& gridView_;
        const std::map<IdType,int>& gid2Index_;
        const std::map<int,IdType>& index2GID_;
        std::vector<std::set<int> > sparsity_;
        Matrix& A_;
        const Dune::PDELab::ParallelISTLHelper<GFS>& helper_;
      };
      
      //! Local matrix blocks associated with the global id set
      struct MatEntry
      {
        IdType first;
        BlockType second;
        MatEntry (const IdType& f, const BlockType& s) : first(f),second(s) {}
        MatEntry () {}
      };

      //! A DataHandle class to exchange matrix entries
      class MatEntryExchange
      : public Dune::CommDataHandleIF<MatEntryExchange,MatEntry> {
        typedef typename Matrix::RowIterator RowIterator;
        typedef typename Matrix::ColIterator ColIterator;
      public:
        //! Export type of data for message buffer
        typedef MatEntry DataType;

        /** @brief Returns true if data for given valid codim should be communicated
        */
        bool contains (int dim, int codim) const
        {
          return (codim==dim);
        }

        /** @brief Returns true if size of data per entity of given dim and codim is a constant
        */
        bool fixedsize (int dim, int codim) const
        {
          return false;
        }

        /** @brief How many objects of type DataType have to be sent for a given entity
        */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
          int i = gridView_.indexSet().index(e);
          int n = 0;
          for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
            typename std::map<int,IdType>::const_iterator it = index2GID_.find(j.index());
            if (it != index2GID_.end())
              n++;
            }
          
            return n;
        }

        /** @brief Pack data from user to message buffer
        */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
          int i = gridView_.indexSet().index(e);
          for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
              typename std::map<int,IdType>::const_iterator it=index2GID_.find(j.index());
              if (it != index2GID_.end())
                buff.write(MatEntry(it->second,*j));
            }
        
        }

        /** @brief Unpack data from message buffer to user
        */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
          int i = gridView_.indexSet().index(e);
          for (size_t k = 0; k < n; k++)
          {
            MatEntry m;
            buff.read(m);
            // only add entries corresponding to border entities
            typename std::map<IdType,int>::const_iterator it = gid2Index_.find(m.first);
            if (it != gid2Index_.end())
              if (A_[i].find(it->second) != A_[i].end())
                A_[i][it->second] += m.second;
          }
        }

        /** @brief Constructor
        @param[in] gridView Grid view.
        @param[in] g2i Global to local index map.
        @param[in] i2g Local to global index map.
        @param[in] A Matrix to operate on.
        */
        MatEntryExchange (const GridView& gridView, const std::map<IdType,int>& g2i,
                          const std::map<int,IdType>& i2g,
                          Matrix& A)
          : gridView_(gridView), gid2Index_(g2i), index2GID_(i2g), A_(A)
        {}

      private:
        const GridView& gridView_;
        const std::map<IdType,int>& gid2Index_;
        const std::map<int,IdType>& index2GID_;
        Matrix& A_;
      };

      /** @brief communicates values for the sparsity pattern of the new matrix.
          @param A Matrix to operate on.
          @param helper ParallelelISTLHelper.
      */
      void getextendedmatrix (Matrix& A,const Dune::PDELab::ParallelISTLHelper<GFS>& helper)
      {
        if (gridView_.comm().size() > 1) {
          Matrix tmp(A);
          std::size_t nnz=0;        
          // get entries from other processes
          MatPatternExchange datahandle(gridView_, gid2Index_, index2GID_, A, helper);
          gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
          std::vector<std::set<int> >& sparsity = datahandle.sparsity();
          // add own entries, count number of nonzeros
          for (RowIterator i = A.begin(); i != A.end(); ++i){
            for (ColIterator j = A[i.index()].begin(); j != A[i.index()].end(); ++j){
              if (sparsity[i.index()].find(j.index()) == sparsity[i.index()].end())
                sparsity[i.index()].insert(j.index());
            }
            nnz += sparsity[i.index()].size();
          }
          A.setSize(tmp.N(), tmp.N(), nnz);
          A.setBuildMode(Matrix::row_wise);
          typename Matrix::CreateIterator citer = A.createbegin();
          typedef typename std::vector<std::set<int> >::const_iterator Iter;
          for (Iter i = sparsity.begin(), end = sparsity.end(); i!=end; ++i, ++citer){
            typedef typename std::set<int>::const_iterator SIter;
            for (SIter si = i->begin(), send = i->end(); si!=send; ++si)
              citer.insert(*si);
          }
          // set matrix old values
          A = 0;
          for (RowIterator i = tmp.begin(); i != tmp.end(); ++i)
            for (ColIterator j = tmp[i.index()].begin(); j != tmp[i.index()].end(); ++j){
              A[i.index()][j.index()] = tmp[i.index()][j.index()];
            }
        }
      }

      /** @brief Sums up the entries corresponding to border vertices.
      @param A Matrix to operate on.
      */
      void sumEntries (Matrix& A)
      {
        if (gridView_.comm().size() > 1)
        {
          MatEntryExchange datahandle(gridView_, gid2Index_, index2GID_, A);
          gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
        }
      }

      private:
      const GridView& gridView_;
      std::map<IdType,int> gid2Index_;
      std::map<int,IdType> index2GID_;
    };
      
    template<class GO, 
             template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_NOVLP_BASE_PREC
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      typedef Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;
      
    public:
      /*! \brief Constructor.

        \param[in] gfs_ a grid function space
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] steps_ number of preconditioner steps to apply as inner iteration
        \param[in] verbose_ print messages if true
      */
      explicit ISTLBackend_NOVLP_BASE_PREC (const GFS& gfs_, unsigned maxiter_ = 5000, unsigned steps_ = 5, int verbose_ = 1)
        : gfs(gfs_), phelper(gfs_,verbose_),
          maxiter(maxiter_), steps(steps_), verbose(verbose_)
      {}

      /*! \brief Compute global norm of a vector.

        \param[in] v the given vector
      */
      template<class Vector>
      typename Vector::ElementType norm (const Vector& v) const
      {
        Vector x(v); // make a copy because it has to be made consistent
        typedef NonoverlappingScalarProduct<GFS,Vector> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }

      /*! \brief Solve the given linear system.

        \param[in] A the given matrix
        \param[out] z the solution vector to be computed
        \param[in] r right hand side
        \param[in] reduction to be achieved
      */
      template<class M, class V, class W>
      void apply(M& A, V& z, W& r, typename V::ElementType reduction)
      { 
        typedef typename Dune::PDELab::CommSelector<96,Dune::MPIHelper::isFake>::type Comm;
        typedef typename M::BaseT MatrixType;
        MatrixType& mat=A.base();
        typedef typename Dune::PDELab::BlockProcessor<GFS>::template AMGVectorTypeSelector<V>::Type VectorType;
#if HAVE_MPI
        Comm oocc(gfs.gridView().comm(),Dune::SolverCategory::nonoverlapping);
        typedef VertexExchanger<GO,MatrixType> Exchanger;
        Exchanger exchanger(gfs.gridView());
        exchanger.getextendedmatrix(mat,phelper);
        exchanger.sumEntries(mat);
        phelper.createIndexSetAndProjectForAMG(mat, oocc);
        typedef Preconditioner<MatrixType,VectorType,VectorType,1> Smoother;
        Smoother smoother(mat, steps, 1.0);
        typedef Dune::NonoverlappingSchwarzScalarProduct<VectorType,Comm> PSP;
        PSP psp(oocc);      
        typedef Dune::NonoverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> Operator;
        Operator oop(mat,oocc);
        typedef Dune::NonoverlappingBlockPreconditioner<Comm, Smoother> ParSmoother;
        ParSmoother parsmoother(smoother, oocc);
#else
        Comm oocc(gfs.gridView().comm());
        typedef Preconditioner<MatrixType,VectorType,VectorType,1> ParSmoother;
        ParSmoother parsmoother(mat, steps, 1.0);
        typedef Dune::SeqScalarProduct<VectorType> PSP;
        PSP psp;
        typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
        Operator oop(mat);
#endif
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        Solver<VectorType> solver(oop,psp,parsmoother,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        //make r consistent
        if (gfs.gridView().comm().size()>1){
          Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,r);
          gfs.gridView().communicate(adddh,
                                     Dune::InteriorBorder_InteriorBorder_Interface,
                                     Dune::ForwardCommunication);
        }
        
        solver.apply(z,r,stat);
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }

      /*! \brief Return access to result data. */
      const Dune::PDELab::LinearSolverResult<double>& result() const
      {
        return res;
      }

    private:
      const GFS& gfs;
      PHELPER phelper;
      Dune::PDELab::LinearSolverResult<double> res;
      unsigned maxiter;
      unsigned steps;
      int verbose;
    };

    //! \addtogroup PDELab_novlpsolvers Nonoverlapping Solvers
    //! \{
    
  /**
   * @brief Nonoverlapping parallel BiCGSTAB solver preconditioned by block SSOR.
   * @tparam GO The type of the grid operator 
   * (or the fakeGOTraits class for the old grid operator space).
   *
   * The solver uses a NonoverlappingBlockPreconditioner with underlying
   * sequential SSOR preconditioner. The crucial step is to add up the matrix entries
   * corresponding to the border vertices on each process. This is achieved by
   * performing a VertexExchanger::sumEntries(Matrix&) before constructing the
   * sequential SSOR.
   */
  template<class GO>
  class ISTLBackend_NOVLP_BCGS_SSORk
    : public ISTLBackend_NOVLP_BASE_PREC<GO,Dune::SeqSSOR, Dune::BiCGSTABSolver>
  {
    typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      
  public:
    /*! \brief make a linear solver object
        
      \param[in] gfs_ a grid function space
      \param[in] maxiter_ maximum number of iterations to do
      \param[in] steps_ number of SSOR steps to apply as inner iteration
      \param[in] verbose_ print messages if true
    */
    explicit ISTLBackend_NOVLP_BCGS_SSORk (const GFS& gfs_, unsigned maxiter_=5000,
                                           int steps_=5, int verbose_=1)
      : ISTLBackend_NOVLP_BASE_PREC<GO,Dune::SeqSSOR, Dune::BiCGSTABSolver>(gfs_, maxiter_, steps_, verbose_)
    {}
  };

  /**
   * @brief Nonoverlapping parallel CG solver preconditioned by block SSOR.
   */
  template<class GO>
  class ISTLBackend_NOVLP_CG_SSORk
    : public ISTLBackend_NOVLP_BASE_PREC<GO,Dune::SeqSSOR, Dune::CGSolver>
  {
    typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      
  public:
    /*! \brief make a linear solver object
        
      \param[in] gfs_ a grid function space
      \param[in] maxiter_ maximum number of iterations to do
      \param[in] steps_ number of SSOR steps to apply as inner iteration
      \param[in] verbose_ print messages if true
    */
    explicit ISTLBackend_NOVLP_CG_SSORk (const GFS& gfs_, unsigned maxiter_=5000,
                                         int steps_=5, int verbose_=1)
      : ISTLBackend_NOVLP_BASE_PREC<GO,Dune::SeqSSOR, Dune::CGSolver>(gfs_, maxiter_, steps_, verbose_)
    {}
  };
    //! \} Nonoverlapping Solvers
    //! \} group Backend
    
    template<class GO,int s, template<class,class,class,int> class Preconditioner,
             template<class> class Solver>
    class ISTLBackend_AMG_NOVLP : public Dune::PDELab::LinearResultStorage
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      typedef typename Dune::PDELab::ParallelISTLHelper<GFS> PHELPER;
      typedef typename GO::Traits::Jacobian M;
      typedef typename M::BaseT MatrixType;
      typedef typename GO::Traits::Domain V;
      typedef typename Dune::PDELab::BlockProcessor<GFS>::template AMGVectorTypeSelector<V>::Type VectorType;
      typedef typename Dune::PDELab::CommSelector<s,Dune::MPIHelper::isFake>::type Comm;
#if HAVE_MPI
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> Smoother;
      typedef Dune::NonoverlappingBlockPreconditioner<Comm,Smoother> ParSmoother;
      typedef Dune::NonoverlappingSchwarzOperator<MatrixType,VectorType,VectorType,Comm> Operator;
#else
      typedef Preconditioner<MatrixType,VectorType,VectorType,1> ParSmoother;
      typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#endif
      typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments SmootherArgs;
      typedef Dune::Amg::AMG<Operator,VectorType,ParSmoother,Comm> AMG;
      typedef Dune::Amg::Parameters Parameters;

    public:
      ISTLBackend_AMG_NOVLP(const GFS& gfs_, unsigned maxiter_=5000, 
                            int verbose_=1, bool reuse_=false,
                            bool usesuperlu_=true)
        : gfs(gfs_), phelper(gfs,verbose_), maxiter(maxiter_),
          params(15,2000,1.2,1.6,Dune::Amg::noAccu/*Dune::Amg::atOnceAccu*/),
          verbose(verbose_), reuse(reuse_), firstapply(true),
          usesuperlu(usesuperlu_)
      {
        params.setDefaultValuesIsotropic(GFS::Traits::GridViewType::Traits::Grid::dimension);
        params.setDebugLevel(verbose_);
#if !HAVE_SUPERLU
        if (gfs.gridView().comm().rank() == 0 && usesuperlu == true)
          {
            std::cout << "WARNING: You are using AMG without SuperLU!"
                      << " Please consider installing SuperLU," 
                      << " or set the usesuperlu flag to false"
                      << " to suppress this warning." << std::endl;
          }
#endif
      }

       /*! \brief set AMG parameters

        \param[in] params_ a parameter object of Type Dune::Amg::Parameters
      */     
      void setParameters(const Parameters& params_)
      {
        params = params_;
      }
        
      void setparams(Parameters params_)
      {
        params = params_;
      }

      /** 
       * @brief Get the parameters describing the behaviuour of AMG.
       *
       * The returned object can be adjusted to ones needs and then can be
       * reset using setParameters.
       * @return The object holding the parameters of AMG.
       */
      const Parameters& parameters() const
      {
        return params;
      }
      
      /*! \brief compute global norm of a vector

        \param[in] v the given vector
      */
      typename V::ElementType norm (const V& v) const
      {
        V x(v); // make a copy because it has to be made consistent
        typedef NonoverlappingScalarProduct<GFS,V> PSP;
        PSP psp(gfs,phelper);
        psp.make_consistent(x);
        return psp.norm(x);
      }
      
      void apply(MatrixType& A, V& z, V& r, typename V::ElementType reduction)
      {
        Dune::Timer watch;
        MatrixType& mat=A;
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
          Dune::Amg::FirstDiagonal> > Criterion;
#if HAVE_MPI
        Comm oocc(gfs.gridView().comm(),Dune::SolverCategory::nonoverlapping);
        typedef VertexExchanger<GO,MatrixType> Exchanger;
        Exchanger exchanger(gfs.gridView());
        exchanger.getextendedmatrix(mat,phelper);
        exchanger.sumEntries(mat);        
        phelper.createIndexSetAndProjectForAMG(mat, oocc);
        Dune::NonoverlappingSchwarzScalarProduct<VectorType,Comm> sp(oocc);
        Operator oop(mat, oocc);
#else
        Comm oocc(gfs.gridView().comm());
        Operator oop(mat);
        Dune::SeqScalarProduct<VectorType> sp;
#endif
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;
        //use noAccu or atOnceAccu
        Criterion criterion(params);
        stats.tprepare=watch.elapsed();
        watch.reset();
        
        int verb=0;
        if (gfs.gridView().comm().rank()==0) verb=verbose;
        //only construct a new AMG if the matrix changes
        if (reuse==false || firstapply==true){
          amg.reset(new AMG(oop, criterion, smootherArgs, oocc));
          firstapply = false;
          stats.tsetup = watch.elapsed();
          stats.levels = amg->maxlevels();
          stats.directCoarseLevelSolver=amg->usesDirectCoarseLevelSolver();
        }
        
        Dune::InverseOperatorResult stat;
        // make r consistent
        if (gfs.gridView().comm().size()>1) {
          Dune::PDELab::AddDataHandle<GFS,V> adddh(gfs,r);
          gfs.gridView().communicate(adddh,
                                     Dune::InteriorBorder_InteriorBorder_Interface,
                                     Dune::ForwardCommunication);
        }
        watch.reset();
        Solver<VectorType> solver(oop,sp,*amg,reduction,maxiter,verb);
        solver.apply(Dune::PDELab::BlockProcessor<GFS>::getVector(z),Dune::PDELab::BlockProcessor<GFS>::getVector(r),stat);
        stats.tsolve= watch.elapsed();
        res.converged  = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed    = stat.elapsed;
        res.reduction  = stat.reduction;
        res.conv_rate  = stat.conv_rate;
      }
      
      /** 
       * @brief Get statistics of the AMG solver (no of levels, timings). 
       * @return statistis of the AMG solver. 
       */
      const ISTLAMGStatistics& statistics() const
      {
        return stats;
      }
      
    private:
      const GFS& gfs;
      PHELPER phelper;
      unsigned maxiter;
      Parameters params;
      int verbose;
      bool reuse;
      bool firstapply;
      bool usesuperlu;
      Dune::shared_ptr<AMG> amg;
      ISTLAMGStatistics stats;
    };

    template<class GO, int s=96>
    class ISTLBackend_NOVLP_CG_AMG_SSOR
      : public ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::CGSolver>
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      
    public:
      ISTLBackend_NOVLP_CG_AMG_SSOR(const GFS& gfs_, unsigned maxiter_=5000, 
                                    int verbose_=1, bool reuse_=false,
                                    bool usesuperlu_=true)
        : ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::CGSolver>
          (gfs_, maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };


    template<class GO, int s=96>
    class ISTLBackend_NOVLP_BCGS_AMG_SSOR
      : public ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::BiCGSTABSolver>
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      
    public:
      ISTLBackend_NOVLP_BCGS_AMG_SSOR(const GFS& gfs_, unsigned maxiter_=5000,
                                      int verbose_=1, bool reuse_=false,
                                      bool usesuperlu_=true)
        : ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::BiCGSTABSolver>
          (gfs_, maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

    template<class GO, int s=96>
    class ISTLBackend_NOVLP_LS_AMG_SSOR
      : public ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::LoopSolver>
    {
      typedef typename GO::Traits::TrialGridFunctionSpace GFS;
      
    public:
      ISTLBackend_NOVLP_LS_AMG_SSOR(const GFS& gfs_, unsigned maxiter_=5000, 
                                    int verbose_=1, bool reuse_=false,
                                    bool usesuperlu_=true)
        : ISTLBackend_AMG_NOVLP<GO, s, Dune::SeqSSOR, Dune::LoopSolver>
          (gfs_, maxiter_, verbose_, reuse_, usesuperlu_)
      {}
    };

} // namespace Dumux

#endif
