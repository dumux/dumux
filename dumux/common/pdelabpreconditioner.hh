// $Id$: preconditionerpdelab.hh 3728 2010-06-10 15:44:39Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
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
 * \brief A parallel solver for nonoverlapping grid partitions.
 */
#ifndef DUMUX_PDELAB_PRECONDITIONER_HH
#define DUMUX_PDELAB_PRECONDITIONER_HH

#include <dumux/common/pardiso.hh>
#include <dumux/common/propertysystem.hh>

#include <dune/pdelab/backend/istlsolverbackend.hh>

namespace Dumux {
// forward declaration of property tags
namespace Properties {
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Model);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(NumEq);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(VertexMapper);
NEW_PROP_TAG(GridOperatorSpace);
NEW_PROP_TAG(JacobianMatrix);
NEW_PROP_TAG(ReferenceElements);
NEW_PROP_TAG(GridFunctionSpace);
NEW_PROP_TAG(ConstraintsTrafo);
};


namespace PDELab {

/*!
 * \brief exchanges matrix entries for parallel computations
 */
template<class TypeTag>
class Exchanger
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        dim = GridView::dimension
    };
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) JacobianMatrix;

    typedef typename JacobianMatrix::block_type BlockType;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename Grid::Traits::GlobalIdSet IDS;
    typedef typename IDS::IdType IdType;
    typedef typename GridView::ctype CoordScalar;

    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

public:

    Exchanger(const Problem& problem)
        : gridView_(problem.gridView()), vertexMapper_(problem.vertexMapper()), borderIndices_(0)
    {
        gid2Index_.clear();
        index2GID_.clear();

        VertexIterator vertexEndIt = gridView_.template end<dim>();
        for (VertexIterator vertexIt = gridView_.template begin<dim>(); vertexIt != vertexEndIt; ++vertexIt)
        {
//            std::cout << gridView_.comm().rank() << ": node " << vertexMapper_.map(*vertexIt)
//                    << " at (" << vertexIt->geometry().corner(0) << ") is of type "
//                    << vertexIt->partitionType() << ", GID = "
//                    << gridView_.grid().globalIdSet().id(*vertexIt) << std::endl;
            if (vertexIt->partitionType() == Dune::BorderEntity)
            {
                int localIdx = vertexMapper_.map(*vertexIt);
                IdType globalIdx = gridView_.grid().globalIdSet().id(*vertexIt);

                std::pair<IdType,int> g2iPair(globalIdx, localIdx);
                gid2Index_.insert(g2iPair);

                std::pair<int,IdType> i2gPair(localIdx, globalIdx);
                index2GID_.insert(i2gPair);

                borderIndices_.push_back(localIdx);
            }
        }
    }

    /*!
     * \brief matrix entry for the MatEntryExchange class
     */
    struct MatEntry
    {
        IdType first;
        BlockType second;
        MatEntry (const IdType& f, const BlockType& s) : first(f),second(s) {}
        MatEntry () {}
    };

    /*!
     * \brief A DataHandle class to exchange matrix entries
     */
    class MatEntryExchange
        : public Dune::CommDataHandleIF<MatEntryExchange,MatEntry>
    {
        typedef typename JacobianMatrix::RowIterator RowIterator;
        typedef typename JacobianMatrix::ColIterator ColIterator;
    public:
        //! export type of data for message buffer
        typedef MatEntry DataType;

    //! returns true if data for this codim should be communicated
    bool contains (int dim, int codim) const
    {
        return (codim==dim);
    }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedsize (int dim, int codim) const
    {
        return false;
    }

    /*! how many objects of type DataType have to be sent for a given entity

          Note: Only the sender side needs to know this size.
     */
    template<class EntityType>
    size_t size (EntityType& e) const
    {
//        std::cout << gridView_.comm().rank() << ": begin loop over vertices.\n";
//        VertexIterator vertexEndIt = gridView_.template end<dim>();
//        for (VertexIterator vertexIt = gridView_.template begin<dim>(); vertexIt != vertexEndIt; ++vertexIt)
//        {
//            std::cout << gridView_.comm().rank() << ": node " << vertexMapper_.map(*vertexIt)
//                    << " at (" << vertexIt->geometry().corner(0) << ") is of type "
//                    << vertexIt->partitionType() << ", GID = "
//                    << gridView_.grid().globalIdSet().id(*vertexIt) << std::endl;
//        }
//        std::cout << gridView_.comm().rank() << ": end loop over vertices.\n";
//        std::cout.flush();
//        std::cout << gridView_.comm().rank() << ": node " << vertexMapper_.map(e)
//                << " on level " << e.level() << " at (" << e.geometry().corner(0) << ") is of type "
//                << e.partitionType() << ", GID = "
//                << gridView_.grid().globalIdSet().id(e) << std::endl;;
        int i = vertexMapper_.map(e);
        int n = 0;
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                // only count those entries corresponding to border entities
                typename std::map<int,IdType>::const_iterator it = index2GID_.find(j.index());
                if (it != index2GID_.end())
                    n++;
            }
//        std::cout << gridView_.comm().rank() << ": node " << i << " has sending size " << n << std::endl;
        return n;
    }

    //! pack data from user to message buffer
    template<class MessageBuffer, class EntityType>
    void gather (MessageBuffer& buff, const EntityType& e) const
    {
            int i = vertexMapper_.map(e);
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                // only send those entries corresponding to border entities
                typename std::map<int,IdType>::const_iterator it=index2GID_.find(j.index());
                if (it != index2GID_.end()) {
                    buff.write(MatEntry(it->second,*j));
//                    std::cout << gridView_.comm().rank() << ": node " << i << " gathers (" << it->second
//                            << ", " << *j << ") for j = " << j.index() << std::endl;
                }
            }
    }

    /*! unpack data from message buffer to user

          n is the number of objects sent by the sender
     */
    template<class MessageBuffer, class EntityType>
    void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
    {
           int i = vertexMapper_.map(e);
            for (size_t k = 0; k < n; k++)
            {
                MatEntry m;
                buff.read(m);
                // only add entries corresponding to border entities
                typename std::map<IdType,int>::const_iterator it = gid2Index_.find(m.first);
                if (it != gid2Index_.end())
                {
//                    std::cout << gridView_.comm().rank() << ": node " << i << " adds " << m.second
//                            << " to j = " << it->second << ", GID = " << m.first << std::endl;
                    if (A_[i].find(it->second) != A_[i].end())
                        A_[i][it->second] += m.second;
                }
            }
    }

    //! constructor
    MatEntryExchange (const GridView& gridView, const std::map<IdType,int>& g2i,
            const std::map<int,IdType>& i2g,
            const VertexMapper& vm,
            JacobianMatrix& A)
            : gridView_(gridView), gid2Index_(g2i), index2GID_(i2g), vertexMapper_(vm), A_(A)
            {}

private:
    const GridView& gridView_;
    const std::map<IdType,int>& gid2Index_;
    const std::map<int,IdType>& index2GID_;
    const VertexMapper& vertexMapper_;
    JacobianMatrix& A_;
    };

    void sumEntries (JacobianMatrix& A)
    {
      if (gridView_.comm().size() > 1)
      {
          MatEntryExchange datahandle(gridView_, gid2Index_, index2GID_, vertexMapper_, A);
          gridView_.communicate(datahandle,
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);
      }
    }

    const std::vector<int>& borderIndices()
    {
       return borderIndices_;
    }

private:
    const GridView& gridView_;
    std::map<IdType,int> gid2Index_;
    std::map<int,IdType> index2GID_;
    const VertexMapper& vertexMapper_;
    std::vector<int> borderIndices_;
};

/*!
 * \brief Wrapper for a sequential preconditioner
 */
template<class CC, class GFS, class P>
class NonoverlappingWrappedPreconditioner
  : public Dune::Preconditioner<typename P::domain_type,typename P::range_type>
{
public:
  //! \brief The domain type of the preconditioner.
  typedef typename P::domain_type domain_type;
  //! \brief The range type of the preconditioner.
  typedef typename P::range_type range_type;

  // define the category
  enum {
    //! \brief The category the preconditioner is part of.
    category=Dune::SolverCategory::nonoverlapping
  };

      //! Constructor.
  NonoverlappingWrappedPreconditioner (const GFS& gfs_, P& prec_, const CC& cc_,
                                       const std::vector<int>& borderIndices,
                                       const Dune::PDELab::ParallelISTLHelper<GFS>& helper_)
    : gfs(gfs_), prec(prec_), cc(cc_), borderIndices_(borderIndices), helper(helper_)
  {}

  /*!
    \brief \copybrief Dune::SeqPardiso::pre

    \copydetails Dune::SeqPardiso::pre
  */
  virtual void pre (domain_type& x, range_type& b)
  {
    prec.pre(x,b);
  }

  /*!
    \brief \copybrief Dune::SeqPardiso::apply(X&,const Y&)

    \copydetails Dune::SeqPardiso::apply(X&,const Y&)
  */
  virtual void apply (domain_type& v, const range_type& d)
  {
    range_type dd(d);
    set_constrained_dofs(cc,0.0,dd);
    prec.apply(v,dd);

    Dune::PDELab::AddDataHandle<GFS,domain_type> adddh(gfs,v);
    if (gfs.gridview().comm().size()>1)
      gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

    for (int k = 0; k < borderIndices_.size(); k++)
        v[borderIndices_[k]] *= 0.5;
  }

  /*!
    \brief \copybrief Dune::SeqPardiso::post(X&)

    \copydetails Dune::SeqPardiso::post(X&)
  */
  virtual void post (domain_type& x)
  {
    prec.post(x);
  }

private:
  const GFS& gfs;
  P& prec;
  const CC& cc;
  const std::vector<int>& borderIndices_;
    const Dune::PDELab::ParallelISTLHelper<GFS>& helper;
};

/*!
 * \brief backend for an ISTL parallel ILU preconditioned BiCGSTAB solver
 */
template<class TypeTag>
class ISTLBackend_NoOverlap_BCGS_ILU
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConstraintsTrafo)) ConstraintsTrafo;
    typedef Dune::PDELab::ParallelISTLHelper<GridFunctionSpace> PHELPER;

public:
    /*! \brief make a linear solver object

    \param[in] problem The Dumux problem
    \param[in] maxiter_ Maximum number of iterations to do
    \param[in] verbose_ Verbosity level
    */
    explicit ISTLBackend_NoOverlap_BCGS_ILU (Problem& problem, unsigned maxiter_=5000, int verbose_=1)
        : gfs(problem.model().jacobianAssembler().gridFunctionSpace()),
          phelper(gfs),
          maxiter(maxiter_),
          verbose(verbose_),
          constraintsTrafo_(problem.model().jacobianAssembler().constraintsTrafo()),
          exchanger_(problem)
    {}

    /*! \brief compute global norm of a vector

    \param[in] v the given vector
    */
    template<class Vector>
    typename Vector::ElementType norm (const Vector& v) const
    {
        Vector x(v); // make a copy because it has to be made consistent
        typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,Vector> PSP;
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
    template<class JacobianMatrix, class SolVector, class RhsVector>
    void apply(const JacobianMatrix& A, SolVector& z, RhsVector& r, typename SolVector::ElementType reduction)
    {
        typedef Dune::SeqILU0<JacobianMatrix,SolVector,RhsVector> SeqPreCond;
        JacobianMatrix B(A);
        exchanger_.sumEntries(B);
        SeqPreCond seqPreCond(B, 1.0);

        typedef Dune::PDELab::NonoverlappingOperator<GridFunctionSpace,JacobianMatrix,SolVector,RhsVector> POP;
        POP pop(gfs,A,phelper);
        typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,SolVector> PSP;
        PSP psp(gfs,phelper);
        typedef NonoverlappingWrappedPreconditioner<ConstraintsTrafo, GridFunctionSpace, SeqPreCond> ParPreCond;
        ParPreCond parPreCond(gfs, seqPreCond, constraintsTrafo_, exchanger_.borderIndices(), phelper);

        int verb=0;
        if (gfs.gridview().comm().rank()==0) verb=verbose;
        Dune::BiCGSTABSolver<SolVector> solver(pop,psp,parPreCond,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed = stat.elapsed;
        res.reduction = stat.reduction;
    }

    /*! \brief Return access to result data */
    const Dune::PDELab::LinearSolverResult<double>& result() const
                      {
        return res;
                      }

private:
    const GridFunctionSpace& gfs;
    PHELPER phelper;
    Dune::PDELab::LinearSolverResult<double> res;
    unsigned maxiter;
    int verbose;
    const ConstraintsTrafo& constraintsTrafo_;
    Exchanger<TypeTag> exchanger_;
};

/*!
 * \brief backend for an ISTL parallel Pardiso preconditioned loop solver
 */
template<class TypeTag>
class ISTLBackend_NoOverlap_Loop_Pardiso
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConstraintsTrafo)) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) JacobianMatrix;
    typedef Dune::PDELab::ParallelISTLHelper<GridFunctionSpace> PHELPER;

public:
    /*! \brief make a linear solver object

    \param[in] problem The Dumux problem
    \param[in] maxiter_ Maximum number of iterations to do
    \param[in] verbose_ Verbosity level
    */
    explicit ISTLBackend_NoOverlap_Loop_Pardiso (Problem& problem, unsigned maxiter_=5000, int verbose_=1)
        : gfs(problem.model().jacobianAssembler().gridFunctionSpace()),
          phelper(gfs),
          maxiter(maxiter_),
          verbose(verbose_),
          constraintsTrafo_(problem.model().jacobianAssembler().constraintsTrafo()),
          exchanger_(problem)
    {}

    /*! \brief compute global norm of a vector

    \param[in] v the given vector
    */
    template<class Vector>
    typename Vector::ElementType norm (const Vector& v) const
    {
        Vector x(v); // make a copy because it has to be made consistent
        typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,Vector> PSP;
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
    template<class M, class SolVector, class RhsVector>
    void apply(M& A, SolVector& z, RhsVector& r, typename SolVector::ElementType reduction)
    {
        typedef Dune::SeqPardiso<JacobianMatrix,SolVector,RhsVector> SeqPreCond;
        JacobianMatrix B(A);
        exchanger_.sumEntries(B);
        SeqPreCond seqPreCond(B);

        typedef Dune::PDELab::NonoverlappingOperator<GridFunctionSpace,JacobianMatrix,SolVector,RhsVector> POP;
        POP pop(gfs,A,phelper);
        typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,SolVector> PSP;
        PSP psp(gfs,phelper);
        typedef NonoverlappingWrappedPreconditioner<ConstraintsTrafo, GridFunctionSpace, SeqPreCond> ParPreCond;
        ParPreCond parPreCond(gfs, seqPreCond, constraintsTrafo_, exchanger_.borderIndices(), phelper);

        //        typedef Dune::PDELab::NonoverlappingRichardson<GridFunctionSpace,SolVector,RhsVector> PRICH;
        //        PRICH prich(gfs,phelper);
        int verb=0;
        if (gfs.gridview().comm().rank()==0) verb=verbose;
        Dune::LoopSolver<SolVector> solver(pop,psp,parPreCond,reduction,maxiter,verb);
//        Dune::BiCGSTABSolver<SolVector> solver(pop,psp,parPreCond,reduction,maxiter,verb);
        Dune::InverseOperatorResult stat;
        solver.apply(z,r,stat);
        res.converged = stat.converged;
        res.iterations = stat.iterations;
        res.elapsed = stat.elapsed;
        res.reduction = stat.reduction;
    }

    /*! \brief Return access to result data */
    const Dune::PDELab::LinearSolverResult<double>& result() const
                      {
        return res;
                      }

private:
    const GridFunctionSpace& gfs;
    PHELPER phelper;
    Dune::PDELab::LinearSolverResult<double> res;
    unsigned maxiter;
    int verbose;
    const ConstraintsTrafo& constraintsTrafo_;
    Exchanger<TypeTag> exchanger_;
};



} // namespace PDELab
} // namespace Dumux

#endif
