// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_CROPERATOR_HH
#define DUNE_CROPERATOR_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>

//#include<dune/grid/common/datahandleif.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/operators/localstiffness.hh>
#include"dumux/functions/CRfunction.hh"
#include"dumux/shapefunctions/CRshapefunctions.hh"

/**
 * @file
 * @brief  defines a class for piecewise linear finite element functions
 * @author Peter Bastian
 */

/*! @defgroup DISC_Operators Operators
  @ingroup DISC
  @brief

  @section D1 Introduction
  <!--=================-->

  To be written
*/

namespace Dune
{
/** @addtogroup DISC_Operators
 *
 * @{
 */
/**
 * @brief defines a class for Crozieux-Raviart piecewise linear finite element functions
 *
 */

/*! @brief A class for mapping a CR function to a CR function

  This class sets up a compressed row storage matrix with connectivity for CR elements.

  This class does not fill any entries into the matrix.

  The template parameter TypeTag describes what kind of Assembler we are. There two choices:
  <dt>LevelTag</dt> We assemble on a grid level.
  <dt>LeafTag</dt> We assemble on the leaf entities of the grid
*/
template<class G, class RT, class GV, class LC, int m=1>
class CROperatorBase
{
public:
    // export type used to store the matrix
    typedef FieldMatrix<RT,m,m> BlockType;
    typedef BCRSMatrix<BlockType> RepresentationType;

    // mapper: one data element per vertex
    template<int dim>
    struct CRLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim-1;
        }
    };

    // mapper: one data element in every entity
    template<int dim>
    struct AllLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return true;
        }
    };

    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename GV::IndexSet IS;
    typedef typename G::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename G::template Codim<0>::EntityPointer EEntityPointer;
    typedef typename G::Traits::GlobalIdSet IDS;
    typedef typename IDS::IdType IdType;
    typedef std::set<IdType> GIDSet;
    typedef MultipleCodimMultipleGeomTypeMapper<GV,CRLayout> EM;
    typedef MultipleCodimMultipleGeomTypeMapper<GV,AllLayout> AM;

private:

    //! a function to approximately compute the number of nonzeros
    int nnz (const IS& is)
    {
        return (4*n - 1)*is.size(1);
    }


    // return number of rows/columns
    int size () const
    {
        return facemapper.size();
    }

    struct MatEntry
    {
        IdType first;
        BlockType second;
        MatEntry (const IdType& f, const BlockType& s) : first(f),second(s) {}
        MatEntry () {}
    };

public:

    CROperatorBase (const G& g, const GV& gv, LC lcomm)
        : grid(g), gridview(gv), is(gv.indexSet()), lc(lcomm), facemapper(gv), allmapper(gv),
          A(size(), size(), nnz(is), RepresentationType::random)
    {
        // be verbose
        std::cout << g.comm().rank() << ": " << "making " << size() << "x"
                  << size() << " matrix with " << nnz(is) << " nonzeros" << std::endl;

        // set size of all rows to zero
        for (unsigned int i = 0; i < is.size(n); i++)
            A.setrowsize(i,0);

        // build needs a flag for all entities of all codims
        std::vector<bool> visited(allmapper.size());
        for (int i = 0; i < allmapper.size(); i++)
            visited[i] = false;

        // LOOP 1 : Compute row sizes
        watch.reset();
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
        {
            Dune::GeometryType gt = it->geometry().type();
            const typename Dune::ReferenceElementContainer<DT,n>::value_type&
                refelem = ReferenceElements<DT,n>::general(gt);

            // faces, c=1
            for (int i = 0; i < refelem.size(1); i++)
            {
                int index = allmapper.map(*it, i,1);
                int alpha = facemapper.map(*it, i,1);
                //std::cout << "index = " << index << ", alpha = " << alpha << std::endl;
                if (!visited[index])
                {
                    A.incrementrowsize(alpha);
                    visited[index] = true;
                    // printf("increment row %04d\n",alpha);
                }
                for (int k = 0; k < refelem.size(1)-1; k++) {
                    A.incrementrowsize(alpha);
                    // printf("increment row %04d\n",alpha);
                }

            }

        }

        // now the row sizes have been set
        A.endrowsizes();
        std::cout << "=== CROperatorBase compute row sizes " <<  watch.elapsed() << std::endl;

        // clear the flags for the next round, actually that is not necessary because addindex takes care of this
        for (int i = 0; i < allmapper.size(); i++)
            visited[i] = false;

        // LOOP 2 : insert the nonzeros
        watch.reset();
        for (Iterator it = gridview.template begin<0>(); it!=eendit; ++it)
        {
            Dune::GeometryType gt = it->geometry().type();
            const typename Dune::ReferenceElementContainer<DT,n>::value_type&
                refelem = ReferenceElements<DT,n>::general(gt);
            //           std::cout << "ELEM " << GeometryName(gt) << std::endl;

            // faces, c=1
            for (int i = 0; i < refelem.size(1); i++)
            {
                int index = allmapper.map(*it, i, 1);
                int alpha = facemapper.map(*it, i, 1);
                if (!visited[index])
                {
                    A.addindex(alpha,alpha);
                    visited[index] = true;
                }
                for (int k = 0; k < refelem.size(1); k++)
                    if (k != i) {
                        int beta = facemapper.map(*it, k, 1);
                        A.addindex(alpha, beta);
                        //std::cout << "alpha = " << alpha << ", added beta = " << beta << std::endl;
                    }

            }

        }

        // now the matrix is ready for use
        A.endindices();
        std::cout << "=== CROperatorBase index insertion " <<  watch.elapsed() << std::endl;
    }

    //! return const reference to operator matrix
    const RepresentationType& operator* () const
    {
        return A;
    }

    //! return reference to operator matrix
    RepresentationType& operator* ()
    {
        return A;
    }

protected:
    Timer watch;
    const G& grid;
    const GV& gridview;
    const IS& is;
    LC lc;
    EM facemapper;
    AM allmapper;
    RepresentationType A;
};




/*! @brief Extends CROperatorBase by a generic methods to assemble global stiffness matrix from local stiffness matrices
 *
 *
 * The template parameter TypeTag describes what kind of Assembler we are. There two choices:
 * <dt>LevelTag</dt> We assemble on a grid level.
 * <dt>LeafTag</dt> We assemble on the leaf entities of the grid
 */
template<class G, class RT, class GV, class LC, int m>
class CROperatorAssembler : public CROperatorBase<G,RT,GV,LC,m>
{
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename GV::IndexSet IS;
    typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
    typedef typename G::template Codim<0>::EntityPointer EEntityPointer;
    typedef typename CRFunction<G,RT,GV,LC,m>::RepresentationType VectorType;
    typedef typename VectorType::block_type VBlockType;
    typedef typename CROperatorBase<G,RT,GV,LC,m>::RepresentationType MatrixType;
    typedef typename MatrixType::block_type MBlockType;
    typedef typename MatrixType::RowIterator rowiterator;
    typedef typename MatrixType::ColIterator coliterator;
    typedef array<BoundaryConditions::Flags,m> BCBlockType;     // componentwise boundary conditions
    typedef Dune::BlockVector< Dune::FieldVector<double,1> > SatType;

public:
    CROperatorAssembler (const G& g, const GV& gridview, LC lcomm)
        : CROperatorBase<G,RT,GV,LC,m>(g,gridview,lcomm)
    {    }




    /*! @brief Assemble global stiffness matrix

      This method takes an object that can compute local stiffness matrices and
      assembles the global linear system Au=f.

      @param[in] loc    the local assembler providing element stiffness and boundary conditions for all elements
      @param[in,out] u  solution, contains initial values on input, Dirichlet values are set. The
      type of boundary condition for a node is inferred from the values returned
      by the local assembler. A node is of Neumann type if all elements referring
      to that node report a Neumann boundary condition, it is set to Dirichlet
      if a least one element reports a process or Dirichlet boundary condition. The difference
      between process and Dirichlet is that process always denotes a homogeneous Dirichlet
      value.
      @param[in] f      right hand side is filled by this method

      Note that the rows corresponding to nodes at the Dirichlet boundary are filled
      with trivial equations of the form \f[1\cdot u_i = f_i \f] where \f$u_i\f$ and \f$f_i\f$ are both set to the
      Dirichlet value at the \f$i\f$th node.

    */
    void assemble (LocalStiffness<GV,RT,m>& loc,
                   CRFunction<G,RT,GV,LC,m>& u,
                   CRFunction<G,RT,GV,LC,m>& f)
    {

        // check size
        if ((*u).N()!=this->A.M() || (*f).N()!=this->A.N())
            DUNE_THROW(MathError,"CROperatorAssembler::assemble(): size mismatch");
        // clear global stiffness matrix and right hand side
        this->watch.reset();
        this->A = 0;
        *f = 0;
        //       std::cout << "=== CROperatorBase clear matrix " <<  this->watch.elapsed() << std::endl;
        // allocate flag vector to hold flags for essential boundary conditions
        std::vector<BCBlockType> essential(this->facemapper.size());
        for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++)
            essential[i].assign(BoundaryConditions::neumann);

        // local to global id mapping (do not ask vertex mapper repeatedly
        int local2Global[Dune::CRShapeFunctionSetContainer<DT,RT,n>::maxsize];

        // run over all leaf elements
        Iterator eendit = this->gridview.template end<0>();
        for (Iterator it = this->gridview.template begin<0>(); it!=eendit; ++it)
        {
            // get access to shape functions for CR elements
            Dune::GeometryType gt = it->geometry().type();
            const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
                sfs = Dune::CRShapeFunctions<DT,RT,n>::general(gt,1);

            // get local to global id map
            for (int k = 0; k < sfs.size(); k++)
            {
                if (sfs[k].codim() != 1) DUNE_THROW(MathError, "expected codim == dim");
                int alpha = this->facemapper.map(*it, k, 1);
                local2Global[k] = alpha;
                //FieldVector<double,n> global = (*it).geometry().global(sfs[k].position());
                //std::cout << "local = " << sfs[k].position() << ", global = " << global << " -> " << alpha << std::endl;
            }

            // build local stiffness matrix for CR elements
            // inludes rhs and boundary condition information
            loc.assemble(*it, 1); // assemble local stiffness matrix


            // accumulate local matrix into global matrix for non-hanging nodes
            for (int i=0; i<sfs.size(); i++) // loop over rows, i.e. test functions
            {
                // accumulate matrix
                for (int j=0; j<sfs.size(); j++)
                {
                    // the standard entry
                    this->A[local2Global[i]][local2Global[j]] += loc.mat(i,j);
                }

                // essential boundary condition and rhs
                for (int comp=0; comp<m; comp++)
                {
                    if (loc.bc(i)[comp]>essential[local2Global[i]][comp])
                    {
                        essential[local2Global[i]][comp] = loc.bc(i)[comp];
                        (*f)[local2Global[i]][comp] = loc.rhs(i)[comp];
                    }
                    if (essential[local2Global[i]][comp]==BoundaryConditions::neumann)
                        (*f)[local2Global[i]][comp] += loc.rhs(i)[comp];
                }
            }

        }

        // put in essential boundary conditions
        rowiterator endi=this->A.end();
        for (rowiterator i=this->A.begin(); i!=endi; ++i)
        {
            // muck up extra rows
            if ((int) i.index() >= (int) this->facemapper.size())
            {
                coliterator endj=(*i).end();
                for (coliterator j=(*i).begin(); j!=endj; ++j)
                {
                    (*j) = 0;
                    if (j.index()==i.index())
                        for (int comp=0; comp<m; comp++)
                            (*j)[comp][comp] = 1;
                }
                (*f)[i.index()] = 0;
                continue;
            }

            // insert dirichlet ans processor boundary conditions
            for (int icomp=0; icomp<m; icomp++)
                if (essential[i.index()][icomp]!=BoundaryConditions::neumann)
                {
                    coliterator endj=(*i).end();
                    for (coliterator j=(*i).begin(); j!=endj; ++j)
                        if (j.index()==i.index())
                        {
                            for (int jcomp=0; jcomp<m; jcomp++)
                                if (icomp==jcomp)
                                    (*j)[icomp][jcomp] = 1;
                                else
                                    (*j)[icomp][jcomp] = 0;
                        }
                        else
                        {
                            for (int jcomp=0; jcomp<m; jcomp++)
                                (*j)[icomp][jcomp] = 0;
                        }
                    (*u)[i.index()][icomp] = (*f)[i.index()][icomp];
                }
        }
    }

    void preMark ()
    {
        marked.resize(this->facemapper.size());
        for (std::size_t i=0; i<marked.size(); i++) marked[i] = false;
        return;
    }

    void postMark (G& g)
    {
        // run over all leaf elements
        int extra=0;
        Iterator eendit = this->gridview.template end<0>();
        for (Iterator it = this->gridview.template begin<0>(); it!=eendit; ++it)
        {
            // get access to shape functions for CR elements
            Dune::GeometryType gt = it->geometry().type();
            //          if (gt!=Dune::simplex && gt!=Dune::triangle && gt!=Dune::tetrahedron) continue;

            const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
                sfs=Dune::CRShapeFunctions<DT,RT,n>::general(gt,1);

            // count nodes with mark
            int count=0;
            for (int k=0; k<sfs.size(); k++)
            {
                int alpha = this->facemapper.template map<n>(*it, k);
                if (marked[alpha]) count++;
            }

            // refine if a marked edge exists
            if (count>0) {
                extra++;
                g.mark(1,it);
            }
        }

        //       std::cout << "placed " << extra << " extra marks" << std::endl;
        marked.clear();
        return;
    }

    void mark (G& g, EEntityPointer& it)
    {
        // refine this element
        g.mark(1,it);

        // check geom type, exit if not simplex
        Dune::GeometryType gt = it->geometry().type();
        //if (gt!=Dune::simplex && gt!=Dune::triangle && gt!=Dune::tetrahedron) return;
        assert(it->isLeaf());

        return;

    }


private:
    std::vector<bool> marked;

    template<class M, class K>
    void accumulate (M& A, const K& alpha, const M& B)
    {
        for (std::size_t i=0; i<A.N(); i++)
            for (std::size_t j=0; j<A.M(); j++)
                A[i][j] += alpha*B[i][j];
    }
};


/*! @brief Leafwise assembler

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:

  - G    A grid type
  - RT   The field type used in the elements of the stiffness matrix
  - m    number of degrees of freedom per node (system size)
*/
template<class G, class RT, int m>
class LeafCROperatorAssembler : public CROperatorAssembler<G,RT,typename G::LeafGridView,LeafCommunicate<G>,m>
{
public:
    LeafCROperatorAssembler (const G& grid)
        : CROperatorAssembler<G,RT,typename G::LeafGridView,LeafCommunicate<G>,m>(grid,grid.leafView(),LeafCommunicate<G>(grid))
    {}
};


/*! @brief Levelwise assembler

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:

  - G    A grid type
  - RT   The field type used in the elements of the stiffness matrix
  - m    number of degrees of freedom per node (system size)
*/
template<class G, class RT, int m>
class LevelCROperatorAssembler : public CROperatorAssembler<G,RT,typename G::LevelGridView,LevelCommunicate<G>,m>
{
public:
    LevelCROperatorAssembler (const G& grid, int level)
        : CROperatorAssembler<G,RT,typename G::LevelGridView,LevelCommunicate<G>,m>(grid,grid.levelView(level),LevelCommunicate<G>(grid,level))
    {}
};


/** @} */

}
#endif
