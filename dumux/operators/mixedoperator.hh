// $Id$

#ifndef DUNE_MIXEDOPERATOR_HH
#define DUNE_MIXEDOPERATOR_HH

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
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include"dumux/functions/mixedfunction.hh"
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

/*! @brief A class for mapping a Mixed function to a Mixed function

  This class sets up a compressed row storage matrix with connectivity for Mixed elements.

  This class does not fill any entries into the matrix.

  The template parameter TypeTag describes what kind of Assembler we are. There two choices:
  <dt>LevelTag</dt> We assemble on a grid level.
  <dt>LeafTag</dt> We assemble on the leaf entities of the grid
*/
template<typename TypeTag, class G, class RT, class GV, class LC, int m=1>
class MixedOperatorBase
{
public:
    // export type used to store the matrix
    typedef FieldMatrix<RT,m,m> BlockType;
    typedef BCRSMatrix<BlockType> RepresentationType;

    // mapper: one data element per face
    template<int dim>
    struct FaceLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim-1;
        }
    };

    // mapper: one data element per element
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    // mapper: one data element per element
    template<int dim>
    struct ElementAndFaceLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return (gt.dim() == dim || gt.dim() == dim-1);
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

    typedef typename G::ctype Scalar;
    enum {dim=G::dimension};
    typedef typename GV::IndexSet IS;
    typedef typename G::template Codim<0>::Entity Element;
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef typename G::template Codim<0>::EntityPointer ElementPointer;
    typedef typename G::Traits::GlobalIdSet IDS;
    typedef typename IDS::IdType IdType;
    typedef std::set<IdType> GIDSet;
    typedef MultipleCodimMultipleGeomTypeMapper<G,IS,FaceLayout> FaceMapper;
    typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> ElementMapper;
    typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementAndFaceLayout> ElementAndFaceMapper;
    typedef MultipleCodimMultipleGeomTypeMapper<G,IS,AllLayout> AllMapper;

private:

    //! a function to approximately compute the number of nonzeros
    int nnz (const IS& is)
    {
        return (((4*dim - 1) + 2)*is.size(1) + (3*dim + 1)*is.size(0));
    }


    // return number of rows/columns
    int size () const
    {
        return elementAndFaceMapper.size();
    }

public:

    MixedOperatorBase (const G& g, const GV& gv, LC lcomm)
        : grid(g), gridview(gv), is(gv.indexSet()), lc(lcomm), faceMapper(g,is), elementMapper(g, is), elementAndFaceMapper(g, is), allMapper(g,is),
          A(size(), size(), nnz(is), RepresentationType::random)
    {
        // be verbose
        std::cout << g.comm().rank() << ": " << "making " << size() << "x"
                  << size() << " matrix with " << nnz(is) << " nonzeros" << std::endl;

        // set size of all rows to zero
        for (int i = 0; i < size(); i++)
            A.setrowsize(i,0);

        // build needs a flag for all entities of all codims
        std::vector<bool> visited(allMapper.size());
        for (int i = 0; i < allMapper.size(); i++)
            visited[i] = false;

        // LOOP 1 : Compute row sizes
        watch.reset();
        ElementIterator endEIt = grid.template leafend<0>();
        for (ElementIterator eIt = grid.template leafbegin<0>(); eIt != endEIt; ++eIt)
        {
            const Element& element = *eIt;

            int eIdx = elementAndFaceMapper.map(element);
            int nFaces = element.template count<1>();
            A.incrementrowsize(eIdx, nFaces + 1);

            for (int i = 0; i < nFaces; i++)
            {
                int index = allMapper.template map<1>(element, i);
                int fIdx = elementAndFaceMapper.template map<1>(element, i);

                if (!visited[index])
                {
                    A.incrementrowsize(fIdx);
                    visited[index] = true;
                }
                A.incrementrowsize(fIdx, nFaces);
            }
        }

        // now the row sizes have been set
        A.endrowsizes();
        std::cout << "=== MixedOperatorBase compute row sizes " <<  watch.elapsed() << std::endl;

        // clear the flags for the next round, actually that is not necessary because addindex takes care of this
        for (int i = 0; i < allMapper.size(); i++)
            visited[i] = false;

        // LOOP 2 : insert the nonzeros
        watch.reset();
        for (ElementIterator eIt = grid.template leafbegin<0>(); eIt!=endEIt; ++eIt)
        {
            const Element& element = *eIt;

            int eIdx = elementAndFaceMapper.map(element);
            A.addindex(eIdx, eIdx);
            int nFaces = element.template count<1>();

            for (int i = 0; i < nFaces; i++)
            {
                int index = allMapper.template map<1>(element, i);
                int fIdx = elementAndFaceMapper.template map<1>(element, i);
                if (!visited[index])
                {
                    A.addindex(fIdx,fIdx);
                    visited[index] = true;
                }
                for (int k = 0; k < nFaces; k++)
                    if (k != i) {
                        int fIdx2 = elementAndFaceMapper.template map<1>(element, k);
                        A.addindex(fIdx, fIdx2);
                    }

                A.addindex(fIdx, eIdx);
                A.addindex(eIdx, fIdx);
            }
        }
        A.endindices();
        std::cout << "=== MixedOperatorBase index insertion " <<  watch.elapsed() << std::endl;
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
    FaceMapper faceMapper;
    ElementMapper elementMapper;
    ElementAndFaceMapper elementAndFaceMapper;
    AllMapper allMapper;
    RepresentationType A;
};




/*! @brief Extends MixedOperatorBase by a generic methods to assemble global stiffness matrix from local stiffness matrices
 *
 *
 * The template parameter TypeTag describes what kind of Assembler we are. There two choices:
 * <dt>LevelTag</dt> We assemble on a grid level.
 * <dt>LeafTag</dt> We assemble on the leaf entities of the grid
 */
template<typename TypeTag, class G, class RT, class GV, class LC, int m>
class MixedOperatorAssembler : public MixedOperatorBase<TypeTag,G,RT,GV,LC,m>
{
    typedef typename G::ctype Scalar;
    enum {n=G::dimension};
    typedef typename G::template Codim<0>::Entity Element;
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef typename GV::IndexSet IS;
    typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
    typedef typename G::template Codim<0>::EntityPointer ElementPointer;
    typedef typename MixedFunction<G,RT,GV,LC,m>::RepresentationType VectorType;
    typedef typename VectorType::block_type VBlockType;
    typedef typename MixedOperatorBase<TypeTag,G,RT,GV,LC,m>::RepresentationType MatrixType;
    typedef typename MatrixType::block_type MBlockType;
    typedef typename MatrixType::RowIterator rowiterator;
    typedef typename MatrixType::ColIterator coliterator;
    typedef array<BoundaryConditions::Flags,m> BCBlockType;     // componentwise boundary conditions
    typedef Dune::BlockVector< Dune::FieldVector<double,1> > SatType;

public:
    MixedOperatorAssembler (const G& g, const GV& gridview, LC lcomm)
        : MixedOperatorBase<TypeTag,G,RT,GV,LC,m>(g,gridview,lcomm)
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
    template<class LocalStiffness>
    void assemble (LocalStiffness& loc,
                   MixedFunction<G,RT,GV,LC,m>& u,
                   MixedFunction<G,RT,GV,LC,m>& f)
    {

        // check size
        //        if ((*u).N()!=this->A.M() || (*f).N()!=this->A.N())
        //            DUNE_THROW(MathError,"MixedOperatorAssembler::assemble(): size mismatch");

        // clear global stiffness matrix and right hand side
        this->watch.reset();
        this->A = 0;
        *f = 0;
        //       std::cout << "=== MixedOperatorBase clear matrix " <<  this->watch.elapsed() << std::endl;
        // allocate flag vector to hold flags for essential boundary conditions
        std::vector<BCBlockType> essential(this->elementAndFaceMapper.size());
        for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++)
            essential[i].assign(BoundaryConditions::neumann);

        // local to global id mapping (do not ask vertex mapper repeatedly
        int local2Global[Dune::CRShapeFunctionSetContainer<Scalar,RT,n>::maxsize+1];

        // run over all leaf elements
        ElementIterator endEIt = this->grid.template leafend<0>();
        for (ElementIterator eIt = this->grid.template leafbegin<0>(); eIt!=endEIt; ++eIt)
        {
            const Element& element = *eIt;

            int nFaces = element.template count<1>();
            int nDOF = nFaces + 1;

            // get local to global id map
            for (int k = 0; k < nFaces; k++)
                local2Global[k] = this->elementAndFaceMapper.template map<1>(element, k);

            local2Global[nFaces] = this->elementAndFaceMapper.map(element);

            // build local stiffness matrix for Mixed elements
            // inludes rhs and boundary condition information
            loc.assemble(element, 1); // assemble local stiffness matrix


            // accumulate local matrix into global matrix for non-hanging nodes
            for (int i=0; i<nDOF; i++) // loop over rows, i.e. test functions
            {
                // accumulate matrix
                for (int j=0; j<nDOF; j++)
                {
                    // the standard entry
                    this->A[local2Global[i]][local2Global[j]] += loc.mat(i,j);
                }

                // essential boundary condition and rhs
                for (int comp = 0; comp < m; comp++)
                {
                    if (loc.bc(i)[comp] > essential[local2Global[i]][comp])
                    {
                        essential[local2Global[i]][comp] = loc.bc(i)[comp];
                        (*f)[local2Global[i]][comp] = loc.rhs(i)[comp];
                    }
                    if (essential[local2Global[i]][comp] == BoundaryConditions::neumann)
                        (*f)[local2Global[i]][comp] += loc.rhs(i)[comp];
                }
            }

        }

        // put in essential boundary conditions
        rowiterator endi = this->A.end();
        for (rowiterator i = this->A.begin(); i!=endi; ++i)
        {
            // insert dirichlet ans processor boundary conditions
            for (int icomp = 0; icomp < m; icomp++)
                if (essential[i.index()][icomp] != BoundaryConditions::neumann
                    || (i.index() == 0 && icomp == 0)) // the pressure is set to zero in the first element
                {
                    coliterator endj = (*i).end();
                    for (coliterator j = (*i).begin(); j != endj; ++j)
                        if (j.index() == i.index())
                        {
                            for (int jcomp = 0; jcomp < m; jcomp++)
                                if (icomp == jcomp)
                                    (*j)[icomp][jcomp] = 1;
                                else
                                    (*j)[icomp][jcomp] = 0;
                        }
                        else
                        {
                            for (int jcomp = 0; jcomp < m; jcomp++)
                                (*j)[icomp][jcomp] = 0;
                        }
                    (*u)[i.index()][icomp] = (*f)[i.index()][icomp];
                }
        }
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
class LeafMixedOperatorAssembler : public MixedOperatorAssembler<LeafTag,G,RT,typename G::LeafGridView,LeafCommunicate<G>,m>
{
public:
    LeafMixedOperatorAssembler (const G& grid, bool overlap = false)
        : MixedOperatorAssembler<LeafTag,G,RT,typename G::LeafGridView,LeafCommunicate<G>,m>(grid,grid.leafView(),LeafCommunicate<G>(grid))
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
class LevelMixedOperatorAssembler : public MixedOperatorAssembler<LevelTag,G,RT,typename G::LevelGridView,LevelCommunicate<G>,m>
{
public:
    LevelMixedOperatorAssembler (const G& grid, int level)
        : MixedOperatorAssembler<LevelTag,G,RT,typename G::LevelGridView,LevelCommunicate<G>,m>(grid,grid.levelView(level),LevelCommunicate<G>(grid,level))
    {}
};


/** @} */

}
#endif
