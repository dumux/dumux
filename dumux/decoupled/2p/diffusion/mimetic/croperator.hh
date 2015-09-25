// $Id$
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
#ifndef DUMUX_CROPERATOR_HH
#define DUMUX_CROPERATOR_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<cassert>
#include<stdio.h>
#include<stdlib.h>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>

#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dumux/common/boundaryconditions.hh>
#include<dumux/decoupled/2p/diffusion/mimetic/localstiffness.hh>

/**
 * @file
 * @brief  defines a class for piecewise linear finite element functions
 * @author Peter Bastian
 */

namespace Dumux
{
/*!
 * \ingroup Mimetic2p
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
/*! @brief Extends CROperatorBase by a generic methods to assemble global stiffness matrix from local stiffness matrices
 *
 *
 * The template parameter TypeTag describes what kind of Assembler we are. There two choices:
 * <dt>LevelTag</dt> We assemble on a grid level.
 * <dt>LeafTag</dt> We assemble on the leaf entities of the grid
 */
template<class Scalar, class GridView>
class CROperatorAssembler
{
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

    typedef typename GridView::Grid Grid;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Entity;
    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef typename GridView::IndexSet IS;
    typedef typename Grid::template Codim<0>::EntityPointer EEntityPointer;
    typedef Dune::FieldMatrix<Scalar,1,1> BlockType;
    typedef Dune::BCRSMatrix<BlockType> MatrixType;
    typedef typename MatrixType::block_type MBlockType;
    typedef typename MatrixType::RowIterator rowiterator;
    typedef typename MatrixType::ColIterator coliterator;
    typedef Dune::array<BoundaryConditions::Flags,1> BCBlockType;     // componentwise boundary conditions
    typedef Dune::BlockVector< Dune::FieldVector<double,1> > SatType;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,CRLayout> EM;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,AllLayout> AM;

    // return number of rows/columns
    int size () const
    {
        return faceMapper_.size();
    }

    //! a function to approximately compute the number of nonzeros
    int nnz (const IS& is)
    {
        return (4*dim - 1)*is.size(1);
    }

public:
    typedef MatrixType RepresentationType;

    CROperatorAssembler (const GridView& gridview)
    : gridView_(gridview), is_(gridView_.indexSet()), faceMapper_(gridView_), allMapper_(gridView_),
      A_(size(), size(), nnz(is_), RepresentationType::random)
    {
        assert(nnz(is_) != 0);

        // set size of all rows to zero
        for (unsigned int i = 0; i < is_.size(dim); i++)
            A_.setrowsize(i,0);

        // build needs a flag for all entities of all codims
        std::vector<bool> visited(allMapper_.size());
        for (int i = 0; i < allMapper_.size(); i++)
            visited[i] = false;

        // LOOP 1 : Compute row sizes
        Iterator eendit = gridView_.template end<0>();
        for (Iterator it = gridView_.template begin<0>(); it != eendit; ++it)
        {
            Dune::GeometryType gt = it->geometry().type();
            const typename Dune::GenericReferenceElementContainer<Scalar,dim>::value_type&
                refelem = Dune::GenericReferenceElements<Scalar,dim>::general(gt);

            // faces, c=1
            for (int i = 0; i < refelem.size(1); i++)
            {
                int index = allMapper_.map(*it, i,1);
                int alpha = faceMapper_.map(*it, i,1);
                //std::cout << "index = " << index << ", alpha = " << alpha << std::endl;
                if (!visited[index])
                {
                    A_.incrementrowsize(alpha);
                    visited[index] = true;
                    // printf("increment row %04d\n",alpha);
                }
                for (int k = 0; k < refelem.size(1)-1; k++) {
                    A_.incrementrowsize(alpha);
                    // printf("increment row %04d\n",alpha);
                }

            }

        }

        // now the row sizes have been set
        A_.endrowsizes();

        // clear the flags for the next round, actually that is not necessary because addindex takes care of this
        for (int i = 0; i < allMapper_.size(); i++)
            visited[i] = false;

        // LOOP 2 : insert the nonzeros
        for (Iterator it = gridView_.template begin<0>(); it!=eendit; ++it)
        {
            Dune::GeometryType gt = it->geometry().type();
            const typename Dune::GenericReferenceElementContainer<Scalar,dim>::value_type&
                refelem = Dune::GenericReferenceElements<Scalar,dim>::general(gt);
            //           std::cout << "ELEM " << GeometryName(gt) << std::endl;

            // faces, c=1
            for (int i = 0; i < refelem.size(1); i++)
            {
                int index = allMapper_.map(*it, i, 1);
                int alpha = faceMapper_.map(*it, i, 1);
                if (!visited[index])
                {
                    A_.addindex(alpha,alpha);
                    visited[index] = true;
                }
                for (int k = 0; k < refelem.size(1); k++)
                    if (k != i) {
                        int beta = faceMapper_.map(*it, k, 1);
                        A_.addindex(alpha, beta);
                        //std::cout << "alpha = " << alpha << ", added beta = " << beta << std::endl;
                    }

            }

        }

        // now the matrix is ready for use
        A_.endindices();
    }

    //! return const reference to operator matrix
    const RepresentationType& operator* () const
    {
        return A_;
    }

    //! return reference to operator matrix
    RepresentationType& operator* ()
    {
        return A_;
    }

    /*! @brief Assemble global stiffness matrix

      This method takes an object that can compute local stiffness matrices and
      assembles the global linear system Au=f.

      @param[in] loc the local assembler providing element stiffness and boundary conditions for all elements
      @param[in,out] u solution, contains initial values on input, Dirichlet values are set. The
      type of boundary condition for a node is inferred from the values returned
      by the local assembler. A node is of Neumann type if all elements referring
      to that node report a Neumann boundary condition, it is set to Dirichlet
      if a least one element reports a process or Dirichlet boundary condition. The difference
      between process and Dirichlet is that process always denotes a homogeneous Dirichlet
      value.
      @param[in] f right hand side is filled by this method

      Note that the rows corresponding to nodes at the Dirichlet boundary are filled
      with trivial equations of the form \f[1\cdot u_i = f_i \f] where \f$u_i\f$ and \f$f_i\f$ are both set to the
      Dirichlet value at the \f$i\f$th node.

    */
    template <class LocalStiffness, class Vector>
    void assemble (LocalStiffness& loc, Vector& u, Vector& f)
    {

        // check size
        if (u.N()!=A_.M() || f.N()!=A_.N())
            DUNE_THROW(Dune::MathError,"CROperatorAssembler::assemble(): size mismatch");
        // clear global stiffness matrix and right hand side
        A_ = 0;
        f = 0;

        // allocate flag vector to hold flags for essential boundary conditions
        std::vector<BCBlockType> essential(faceMapper_.size());
        for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++)
                essential[i][0] = BoundaryConditions::neumann;

        // local to global id mapping (do not ask vertex mapper repeatedly
        int local2Global[2*GridView::dimension];

        // run over all leaf elements
        Iterator eendit = gridView_.template end<0>();
        for (Iterator it = gridView_.template begin<0>(); it!=eendit; ++it)
        {
            unsigned int numFaces = it->template count<1>();

            // get local to global id map
            for (int k = 0; k < numFaces; k++)
            {
                int alpha = faceMapper_.map(*it, k, 1);
                local2Global[k] = alpha;
            }

            // build local stiffness matrix for CR elements
            // inludes rhs and boundary condition information
            loc.assemble(*it, 1); // assemble local stiffness matrix


            // accumulate local matrix into global matrix for non-hanging nodes
            for (int i=0; i<numFaces; i++) // loop over rows, i.e. test functions
            {
                // accumulate matrix
                for (int j=0; j<numFaces; j++)
                {
                    // the standard entry
                    A_[local2Global[i]][local2Global[j]] += loc.mat(i,j);
                }

                // essential boundary condition and rhs
                if (loc.bc(i)[0]>essential[local2Global[i]][0])
                {
                    essential[local2Global[i]][0] = loc.bc(i)[0];
                    f[local2Global[i]][0] = loc.rhs(i)[0];
                }
                if (essential[local2Global[i]][0]==BoundaryConditions::neumann)
                    f[local2Global[i]][0] += loc.rhs(i)[0];
            }

        }

        // put in essential boundary conditions
        rowiterator endi=A_.end();
        for (rowiterator i=A_.begin(); i!=endi; ++i)
        {
            // muck up extra rows
            if ((int) i.index() >= (int) faceMapper_.size())
            {
                coliterator endj=(*i).end();
                for (coliterator j=(*i).begin(); j!=endj; ++j)
                {
                    (*j) = 0;
                    if (j.index()==i.index())
                        (*j)[0][0] = 1;
                }
                f[i.index()] = 0;
                continue;
            }

            // insert dirichlet ans processor boundary conditions
            if (essential[i.index()][0]!=BoundaryConditions::neumann)
            {
                coliterator endj=(*i).end();
                for (coliterator j=(*i).begin(); j!=endj; ++j)
                    if (j.index()==i.index())
                    {
                        (*j)[0][0] = 1;
                    }
                    else
                    {
                        (*j)[0][0] = 0;
                    }
                u[i.index()][0] = f[i.index()][0];
            }
        }
    }

protected:
    const GridView& gridView_;
    const IS& is_;
    EM faceMapper_;
    AM allMapper_;
    RepresentationType A_;
};

}

#endif
