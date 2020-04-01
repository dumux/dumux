// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SequentialTwoPModel
 * \brief Defines a class for Crozieux-Raviart piecewise linear finite element functions.
 */
#ifndef DUMUX_CROPERATOR2P_HH
#define DUMUX_CROPERATOR2P_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<cassert>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/geometry/type.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>

#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/common/boundaryconditions.hh>
#include "localstiffness.hh"

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Extends CROperatorBase by a generic methods to assemble global stiffness matrix from local stiffness matrices.
 *
 * A class for mapping a CR function to a CR function
 * This class sets up a compressed row storage matrix with connectivity for CR elements.
 * This class does not fill any entries into the matrix.
 *
 * The template parameter TypeTag describes what kind of Assembler we are. There two choices:
 * <dt>LevelTag</dt> We assemble on a grid level.
 * <dt>LeafTag</dt> We assemble on the leaf entities of the grid
 *
 * \tparam TypeTag The problem Type Tag
 */
template<class TypeTag>
class CROperatorAssemblerTwoP
{
    template<int dim>
    struct FaceLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim-1;
        }
    };
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    enum {dim=GridView::dimension};
    using IS = typename GridView::IndexSet;
    using BlockType = Dune::FieldMatrix<Scalar, 1, 1>;
    using MatrixType = Dune::BCRSMatrix<BlockType>;
    using MBlockType = typename MatrixType::block_type;
    using rowiterator = typename MatrixType::RowIterator;
    using coliterator = typename MatrixType::ColIterator;
    using BCBlockType = std::array<BoundaryConditions::Flags, 1>;     // componentwise boundary conditions
    using SatType = Dune::BlockVector< Dune::FieldVector<double, 1> >;
    using FaceMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum
    {
        pressureEqIdx = Indices::pressureEqIdx,
    };

    //! A function to approximately compute the number of nonzeros
    int nnz()
    {
        return (4*dim - 1)*size_;
    }

public:
    using RepresentationType = MatrixType;

    CROperatorAssemblerTwoP (const GridView& gridview)
    : gridView_(gridview)
    , faceMapper_(gridView_, Dune::mcmgLayout(Dune::Codim<1>()))
    , size_(faceMapper_.size())
    , A_(size_, size_, nnz()
    , RepresentationType::random)
    {}

    //! Initialize the CR operator assembler
    void initialize()
    {
        faceMapper_.update();

        size_ = faceMapper_.size();
        A_.setSize(size_, size_,  nnz());

        assert(nnz() != 0);

        // set size of all rows to zero
        for (unsigned int i = 0; i < size_; i++)
            A_.setrowsize(i,0);

        // build needs a flag for all entities of all codims
        std::vector<bool> visited(size_, false);

        // LOOP 1 : Compute row sizes
        for (const auto& element : elements(gridView_))
        {
            int numFaces = element.subEntities(1);

            for (int i = 0; i < numFaces; i++)
            {
                int index = faceMapper_.subIndex(element, i,1);

                if (!visited[index])
                {
                    A_.incrementrowsize(index);
                    visited[index] = true;
                    // std::cout << "increment row " << index << std::endl;
                }
                A_.incrementrowsize(index, numFaces - 1);
                // std::cout << "increment row " << index
                // << " by " << numFaces - 1 << std::endl;
            }
        }

        // now the row sizes have been set
        A_.endrowsizes();

        // clear the flags for the next round, actually that is not necessary because addindex takes care of this
        visited.assign(size_, false);

        // LOOP 2 : insert the nonzeros
        for (const auto& element : elements(gridView_))
        {
            int numFaces = element.subEntities(1);

            for (int i = 0; i < numFaces; i++)
            {
                int indexI = faceMapper_.subIndex(element, i, 1);

                if (!visited[indexI])
                {
                    A_.addindex(indexI,indexI);
                    visited[indexI] = true;
                }
                for (int k = 0; k < numFaces; k++)
                    if (k != i) {
                        int indexJ = faceMapper_.subIndex(element, k, 1);

                        A_.addindex(indexI, indexJ);
                        //std::cout << "indexI = " << indexI << ", added indexJ = " << indexJ << std::endl;
                    }
            }
        }

        // now the matrix is ready for use
        A_.endindices();
    }

    //! Returns const reference to operator matrix
    const RepresentationType& operator* () const
    {
        return A_;
    }

    //! Returns a reference to operator matrix
    RepresentationType& operator* ()
    {
        return A_;
    }

    const FaceMapper& faceMapper()
    {
        return faceMapper_;
    }

    const FaceMapper& faceMapper() const
    {
        return faceMapper_;
    }

    /*!
     * \brief Assembles global stiffness matrix
     *
     * This method takes an object that can compute local stiffness matrices and
     * assembles the global linear system Au=f.
     *
     * \param loc the local assembler providing element stiffness and boundary conditions for all elements
     * \param u solution, contains initial values on input, Dirichlet values are set. The type of boundary condition
     * for a node is inferred from the values returned by the local assembler. A node is of Neumann type if all
     * elements referring to that node report a Neumann boundary condition, it is set to Dirichlet if a least one
     * element reports a process or Dirichlet boundary condition. The difference between process and Dirichlet
     * is that process always denotes a homogeneous Dirichlet value.
     * \param f right hand side is filled by this method
     *
     * Note that the rows corresponding to nodes at the Dirichlet boundary are filled with trivial equations
     * of the form \f[1\cdot u_i = f_i \f] where \f$u_i\f$ and \f$f_i\f$ are both set to the Dirichlet value at the \f$i\f$th node.
     */
    template <class LocalStiffness, class Vector>
    void assemble (LocalStiffness& loc, Vector& u, Vector& f)
    {

        // check size
        if (u.N()!=A_.M() || f.N()!=A_.N())
            DUNE_THROW(Dune::MathError,"CROperatorAssemblerTwoP::assemble(): size mismatch");
        // clear global stiffness matrix and right hand side
        A_ = 0;
        f = 0;

        // allocate flag vector to hold flags for essential boundary conditions
        std::vector<BCBlockType> essential(faceMapper_.size());
        for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++)
                essential[i][0] = BoundaryConditions::neumann;

        // local to global id mapping (do not ask vertex mapper repeatedly
        Dune::FieldVector<int, 2*dim> local2Global(0);

        // run over all leaf elements
        for (const auto& element : elements(gridView_))
        {
            unsigned int numFaces = element.subEntities(1);

            // get local to global id map
            for (unsigned int k = 0; k < numFaces; k++)
            {
                int alpha = faceMapper_.subIndex(element, k, 1);
                local2Global[k] = alpha;
            }

            // build local stiffness matrix for CR elements
            // inludes rhs and boundary condition information
            loc.assemble(element, 1); // assemble local stiffness matrix


            // accumulate local matrix into global matrix for non-hanging nodes
            for (unsigned int i=0; i<numFaces; i++) // loop over rows, i.e. test functions
            {
                // accumulate matrix
                for (unsigned int j=0; j<numFaces; j++)
                {
                    // the standard entry
                    A_[local2Global[i]][local2Global[j]] += loc.mat(i,j);
                }

                // essential boundary condition and rhs
                if (loc.bc(i).isDirichlet(pressureEqIdx))
                {
                    essential[local2Global[i]][0] = BoundaryConditions::dirichlet;
                    f[local2Global[i]][0] = loc.rhs(i)[0];
                }
                else
                    f[local2Global[i]][0] += loc.rhs(i)[0];
            }
        }
        // run over all leaf elements
        for (const auto& element : elements(gridView_))
        {
            unsigned int numFaces = element.subEntities(1);

            // get local to global id map
            for (unsigned int k = 0; k < numFaces; k++)
            {
                int alpha = faceMapper_.subIndex(element, k, 1);
                local2Global[k] = alpha;
            }
            loc.completeRHS(element, local2Global, f);
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
            }
            // insert dirichlet and processor boundary conditions
            else if (essential[i.index()][0]!=BoundaryConditions::neumann)
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
    const GridView gridView_;
    FaceMapper faceMapper_;
    unsigned int size_;
    RepresentationType A_;
};

} // end namespace Dumux

#endif
