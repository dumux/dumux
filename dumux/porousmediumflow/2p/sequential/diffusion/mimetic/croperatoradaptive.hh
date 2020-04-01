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
#ifndef DUMUX_CROPERATOR2PADAPTIVE_HH
#define DUMUX_CROPERATOR2PADAPTIVE_HH

#include<iostream>
#include<vector>
#include<cassert>
#include<stdio.h>
#include<stdlib.h>

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
#include <dumux/common/intersectionmapper.hh>

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
 *
 */
template<class TypeTag>
class CROperatorAssemblerTwoPAdaptive
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum {dim=GridView::dimension};
    using Element = typename GridView::template Codim<0>::Entity;
    using IS = typename GridView::IndexSet;
    using BlockType = Dune::FieldMatrix<Scalar, 1, 1>;
    using MatrixType = Dune::BCRSMatrix<BlockType>;
    using MBlockType = typename MatrixType::block_type;
    using rowiterator = typename MatrixType::RowIterator;
    using coliterator = typename MatrixType::ColIterator;
    using BCBlockType = std::array<BoundaryConditions::Flags, 1>;     // componentwise boundary conditions
    using SatType = Dune::BlockVector< Dune::FieldVector<double, 1> >;
    using IntersectionMapper = Dumux::IntersectionMapper<GridView>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        pressureEqIdx = Indices::pressureEqIdx,
    };

    // return number of rows/columns
    int size () const
    {
        return intersectionMapper_.size();
    }

public:
    using RepresentationType = MatrixType;

    CROperatorAssemblerTwoPAdaptive (const GridView& gridview)
    : gridView_(gridview), is_(gridView_.indexSet()), intersectionMapper_(gridView_)
    {
        A_.setBuildMode(MatrixType::random);
    }

    void initialize()
    {
        adapt();
    }

    void adapt()
    {
        intersectionMapper_.update();
        A_.setSize(size(), size());
        updateMatrix();
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

    const IntersectionMapper& intersectionMapper()
    {
        return intersectionMapper_;
    }

    const IntersectionMapper& intersectionMapper() const
    {
        return intersectionMapper_;
    }

    /*!
     * \brief Assembles global stiffness matrix
     *
     * This method takes an object that can compute local stiffness matrices and
     * assembles the global linear system Au=f.
     *
     * \param loc the local assembler providing element stiffness and boundary conditions for all elements
     * \param[in,out] u solution, contains initial values on input, Dirichlet values are set. The type of boundary
     * condition for a node is inferred from the values returned by the local assembler. A node is of Neumann
     * type if all elements referring to that node report a Neumann boundary condition, it is set to Dirichlet
     * if a least one element reports a process or Dirichlet boundary condition. The difference between process
     * and Dirichlet is that process always denotes a homogeneous Dirichlet value.
     * \param f right hand side is filled by this method
     *
     * Note that the rows corresponding to nodes at the Dirichlet boundary are filled with trivial equations
     * of the form \f[1\cdot u_i = f_i \f] where \f$u_i\f$ and \f$f_i\f$ are both set to the Dirichlet value at the \f$i\f$th node.
     */
    template <class LocalStiffness, class Vector>
    void assemble (LocalStiffness& loc, Vector& u, Vector& f);

    void updateMatrix();

protected:
    const GridView gridView_;
    const IS& is_;
    IntersectionMapper intersectionMapper_;
    RepresentationType A_;
};

template<class TypeTag>
void CROperatorAssemblerTwoPAdaptive<TypeTag>::updateMatrix()
{
    // set size of all rows to zero
    for (int i = 0; i < size(); i++)
        A_.setrowsize(i, 0);

    // build needs a flag for all entities of all codims
    std::vector<bool> visited(size(), false);

    int numElem = gridView_.size(0);

    for (int elemIdx = 0; elemIdx < numElem; elemIdx++)
    {
        int numFaces = intersectionMapper_.size(elemIdx);
        for (int fIdx = 0; fIdx < numFaces; fIdx++)
        {
            int faceIdxGlobal = intersectionMapper_.subIndex(elemIdx, fIdx);
            if (!visited[faceIdxGlobal])
            {
                A_.incrementrowsize(faceIdxGlobal);
                visited[faceIdxGlobal] = true;
            }
            for (int k = 0; k < numFaces-1; k++) {
                A_.incrementrowsize(faceIdxGlobal);
            }

        }

    }

    A_.endrowsizes();

    // clear the flags for the next round, actually that is not necessary because addindex takes care of this
    for (int i = 0; i < size(); i++)
        visited[i] = false;

    for (int elemIdx = 0; elemIdx < numElem; elemIdx++)
    {
        int numFaces = intersectionMapper_.size(elemIdx);
        for (int fIdx = 0; fIdx < numFaces; fIdx++)
        {
            int faceIdxGlobalI = intersectionMapper_.subIndex(elemIdx, fIdx);
            if (!visited[faceIdxGlobalI])
            {
                A_.addindex(faceIdxGlobalI,faceIdxGlobalI);
                visited[faceIdxGlobalI] = true;
            }
            for (int k = 0; k < numFaces; k++)
                if (k != fIdx) {
                    int faceIdxGlobalJ = intersectionMapper_.subIndex(elemIdx, k);
                    A_.addindex(faceIdxGlobalI, faceIdxGlobalJ);
                }

        }

    }

    A_.endindices();
}

template<class TypeTag>
template <class LocalStiffness, class Vector>
void CROperatorAssemblerTwoPAdaptive<TypeTag>::assemble(LocalStiffness& loc, Vector& u, Vector& f)
{

    // check size
    if (u.N()!=A_.M() || f.N()!=A_.N())
        DUNE_THROW(Dune::MathError,"CROperatorAssemblerTwoPAdaptive::assemble(): size mismatch");
    // clear global stiffness matrix and right hand side
    A_ = 0;
    f = 0;

    // allocate flag vector to hold flags for essential boundary conditions
    std::vector<BCBlockType> essential(intersectionMapper_.size());
    for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++)
            essential[i][0] = BoundaryConditions::neumann;

    // local to global id mapping (do not ask vertex mapper repeatedly
    std::vector<int> local2Global(2*dim);

    // run over all leaf elements
    for (const auto& element : elements(gridView_))
    {
        // build local stiffness matrix for CR elements
        // inludes rhs and boundary condition information
        loc.assemble(element, 1); // assemble local stiffness matrix

        int eIdxGlobal = intersectionMapper_.index(element);

        unsigned int numFaces = intersectionMapper_.size(eIdxGlobal);
        local2Global.resize(numFaces);

        for (unsigned int i = 0; i < numFaces; i++)
        {
                int idx = intersectionMapper_.subIndex(element, i);
                local2Global[i] = idx;
        }

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
        int eIdxGlobal = intersectionMapper_.index(element);

        unsigned int numFaces = intersectionMapper_.size(eIdxGlobal);
        local2Global.resize(numFaces);

        for (unsigned int i = 0; i < numFaces; i++)
        {
                int idx = intersectionMapper_.subIndex(element, i);
                local2Global[i] = idx;
        }

        loc.completeRHS(element, local2Global, f);
    }

    // put in essential boundary conditions
    rowiterator endi=A_.end();
    for (rowiterator i=A_.begin(); i!=endi; ++i)
    {
        // muck up extra rows
        if ((int) i.index() >= (int) size())
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

        // insert Dirichlet ans processor boundary conditions
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
} // end namespace Dumux

#endif
