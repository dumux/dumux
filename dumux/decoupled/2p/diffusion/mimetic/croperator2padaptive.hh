// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
#ifndef DUMUX_CROPERATOR2PADAPTIVE_HH
#define DUMUX_CROPERATOR2PADAPTIVE_HH

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
#include<dune/geometry/type.hh>
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
 */

namespace Dumux
{

//forward declaration
template<class GridView>
class IntersectionMapper;

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
template<class TypeTag>
class CROperatorAssemblerTwoPAdaptive
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum {dim=GridView::dimension};
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::IndexSet IS;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef Dune::FieldMatrix<Scalar,1,1> BlockType;
    typedef Dune::BCRSMatrix<BlockType> MatrixType;
    typedef typename MatrixType::block_type MBlockType;
    typedef typename MatrixType::RowIterator rowiterator;
    typedef typename MatrixType::ColIterator coliterator;
    typedef Dune::array<BoundaryConditions::Flags,1> BCBlockType;     // componentwise boundary conditions
    typedef Dune::BlockVector< Dune::FieldVector<double,1> > SatType;
    typedef Dumux::IntersectionMapper<GridView> IntersectionMapper;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        pressEqIdx = Indices::pressEqIdx,
    };

    // return number of rows/columns
    int size () const
    {
        return intersectionMapper_.size();
    }

public:
    typedef MatrixType RepresentationType;

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

    const IntersectionMapper& intersectionMapper()
    {
        return intersectionMapper_;
    }

    const IntersectionMapper& intersectionMapper() const
    {
        return intersectionMapper_;
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
    void assemble (LocalStiffness& loc, Vector& u, Vector& f);

    void updateMatrix();

protected:
    const GridView& gridView_;
    const IS& is_;
    IntersectionMapper intersectionMapper_;
    RepresentationType A_;
};

template<class TypeTag>
void CROperatorAssemblerTwoPAdaptive<TypeTag>::updateMatrix()
{
    // set size of all rows to zero
    for (unsigned int i = 0; i < size(); i++)
        A_.setrowsize(i, 0);

    // build needs a flag for all entities of all codims
    std::vector<bool> visited(size(), false);

    int numElem = gridView_.size(0);

    for (int elemIdx = 0; elemIdx < numElem; elemIdx++)
    {
        int numFaces = intersectionMapper_.size(elemIdx);
        for (int faceIdx = 0; faceIdx < numFaces; faceIdx++)
        {
            int faceIdxGlobal = intersectionMapper_.map(elemIdx, faceIdx);
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
        for (int faceIdx = 0; faceIdx < numFaces; faceIdx++)
        {
            int faceIdxGlobalI = intersectionMapper_.map(elemIdx, faceIdx);
            if (!visited[faceIdxGlobalI])
            {
                A_.addindex(faceIdxGlobalI,faceIdxGlobalI);
                visited[faceIdxGlobalI] = true;
            }
            for (int k = 0; k < numFaces; k++)
                if (k != faceIdx) {
                    int faceIdxGlobalJ = intersectionMapper_.map(elemIdx, k);
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
    ElementIterator eendit = gridView_.template end<0>();
    for (ElementIterator eIt = gridView_.template begin<0>(); eIt!=eendit; ++eIt)
    {
        // build local stiffness matrix for CR elements
        // inludes rhs and boundary condition information
        loc.assemble(*eIt, 1); // assemble local stiffness matrix

        int globalIdx = intersectionMapper_.map(*eIt);

        unsigned int numFaces = intersectionMapper_.size(globalIdx);
        local2Global.resize(numFaces);

        for (int i = 0; i < numFaces; i++)
        {
                int idx = intersectionMapper_.map(*eIt, i);
                local2Global[i] = idx;
        }

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
            if (loc.bc(i).isDirichlet(pressEqIdx))
            {
                essential[local2Global[i]][0] = BoundaryConditions::dirichlet;
                f[local2Global[i]][0] = loc.rhs(i)[0];
            }
            else
                f[local2Global[i]][0] += loc.rhs(i)[0];
        }

    }
    // run over all leaf elements
    for (ElementIterator eIt = gridView_.template begin<0>(); eIt!=eendit; ++eIt)
    {
        int globalIdx = intersectionMapper_.map(*eIt);

        unsigned int numFaces = intersectionMapper_.size(globalIdx);
        local2Global.resize(numFaces);

        for (int i = 0; i < numFaces; i++)
        {
                int idx = intersectionMapper_.map(*eIt, i);
                local2Global[i] = idx;
        }

        loc.completeRHS(*eIt, local2Global, f);
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

template<class GridView>
class IntersectionMapper
{
    // mapper: one data element in every entity
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    typedef typename GridView::Grid Grid;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> ElementMapper;



public:
    IntersectionMapper(const GridView& gridview)
    : gridView_(gridview), elementMapper_(gridView_), size_(gridView_.size(1)), intersectionMapGlobal_(gridView_.size(0)), intersectionMapLocal_(gridView_.size(0))
    {
        ElementIterator eIt = gridView_.template begin<0>();

        int faceIdx = 0;

        IntersectionIterator isItEnd = gridView_.iend(*eIt);
        for (IntersectionIterator isIt = gridView_.ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            int idxInInside = isIt->indexInInside();

            standardLocalIdxMap_[idxInInside] = faceIdx;

            faceIdx++;
        }
    }

    const ElementMapper& elementMapper() const
    {
        return elementMapper_;
    }

    int map(const Element& element) const
    {
        return elementMapper_.map(element);
    }

    int map(int elemIdx, int faceIdx)
    {
        return intersectionMapGlobal_[elemIdx][faceIdx];
    }

    int map(int elemIdx, int faceIdx) const
    {
        return (intersectionMapGlobal_[elemIdx].find(faceIdx))->second;//use find() for const function!
    }

    int map(const Element& element, int faceIdx)
    {
        return intersectionMapGlobal_[map(element)][faceIdx];
    }

    int map(const Element& element, int faceIdx) const
    {
        return intersectionMapGlobal_[map(element)].find(faceIdx)->second;//use find() for const function!
    }

    int maplocal(int elemIdx, int faceIdx)
    {
        return intersectionMapLocal_[elemIdx][faceIdx];
    }

    int maplocal(int elemIdx, int faceIdx) const
    {
        return (intersectionMapLocal_[elemIdx].find(faceIdx))->second;//use find() for const function!
    }


    int maplocal(const Element& element, int faceIdx)
    {
        return intersectionMapLocal_[map(element)][faceIdx];
    }

    int maplocal(const Element& element, int faceIdx) const
    {
        return (intersectionMapLocal_[map(element)].find(faceIdx))->second;//use find() for const function!
    }

    // return number intersections
    unsigned int size () const
    {
        return size_;
    }

    unsigned int size (int elemIdx) const
    {
        return intersectionMapLocal_[elemIdx].size();
    }

    unsigned int size (const Element& element) const
    {
        return intersectionMapLocal_[map(element)].size();
    }

    void update()
    {
        elementMapper_.update();

        intersectionMapGlobal_.clear();
        intersectionMapGlobal_.resize(elementMapper_.size());
        intersectionMapLocal_.clear();
        intersectionMapLocal_.resize(elementMapper_.size());

        ElementIterator eItEnd = gridView_.template end<0>();
        for (ElementIterator eIt = gridView_.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            int globalIdx = map(*eIt);

            int faceIdx = 0;
            // run through all intersections with neighbors
            IntersectionIterator isItEnd = gridView_.iend(*eIt);
            for (IntersectionIterator isIt = gridView_.ibegin(*eIt); isIt != isItEnd; ++isIt)
            {
                int indexInInside = isIt->indexInInside();
                intersectionMapLocal_[globalIdx][faceIdx] = indexInInside;

                faceIdx++;
            }
        }

        int globalIntersectionIdx = 0;
        for (ElementIterator eIt = gridView_.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            int globalIdx = map(*eIt);

            int faceIdx = 0;
            // run through all intersections with neighbors
            IntersectionIterator isItEnd = gridView_.iend(*eIt);
            for (IntersectionIterator isIt = gridView_.ibegin(*eIt); isIt != isItEnd; ++isIt)
            {
                if (isIt->neighbor())
                {
                    ElementPointer neighbor = isIt->outside();
                    int globalIdxNeighbor = map(*neighbor);

                    if (eIt->level() > neighbor->level() || (eIt->level() == neighbor->level() && globalIdx < globalIdxNeighbor))
                    {

                        int faceIdxNeighbor = 0;
                        if (size(globalIdxNeighbor) == 2 * dim)
                        {
                            faceIdxNeighbor = standardLocalIdxMap_[isIt->indexInOutside()];
                        }
                        else
                        {
                            IntersectionIterator isItNEnd = gridView_.iend(*neighbor);
                            for (IntersectionIterator isItN = gridView_.ibegin(*neighbor); isItN != isItNEnd; ++isItN)
                            {
                                if (isItN->neighbor())
                                {
                                    if (isItN->outside() == eIt)
                                    {
                                        break;
                                    }
                                }
                                faceIdxNeighbor++;
                            }
                        }

                        intersectionMapGlobal_[globalIdx][faceIdx] = globalIntersectionIdx;
                        intersectionMapGlobal_[globalIdxNeighbor][faceIdxNeighbor] = globalIntersectionIdx;
                        globalIntersectionIdx ++;
                    }
                }
                else
                {
                    intersectionMapGlobal_[globalIdx][faceIdx] = globalIntersectionIdx;
                    globalIntersectionIdx ++;
                }
                faceIdx++;
            }
        }
        size_ = globalIntersectionIdx;
    }

protected:
    const GridView& gridView_;
    ElementMapper elementMapper_;
    unsigned int size_;
    std::vector<std::map<int, int> > intersectionMapGlobal_;
    std::vector<std::map<int, int> > intersectionMapLocal_;
    std::map<int, int> standardLocalIdxMap_;
};

}

#endif
