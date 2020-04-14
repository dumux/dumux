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
 * \brief Base class for assembling local stiffness matrices.
 */
#ifndef DUMUX_LOCAL_STIFFNESS_HH
#define DUMUX_LOCAL_STIFFNESS_HH

#include <array>
#include<iostream>
#include<vector>
#include<set>
#include<map>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/geometry/type.hh>
#include<dune/grid/common/grid.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/matrix.hh>

#include<dumux/common/boundaryconditions.hh>
#include<dumux/porousmediumflow/sequential/properties.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Base class for local assemblers.
 *
 * This class serves as a base class for local assemblers. It provides space
 * and access to the local stiffness matrix. The actual assembling is done
 * in a derived class via the virtual assemble method.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam m number of degrees of freedom per node (system size)
 */
template<class TypeTag, int m>
class LocalStiffness
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    // grid types
    using Entity = typename GridView::template Codim<0>::Entity;
    enum {n=GridView::dimension};

public:
    // types for matrics, vectors and boundary conditions
    using MBlockType = Dune::FieldMatrix<Scalar, m, m>;                      // one entry in the stiffness matrix
    using VBlockType = Dune::FieldVector<Scalar, m>;                        // one entry in the global vectors
    using BCBlockType = std::array<BoundaryConditions::Flags, m>; // componentwise boundary conditions
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    virtual ~LocalStiffness ()
    {
    }

    /*!
     * \brief Assembles local stiffness matrix including boundary conditions for given element and order
     *
     * On exit the following things have been done:
     * - The stiffness matrix for the given entity and polynomial degree has been assembled and is
     * accessible with the mat() method.
     * - The boundary conditions have been evaluated and are accessible with the bc() method.
     *  The boundary conditions are either neumann, process or dirichlet. Neumann indicates
     * that the corresponding node (assuming a nodal basis) is at the Neumann boundary, process
     * indicates that the node is at a process boundary (arising from the parallel decomposition of the mesh).
     * Process boundaries are treated as homogeneous Dirichlet conditions, i.e. the corresponding value
     * in the right hand side is set to 0. Finally, Dirichlet indicates that the node is at the Dirichlet
     * boundary.
     * - The right hand side has been assembled. It contains either the value of the essential boundary
     * condition or the assembled source term and neumann boundary condition.
     * It is accessible via the rhs() method.
     *
     * \param  e a codim 0 entity reference
     * \param  k order of Lagrange basis (default is 1)
     */
    virtual void assemble (const Entity& e, int k=1) = 0;

    /*!
     * \brief assemble local stiffness matrix including boundary conditions for given element and order
     *
     * Unlike the method with only two arguments, this one additionally takes the local solution in order
     * to allow assembly of nonlinear operators.
     *
     * On exit the following things have been done:
     * - The stiffness matrix for the given entity and polynomial degree has been assembled and is
     * accessible with the mat() method.
     * - The boundary conditions have been evaluated and are accessible with the bc() method.
     *  The boundary conditions are either Neumann, process or Dirichlet. Neumann indicates
     * that the corresponding node (assuming a nodal basis) is at the Neumann boundary, process
     * indicates that the node is at a process boundary (arising from the parallel decomposition of the mesh).
     * Process boundaries are treated as homogeneous Dirichlet conditions, i.e. the corresponding value
     * in the right hand side is set to 0. Finally, Dirichlet indicates that the node is at the Dirichlet
     * boundary.
     * - The right hand side has been assembled. It contains either the value of the essential boundary
     * condition or the assembled source term and Neumann boundary condition.
     * It is accessible via the rhs() method.
     *
     * \param  e a codim 0 entity reference
     * \param localSolution The current solution on the entity, which is needed by nonlinear assemblers
     * \param  k order of Lagrange basis (default is 1)
     */
    virtual void assemble (const Entity& e, const Dune::BlockVector<VBlockType>& localSolution, int k=1) = 0;

    /*!
     * \brief assemble only boundary conditions for given element and order
     *
     * On exit the following things have been done:
     * - The boundary conditions have been evaluated and are accessible with the bc() method.
     *  The boundary conditions are either Neumann, process or Dirichlet. Neumann indicates
     * that the corresponding node (assuming a nodal basis) is at the Neumann boundary, process
     * indicates that the node is at a process boundary (arising from the parallel decomposition of the mesh).
     * Process boundaries are treated as homogeneous Dirichlet conditions, i.e. the corresponding value
     * in the right hand side is set to 0. Finally, Dirichlet indicates that the node is at the Dirichlet
     * boundary.
     * - The right hand side has been assembled as far as boundary conditions are concerned.
     * It contains either the value of the essential boundary
     * condition or the assembled Neumann boundary condition.
     * It is accessible via the rhs() method.
     *
     * \param  e a codim 0 entity reference
     * \param  k order of Lagrange basis (default is 1)
     */
    virtual void assembleBoundaryCondition (const Entity& e, int k=1) = 0;

    /*!
     * \brief Prints contents of local stiffness matrix
     *
     * \param s output stream
     * \param width the width
     * \param precision the precision
     */
    void print (std::ostream& s, int width, int precision)
    {
      // set the output format
      s.setf(std::ios_base::scientific, std::ios_base::floatfield);
      int oldprec = s.precision();
      s.precision(precision);

      for (int i=0; i<currentsize(); i++)
        {
          s << "FEM";    // start a new row
          s << " ";      // space in front of each entry
          s.width(4);    // set width for counter
          s << i;        // number of first entry in a line
          for (int j=0; j<currentsize(); j++)
            {
              s << " ";         // space in front of each entry
              s.width(width);   // set width for each entry anew
              s << mat(i,j);     // yeah, the number !
            }
          s << " ";         // space in front of each entry
          s.width(width);   // set width for each entry anew
          s << rhs(i);
          s << " ";         // space in front of each entry
          s.width(width);   // set width for each entry anew
          s << bc(i)[0];
          s << std::endl;// start a new line
        }


      // reset the output format
      s.precision(oldprec);
      s.setf(std::ios_base::fixed, std::ios_base::floatfield);
    }

    /*!
     * \brief Accesses the local stiffness matrix
     *
     * Access elements of the local stiffness matrix. Elements are
     * undefined without prior call to the assemble method.
     *
     * \param i index i
     * \param j index j
     */
    const MBlockType& mat (int i, int j) const
    {
      return A[i][j];
    }

    /*!
     * \brief Accesses right hand side
     *
     * Access elements of the right hand side vector. Elements are
     * undefined without prior call to the assemble method.
     *
     * \param i index i
     */
    const VBlockType& rhs (int i) const
    {
      return b[i];
    }

    /*!
     * \brief Accesses boundary condition for each dof
     *
     * Access boundary condition type for each degree of freedom. Elements are
     * undefined without prior call to the assemble method.
     *
     * \param i index i
     */
     const BoundaryTypes& bc (int i) const
     {
       return bctype[i];
     }

    /*!
     * \brief Sets the current size of the local stiffness matrix
     *
     * \param s size which is to be set
     */
    void setcurrentsize (int s)
    {
          A.setSize(s,s);
          b.resize(s);
          bctype.resize(s);
    }

    /*!
     * \brief Gets the current size of the local stiffness matrix
     */
    int currentsize ()
    {
            return A.N();
    }

protected:
    // assembled data
    Dune::Matrix<MBlockType> A;
    std::vector<VBlockType> b;
    std::vector<BoundaryTypes> bctype;
};

/*!
 * \ingroup SequentialTwoPModel
 * \brief Base class for linear local assemblers
 *
 * This class serves as a base class for linear local assemblers. It provides
 * space and access to the local stiffness matrix. The actual assembling is done
 * in a derived class via the virtual assemble method.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam m number of degrees of freedom per node (system size)
 */
template<class TypeTag, int m>
class LinearLocalStiffness : public LocalStiffness<TypeTag,m>
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    // grid types
    using Entity = typename GridView::template Codim<0>::Entity;
    enum {n=GridView::dimension};

public:
    // types for matrics, vectors and boundary conditions
    using MBlockType = Dune::FieldMatrix<Scalar, m, m>;                      // one entry in the stiffness matrix
    using VBlockType = Dune::FieldVector<Scalar, m>;                        // one entry in the global vectors
    using BCBlockType = std::array<BoundaryConditions::Flags, m>;    // componentwise boundary conditions

    /*! Initialize local stiffness matrix */
    LinearLocalStiffness ()
    {}

    virtual ~LinearLocalStiffness ()
    {
    }

    /*!
     * \brief Assembles local stiffness matrix including boundary conditions for given element and order
     *
     * On exit the following things have been done:
     * - The stiffness matrix for the given entity and polynomial degree has been assembled and is
     * accessible with the mat() method.
     * - The boundary conditions have been evaluated and are accessible with the bc() method.
     * The boundary conditions are either neumann, process or dirichlet. Neumann indicates
     * that the corresponding node (assuming a nodal basis) is at the Neumann boundary, process
     * indicates that the node is at a process boundary (arising from the parallel decomposition of the mesh).
     * Process boundaries are treated as homogeneous Dirichlet conditions, i.e. the corresponding value
     * in the right hand side is set to 0. Finally, Dirichlet indicates that the node is at the Dirichlet
     * boundary.
     * - The right hand side has been assembled. It contains either the value of the essential boundary
     * condition or the assembled source term and neumann boundary condition.
     * It is accessible via the rhs() method.
     *
     * \param  e a codim 0 entity reference
     * \param k order of Lagrange basis (default is 1)
     */
    virtual void assemble (const Entity& e, int k=1) = 0;

    /*!
     * \brief Assembles local stiffness matrix including boundary conditions for given element and order
     *
     * Since this is a base class for linear assemblers, the local solution will be ignored.
     *
     * \param  e a codim 0 entity reference
     * \param localSolution The current solution on the entity, which is needed by nonlinear assemblers
     *  \param  k order of Lagrange basis (default is 1)
     */
    virtual void assemble (const Entity& e, const Dune::BlockVector<VBlockType>& localSolution, int k=1)
    {
        assemble(e,k);
    }
};

} // end namespace Dumux
#endif
