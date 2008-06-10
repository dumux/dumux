// $Id: localstiffness.hh 468 2008-03-27 10:35:12Z sander $

#ifndef DUNE_LOCALSTIFFNESS_HH
#define DUNE_LOCALSTIFFNESS_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/fixedarray.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/grid.hh>
#include<dune/istl/operators.hh>
#include<dune/disc/operators/boundaryconditions.hh>

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
   * @brief base class for assembling local stiffness matrices
   *
   */


  /*! @brief Base class for local assemblers

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:
  - Imp  The implementation of the Interface with Barton-Nackman
  - G    A grid type
  - RT   The field type used in the elements of the stiffness matrix
  - m    number of degrees of freedom per node (system size)
   */
  template<class Imp, class G, class RT, int m>
  class LocalStiffness 
  {
	// grid types
	typedef typename G::ctype DT;
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	enum {n=G::dimension};

  protected:
	enum {SIZE=27};

  public:
	// types for matrics, vectors and boundary conditions
	typedef FieldMatrix<RT,m,m> MBlockType;                      // one entry in the stiffness matrix
	typedef FieldVector<RT,m> VBlockType;                        // one entry in the global vectors
        typedef array<BoundaryConditions::Flags,m> BCBlockType; // componentwise boundary conditions

	/*! initialize local stiffness matrix */
	LocalStiffness ()
	{
	  currentsize_ = 0;
	}

	virtual ~LocalStiffness () 
	{
	}

	//! assemble local stiffness matrix including boundary conditions for given element and order
	/*! On exit the following things have been done:
	  - The stiffness matrix for the given entity and polynomial degree has been assembled and is
        accessible with the mat() method.
	  - The boundary conditions have been evaluated and are accessible with the bc() method. 
        The boundary conditions are either neumann, process or dirichlet. Neumann indicates
        that the corresponding node (assuming a nodal basis) is at the Neumann boundary, process
        indicates that the node is at a process boundary (arising from the parallel decomposition of the mesh).
        Process boundaries are treated as homogeneous Dirichlet conditions, i.e. the corresponding value
        in the right hand side is set to 0. Finally, Dirichlet indicates that the node is at the Dirichlet
        boundary.  
	  - The right hand side has been assembled. It contains either the value of the essential boundary
        condition or the assembled source term and neumann boundary condition. 
		It is accessible via the rhs() method.

	  @param[in]  e    a codim 0 entity reference
	  @param[in]  k    order of Lagrange basis (default is 1)
	 */
        template<class TypeTag>
	void assemble (const Entity& e, int k=1){
	  getImp().template assemble<TypeTag>(e,k);
	}

      /** \brief assemble local stiffness matrix including boundary conditions for given element and order
          
      Unlike the method with only two arguments, this one additionally takes the local solution in order
      to allow assembly of nonlinear operators.

      On exit the following things have been done:
	  - The stiffness matrix for the given entity and polynomial degree has been assembled and is
        accessible with the mat() method.
	  - The boundary conditions have been evaluated and are accessible with the bc() method. 
        The boundary conditions are either neumann, process or dirichlet. Neumann indicates
        that the corresponding node (assuming a nodal basis) is at the Neumann boundary, process
        indicates that the node is at a process boundary (arising from the parallel decomposition of the mesh).
        Process boundaries are treated as homogeneous Dirichlet conditions, i.e. the corresponding value
        in the right hand side is set to 0. Finally, Dirichlet indicates that the node is at the Dirichlet
        boundary.  
	  - The right hand side has been assembled. It contains either the value of the essential boundary
        condition or the assembled source term and neumann boundary condition. 
		It is accessible via the rhs() method.

	  @param[in]  e    a codim 0 entity reference
          @param[in] localSolution The current solution on the entity, which is needed by nonlinear assemblers
	  @param[in]  k    order of Lagrange basis (default is 1)
	 */
        template<class TypeTag>
	void assemble (const Entity& e, const BlockVector<VBlockType>& localSolution, int k=1){
            getImp().template assemble<TypeTag>(e,localSolution, k);
	}
    

	//! assemble only boundary conditions for given element and order
	/*! On exit the following things have been done:
	  - The boundary conditions have been evaluated and are accessible with the bc() method. 
        The boundary conditions are either neumann, process or dirichlet. Neumann indicates
        that the corresponding node (assuming a nodal basis) is at the Neumann boundary, process
        indicates that the node is at a process boundary (arising from the parallel decomposition of the mesh).
        Process boundaries are treated as homogeneous Dirichlet conditions, i.e. the corresponding value
        in the right hand side is set to 0. Finally, Dirichlet indicates that the node is at the Dirichlet
        boundary.  
	  - The right hand side has been assembled as far as boundary conditions are concerned. 
        It contains either the value of the essential boundary
        condition or the assembled neumann boundary condition. 
		It is accessible via the rhs() method.

	  @param[in]  e    a codim 0 entity reference
	  @param[in]  k    order of Lagrange basis (default is 1)
	 */
    template<class TypeTag>
    void assembleBoundaryCondition (const Entity& e, int k=1)
    {
      getImp().template assembleBoundaryCondition<TypeTag>(e,k);
    }
    

	//! print contents of local stiffness matrix
	void print (std::ostream& s, int width, int precision)
	{
	  // set the output format
	  s.setf(std::ios_base::scientific, std::ios_base::floatfield);
	  int oldprec = s.precision();
	  s.precision(precision);

	  for (int i=0; i<currentsize_; i++)
		{
		  s << "FEM";    // start a new row
		  s << " ";      // space in front of each entry
		  s.width(4);    // set width for counter
		  s << i;        // number of first entry in a line
		  for (int j=0; j<currentsize_; j++)
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

	//! access local stiffness matrix
	/*! Access elements of the local stiffness matrix. Elements are
	  undefined without prior call to the assemble method.
	 */
	const MBlockType& mat (int i, int j) const
	{
	  return A[i][j];
	}


	//! access right hand side
	/*! Access elements of the right hand side vector. Elements are
	  undefined without prior call to the assemble method.
	 */
	const VBlockType& rhs (int i) const
	{
	  return b[i];
	}

	//! access boundary condition for each dof
	/*! Access boundary condition type for each degree of freedom. Elements are
	  undefined without prior call to the assemble method.
	 */
 	const BCBlockType& bc (int i) const
 	{
 	  return bctype[i];
 	}
 	
 	const VBlockType& dirichletIdx (int i) const
 	{
 		return dirichletIndex[i];
 	}

	//! set the current size of the local stiffness matrix
	void setcurrentsize (int s)
	{
	  if (s>=SIZE)
 		DUNE_THROW(MathError,"LocalStiffness: increase SIZE");		
	  currentsize_ = s;
	}

	//! get the current size of the local stiffness matrix
	int currentsize ()
	{
	  return currentsize_;
	}

  protected:
	// assembled data
	int currentsize_;
	MBlockType A[SIZE][SIZE];
	VBlockType b[SIZE];
	BCBlockType bctype[SIZE];
	VBlockType dirichletIndex[SIZE];
	
    Imp& getImp()
    {
      return static_cast<Imp&>(*this);
    }
    const Imp& getImp() const
    {
      return static_cast<const Imp&>(*this);
    }
    
  };


  /** @} */

}
#endif
