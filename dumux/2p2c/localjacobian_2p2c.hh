// $Id: localjacobian.hh 404 2006-10-04 08:47:06Z oliver $

#ifndef DUNE_2P2CLOCALJACOBIAN_HH
#define DUNE_2P2CLOCALJACOBIAN_HH

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
#include"dumux/operators/localstiffnessextended.hh"
#include<dune/disc/operators/boundaryconditions.hh>
#include"dumux/fvgeometry/fvelementgeometry.hh"

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
   * @brief base class for assembling local jacobian matrices
   *
   */


  /*! @brief Base class for local assemblers

  This class serves as a base class for local assemblers. It provides
  space and access to the local jacobian matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:
  - Imp  The implementation of the Interface with Barton-Nackman
  - G    A grid type
  - RT   The field type used in the elements of the jacobian matrix
  - m    number of degrees of freedom per node (system size)
   */
  template<class Imp, class G, class RT, int m>
  class LocalJacobian2p2c : public LocalStiffness<Imp, G, RT, m>
  {
    // grid types
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    enum {n=G::dimension};

  public:
    // types for matrics, vectors and boundary conditions
    typedef LocalStiffness<Imp, G, RT, m> LocalStiffness;
    typedef FieldMatrix<RT,m,m> MBlockType;                      // one entry in the stiffness matrix
    typedef FieldVector<RT,m> VBlockType;                        // one entry in the global vectors
    typedef array<BoundaryConditions::Flags,m> BCBlockType; // componentwise boundary conditions
    typedef FVElementGeometry<G> FVElementGeometry;
    enum {SIZE=8};


    template<class TypeTag>
    void assemble (const Entity& e, bool itIsThis, int k = 1)
    {
      fvGeom.update(e);
		
      int size = e.template count<n>();
	  
      // set to Zero 
      for (int i=0; i < size; i++) {
      	this->bctype[i].assign(BoundaryConditions::neumann);
      	this->b[i] = 0;
      	this->def[i] = 0;
      }
  	  
      setLocalSolution(e);

      for (int i=0; i < size; i++)
      	primaryVarSwitch(e, u, i);

      updateStaticData(e, u);
      localDefect<TypeTag>(e, u);
      // assemble boundary conditions 
      assembleBC<TypeTag> (e); 
  	  
      // add to defect 
//       for (int i=0; i < size; i++) {
// 	this->def[i] += this->b[i];
// 	//std::cout << "i = " << ", b[i] = " << this->b[i] << std::endl;
//       }

      VBlockType bTemp[size];
      for (int i=0; i<size; i++) 
      	if (this->bctype[i][0]==BoundaryConditions::neumann) 
      		bTemp[i] = this->def[i];
      	else
      		bTemp[i] = 0;
      
// 	std::cout << "defu         = " << this->def[0] << ", "  << this->def[1] 
// 		  << ", "  << this->def[2] << ", "  << this->def[3] << std::endl;
      
      if (analytic) {
      	analyticJacobian<TypeTag>(e, u);
      }
      else {
      	VBlockType defu[size];
      	for (int i = 0; i < size; i++) 
      		defu[i] = def[i];    	
      	VBlockType uPlusEps[size];
      	VBlockType uMinusEps[size];
    	
	for (int j = 0; j < size; j++) 
	  for (int comp = 0; comp < m; comp++) 
	    {
	      RT eps = std::max(fabs(1e-3*u[j][comp]), 1e-3);
	      for (int i = 0; i < size; i++) {
	      	uPlusEps[i] = u[i];
	      	uMinusEps[i] = u[i];
	      }
	      uPlusEps[j][comp] += eps;
	      uMinusEps[j][comp] -= eps;
	      
	      localDefect<TypeTag>(e, uPlusEps);
	      VBlockType defuPlusEps[size];
	      for (int i = 0; i < size; i++) 
		defuPlusEps[i] = def[i];
	      
	      localDefect<TypeTag>(e, uMinusEps);
	      
// 		std::cout << "defuPlusEps  = " << defuPlusEps[0] << ", "  << defuPlusEps[1] 
// 			  << ", "  << defuPlusEps[2] << ", "  << defuPlusEps[3] << std::endl;
// 		std::cout << "defuMinusEps = " << this->def[0] << ", "  << this->def[1] 
// 			  << ", "  << this->def[2] << ", "  << this->def[3] << std::endl;

	      RT oneByEps = 0.5/eps;
	      for (int i = 0; i < size; i++) 
		for (int compi = 0; compi < m; compi++)
		  this->A[i][j][compi][comp] = oneByEps*(defuPlusEps[i][compi] - def[i][compi]);
	    }
      }
      
      for (int i=0; i<size; i++) 
	if (this->bctype[i][0]==BoundaryConditions::neumann) {
	  this->b[i] = bTemp[i];
	}
      
      return;
    }
    
    template<class TypeTag>
    void localDefect (const Entity& e, const VBlockType* sol)
    {
      this->getImp().template localDefect<TypeTag>(e, sol);
    }
    
    void primaryVarSwitch (const Entity& e, VBlockType* sol, int k)
    {
      this->getImp().primaryVarSwitch(e, sol, k);
    }
    
    void setLocalSolution (const Entity& e)
    {
      this->getImp().setLocalSolution(e);
    }
    
    template<class TypeTag>
    void assembleBC (const Entity& e)
    {
      this->getImp().template assembleBC<TypeTag>(e);
    }
    
    template<class TypeTag>
    void analyticJacobian (const Entity& e, const VBlockType* sol)
    {
	  this->getImp().template analyticJacobian<TypeTag>(e, sol);
    }

    virtual void updateStaticData (const Entity& e, const VBlockType* sol)
    {
	  return;
    }

    virtual void clearVisited ()
    {
	  return;
    }

    //! access defect for each dof
    /*! Access defect for each degree of freedom. Elements are
      undefined without prior call to the assemble method.
    */
    const VBlockType& defect (int i) const
    {
      return def[i];
    }
    
    void setDt(double d)
    {
    	this->getImp().setDt(d);
    }
    
    virtual double getDt()
    {
    	return this->getImp().getDt();
    }
    
    FVElementGeometry fvGeom;
    VBlockType def[SIZE];
    VBlockType u[SIZE];
    bool analytic;
  };


  /** @} */

}
#endif
