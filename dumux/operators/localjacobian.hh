// $Id: localjacobian.hh 404 2006-10-04 08:47:06Z oliver $

#ifndef DUNE_LOCALJACOBIAN_HH
#define DUNE_LOCALJACOBIAN_HH

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
  class LocalJacobian : public LocalStiffness<Imp, G, RT, m>
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

	void computeElementData(const Entity& e) {
		return this->getImp().computeElementData(e);
	}

	virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false) 
	{
		return this->getImp().updateVariableData(e, sol, i, old);
	}

	void updateVariableData(const Entity& e, const VBlockType* sol, bool old = false)
	{
		return this->getImp().updateVariableData(e, sol, old);
	}
    
    template<class TypeTag>
    void assemble (const Entity& e, bool itIsThis, int k = 1)
    {
    	fvGeom.update(e);
		computeElementData(e);
		
    	int size = e.template count<n>();
	  
    	// set to Zero 
    	for (int i=0; i < size; i++) {
    		this->bctype[i].assign(BoundaryConditions::neumann);
    		this->b[i] = 0;
    		this->def[i] = 0;
    	}
  	  
    	setLocalSolution(e);
    	
    	updateStaticData(e, u);
    	bool old = true;
    	updateVariableData(e, uold, old);
    	updateVariableData(e, u);

		localDefect<TypeTag>(e, u);

		VBlockType bTemp[size];
		for (int i=0; i<size; i++)
		{
			for (int equationnumber = 0; equationnumber < m; equationnumber++)
				if (this->bctype[i][equationnumber]==BoundaryConditions::neumann) 
					bTemp[i][equationnumber] = this->def[i][equationnumber];
					else
						bTemp[i][equationnumber] = 0;
		}
      
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
					RT eps = std::max(fabs(1e-5*u[j][comp]), 1e-5);
					for (int i = 0; i < size; i++) {
						uPlusEps[i] = u[i];
						uMinusEps[i] = u[i];
					}
					uPlusEps[j][comp] += eps;
					uMinusEps[j][comp] -= eps;

					updateVariableData(e, uPlusEps, j);
					
					// calculate the defect without taking into account BCs 
					// ASSUMES that BCs do not depend on the solution
					bool withoutBC = false;
					localDefect<TypeTag>(e, uPlusEps, withoutBC);
					VBlockType defuPlusEps[size];
					for (int i = 0; i < size; i++) 
						defuPlusEps[i] = def[i];
	      
					updateVariableData(e, uMinusEps, j);
					localDefect<TypeTag>(e, uMinusEps, withoutBC);

					updateVariableData(e, u, j);
					
					RT oneByEps = 0.5/eps;
					for (int i = 0; i < size; i++) 
						for (int compi = 0; compi < m; compi++)
							this->A[i][j][compi][comp] = oneByEps*(defuPlusEps[i][compi] - def[i][compi]);
				}
		}
		
//		for (int i = 0; i < size; i++) 
//			for (int compi = 0; compi < m; compi++) {
//				for (int j = 0; j < size; j++) {
//					for (int compj = 0; compj < m; compj++)
//						std::cout << std::setw(9) << this->A[i][j][compi][compj] << ", "; 
//					std::cout << "\t";
//				}
//				std::cout << std::endl;
//			}
			
       
		for (int i=0; i<size; i++) 
		{
			for (int equationnumber = 0; equationnumber < m; equationnumber++)
				if (this->bctype[i][equationnumber]==BoundaryConditions::neumann) {
					this->b[i][equationnumber] = bTemp[i][equationnumber];
			}
		}
      
		return;
    }
    
    template<class TypeTag>
    void localDefect (const Entity& e, const VBlockType* sol, bool withBC = true)
    {
      this->getImp().template localDefect<TypeTag>(e, sol, withBC);
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

    virtual void updateStaticData (const Entity& e, VBlockType* sol)
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
    VBlockType uold[SIZE];
    bool analytic;
  };


  /** @} */

}
#endif
