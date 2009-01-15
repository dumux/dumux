// $Id$

#ifndef DUNE_BOX1P2CJACOBIAN_HH
#define DUNE_BOX1P2CJACOBIAN_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include"dumux/functions/p1functionextended.hh"

#include"dumux/operators/boxjacobian.hh"
#include"dumux/1p2c/1p2cproblem.hh"
#include"dumux/io/vtkmultiwriter.hh"


/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for one-phase two-component flow equation
 * @author Karin Erbertseder
 */



namespace Dune
{
  /** @addtogroup DISC_Disc
   *
   * @{
   */
  /**
   * @brief compute local jacobian matrix for the boxfile for one-phase two-component flow equation
   *
   */


  //! Derived class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the one-phase two-component flow equation

	    div j = q; j = -K grad u; in Omega

		u = g on Gamma1; j*n = J on Gamma2.

	Uses box scheme with the Lagrange shape functions.
	It should work for all dimensions and element types.
	All the numbering is with respect to the reference element and the
	Lagrange shape functions

	Template parameters are:

	- Grid  a DUNE grid type
	- RT    type used for return values
	
	A class for computing the transport equation with the primary variables x and p.
	The transport equation is coupled with the continuity equation.
  */

  template<class G, class RT, class BoxFunction = LeafP1FunctionExtended<G, RT, 2> >
  class Box1P2CJacobian
    : public BoxJacobian<Box1P2CJacobian<G,RT,BoxFunction>,G,RT,2,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef Box1P2CJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::MBlockType MBlockType;
    typedef FVElementGeometry<G> FVElementGeometry;

 	enum {konti = 0, transport = 1};	// Solution vector index
	

  public:
    enum {dim=G::dimension};
    enum {m=2};					//number of components
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<RT,dim,dim> FMatrix;
    typedef FieldVector<RT,dim> FVector;

    //! Constructor
    Box1P2CJacobian (OnePTwoCProblem<G,RT>& params, bool levelBoundaryAsDirichlet_, const G& grid,
			      BoxFunction& sol, bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params), vNDat(SIZE), oldVNDat(SIZE)
    {
      this->analytic = false;
    }

	/** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
	 *  @param e entity
	 *  @param sol solution vector
	 *  @param node local node id
	 *  @return storage term
	 */
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol,
    		int node, std::vector<VariableNodeData>& varData)
    {
    	 
   	 VBlockType result;
 
   	 result[0] = 0; //continuity equation has no storage term.
   	 
   	 result[1] = vNDat[node].porosity * vNDat[node].molefraction; //storage term of the transport equation
 
   	 //std::cout << result << " " << node << std::endl;
   	 return result;
    };

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
    	if (old)
    		return computeM(e, sol, node, oldVNDat);
    	else
    		return computeM(e, sol, node, vNDat);
    }

    /** @brief compute diffusive/advective fluxes, loop over subcontrol volume faces
	 *  @param e entity
	 *  @param sol solution vector
	 *  @param face face id
	 *  @return flux term 
     */
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
   	 int i = this->fvGeom.subContVolFace[face].i;
 	 int j = this->fvGeom.subContVolFace[face].j;

 	 // normal vector, value of the area of the scvf
	 const FieldVector<RT,dim> normal(this->fvGeom.subContVolFace[face].normal);

  	 // get global coordinates of nodes i,j
	 const FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
	 const FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;
	 const FieldVector<DT,dim> local_i = this->fvGeom.subContVol[i].local;
	 const FieldVector<DT,dim> local_j = this->fvGeom.subContVol[j].local;
	 
	 VBlockType flux;
	 
	 //calculation of the flux of the continuity equation
	 FieldVector<RT,dim> gradP(0);
	 for (int k = 0; k < this->fvGeom.numVertices; k++) {
	     FieldVector<RT,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
	     grad *= sol[k][0];
	     gradP += grad;
	 }
	 		  
	 FieldVector<RT,dim> KGradP(0);
	 FMatrix Ki(0), Kj(0);
	 
	 // calculate harmonic mean of permeabilities of nodes i and j
	 Ki = this->problem.soil().K(global_i,e,local_i);
	 Kj = this->problem.soil().K(global_j,e,local_j);
	 const FMatrix K = harmonicMeanK(Ki, Kj); //harmonic mean of the permeability
	 
	 K.umv(gradP, KGradP);					 //KGradP = KGradP + gradP
	     	
	 flux[0] = ( KGradP* this->fvGeom.subContVolFace[face].normal) / vNDat[i].viscosity;  
	 
	 //calculation of the flux of the transport equation
	 FieldVector<RT,dim> gradX(0);		 
	 for (int k=0; k<this->fvGeom.numVertices; k++) {
	     FieldVector<RT,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
	     grad *= sol[k][1];
	     gradX += grad;
	 }
	     	
	 const RT T = arithmeticMeanT(i, j, e);	//arithmetic mean of the tortuosity T
	     	
	 const RT P = arithmeticMeanP(i, j, e);	//arithmetic mean of the porosity P
	     	
	 RT outward = KGradP * this->fvGeom.subContVolFace[face].normal;													//calcualte the flux using upwind
	 if (outward < 0)
	    flux[1] = (outward * vNDat[i].molefraction / vNDat[i].viscosity) + (T * P * vNDat[i].diffCoeff * (gradX * this->fvGeom.subContVolFace[face].normal));
	 else
	    flux[1] = (outward * vNDat[j].molefraction / vNDat[j].viscosity) + (T * P * vNDat[i].diffCoeff * (gradX * this->fvGeom.subContVolFace[face].normal)); 
	     	
	 return flux; 
	 }	
 

      /** @brief integrate sources / sinks
       *  @param e entity
       *  @param sol solution vector
       *  @param node local node id
       *  @return source/sink term
       */
      virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
          {
              // ASSUME problem.q already contains \rho.q
              return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
          }


    // harmonic mean of the permeability computed directly
      virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
      {
      	double eps = 1e-20;

      	for (int kx=0; kx<dim; kx++){			
      		for (int ky=0; ky<dim; ky++){
      			if (Ki[kx][ky] != Kj[kx][ky])
      			{
      				Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
      			}
      		}
      	}
     	 return Ki;
      }

      //arithmetic mean of the tortuosity is computed directly
      virtual RT arithmeticMeanT (const int i, const int j, const Entity& e) const 
      {
      	
      	RT Ti, Tj;
      	
      	Ti = this->problem.soil().tortuosity(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
      	Tj = this->problem.soil().tortuosity(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
      	
      	if (Ti != Tj)
      	    {
      	     Ti = (Ti + Tj)/2;
      	     return Ti;
      	    }
      	else{
      		return Ti;
      	}
      }
       
      //arithmetic mean of the porosity is computed directly
      virtual RT arithmeticMeanP (const int i, const int j, const Entity& e) const 
      {
          	
       RT Pi, Pj;
          	
       Pi = this->problem.soil().porosity(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
       Pj = this->problem.soil().porosity(this->fvGeom.subContVol[j].global, e, this->fvGeom.subContVol[j].local);
          	
       if (Pi != Pj)
          {
          Pi = (Pi + Pj)/2;
          return Pi;
          }
          else{
          return Pi;
          }
       }
 

    virtual void clearVisited ()
    {
    	return;
   	}

    // updates old phase state after each time step
    virtual void updatePhaseState ()
    {
       return;
    }

    virtual void resetPhaseState ()
    {
       return;
    }

	  //*********************************************************
	  //*														*
	  //*	Calculation of Data at Elements (elData) 			*
	  //*						 								*
	  //*														*
	  //*********************************************************

    virtual void computeElementData (const Entity& e)
    {
//		 // ASSUMING element-wise constant permeability, evaluate K at the cell center
// 		 elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
//
//		 // ASSUMING element-wise constant porosity
// 		 elData.porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
   	 return;
    }


	  //*********************************************************
	  //*														*
	  //*	Calculation of Data at Nodes that has to be			*
	  //*	determined only once	(sNDat)						*
	  //*														*
	  //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticDataVS (const Entity& e, VBlockType* sol)
    {
   	 return;
    }

    // for initialization of the Static Data (sets porosity)
    virtual void updateStaticData (const Entity& e, VBlockType* sol)
    {
   	 return;
    }


	  //*********************************************************
	  //*														*
	  //*	Calculation of variable Data at Nodes				*
	  //*	(vNDat)												*
	  //*														*
	  //*********************************************************


    struct VariableNodeData
    {
    	RT porosity;
    	RT molefraction;
    	RT viscosity;
    	RT tortuosity;
    	RT diffCoeff;
    	FMatrix K; 
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
	virtual void updateVariableData(const Entity& e, const VBlockType* sol,
			int i, std::vector<VariableNodeData>& varData)
    {
		varData[i].porosity = problem.soil().porosity(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
   		varData[i].viscosity = problem.phase().viscosity();
   		varData[i].tortuosity = problem.soil().tortuosity(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
	   	varData[i].diffCoeff = problem.phase().diffCoeff();
	   	varData[i].molefraction = sol[i][1];
    }

 
 	virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false) 
 		{
 			
 			if (old)
 				updateVariableData(e, sol, i, oldVNDat);
 			else 
 				updateVariableData(e, sol, i, vNDat);
 		}
    
	void updateVariableData (const Entity& e, const VBlockType* sol, bool old = false)
	{
		int size = this->fvGeom.numVertices;
		vNDat.resize(size);
		oldVNDat.resize(size);
		
		for (int i = 0; i < size; i++) 
			updateVariableData(e, sol, i, old);
	}


	struct StaticNodeData
    {
   	 bool visited;
    };




    // parameters given in constructor
    OnePTwoCProblem<G,RT>& problem;
    std::vector<StaticNodeData> sNDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;


  };

}
#endif
