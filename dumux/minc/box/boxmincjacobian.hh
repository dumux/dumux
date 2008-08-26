#ifndef DUNE_BOXMINCJACOBIAN_HH
#define DUNE_BOXMINCJACOBIAN_HH

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
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/functions/p1function.hh>
#include"dumux/operators/boxjacobian.hh"
#include"dumux/minc/mincproblem.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 * @modified Alex Tatomir
 */



namespace Dune
{
  /** @addtogroup DISC_Disc
   *
   * @{
   */
  /**
   * @brief compute local jacobian matrix for conforming finite elements for diffusion equation
   *
   */


  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the 
	diffusion equation

	    div j = q; j = -K grad u; in Omega

		u = g on Gamma1; j*n = J on Gamma2.

	Uses conforming finite elements with the Lagrange shape functions.
	It should work for all dimensions and element types.
	All the numbering is with respect to the reference element and the
	Lagrange shape functions

	Template parameters are:

	- Grid  a DUNE grid type
	- RT    type used for return values 
  */
  template<class G, class RT, class BoxFunction = LeafP1Function<G, RT, 4> > class BoxMincJacobian 
    : public BoxJacobian<BoxMincJacobian<G,RT,BoxFunction>,G,RT,4,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxMincJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,4>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,4>::MBlockType MBlockType;
 	typedef FVElementGeometry<G> FVElementGeometry;
 	enum {pWFIdx = 0, satNFIdx = 1, pWMIdx = 2, satNMIdx = 3};
// 	enum {pWFIdx = 0, satNFIdx = 1};
 	enum {WF = 0, NF = 1, WM = 2, NM = 3};
	
  public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
	// m has to be modified minclensproblem, mincproblem, mincmodel  
    enum {n=G::dimension};
    enum {m=4};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    struct VariableNodeData;
    
    typedef FieldMatrix<RT,n,n> FMatrix;
    typedef FieldVector<RT,n> FVector;
    //! Constructor
    BoxMincJacobian (MincProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,4,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      statNData(this->vertexMapper.size()),
      vNData(SIZE), oldVNData(SIZE)
    {
      this->analytic = false;
    }
/***********************************************************************************************************/
    virtual FMatrix harmonicMeanKMinc (const FVector global_i) const
     {
    	 double eps = 1e-20;

    	 FieldMatrix<RT,n,n> KF, KM;
    	 
    	 KF = this->problem.K1Fracture(global_i); 
    	 KM = this->problem.K1Matrix(global_i);

    	 for (int kx=0; kx<n; kx++){
       	 for (int ky=0; ky<n; ky++){
       		 if (KF[kx][ky] != KM[kx][ky])
       			 {
       				 KF[kx][ky] = 2 / (1/(KF[kx][ky]+eps) + (1/(KM[kx][ky]+eps)));
       			 }
       	 }
 		 }
    	 return KF;
     }
    
    virtual FMatrix harmonicMeanK (const FVector global_i, const FVector global_j) const
     {
    	 double eps = 1e-20;

    	 FieldMatrix<RT,n,n> Ki, Kj;
    	 
    	 Ki = this->problem.K(global_i); 
    	 Kj = this->problem.K(global_j);

    	 for (int kx=0; kx<n; kx++){
       	 for (int ky=0; ky<n; ky++){
       		 if (Ki[kx][ky] != Kj[kx][ky])
       			 {
       				 Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
       			 }
       	 }
 		 }
    	 return Ki;
     }
/*************************************************************************************************************/     

    virtual void clearVisited ()
    {
    	return;
    }

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, 
    		int node, std::vector<VariableNodeData>& varData)
    {
   	 VBlockType result; 
   	 double dummy = varData[node].density[pWFIdx];
   	 result[0] = -varData[node].density[pWFIdx]*elData.porosityFracture*sol[node][satNFIdx];
   	 result[1] = varData[node].density[satNFIdx]*elData.porosityFracture*sol[node][satNFIdx];
//   	 result[2] = sol[node][pWMIdx];
//   	 result[3] = sol[node][satNMIdx];
   	 result[2] = -varData[node].density[pWMIdx]*elData.porosityMatrix*sol[node][satNMIdx];
   	 result[3] = varData[node].density[satNMIdx]*elData.porosityMatrix*sol[node][satNMIdx];
   	 
   	 return result;
    };
    
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false) 
    {
    	if (old)
    		return computeM(e, sol, node, oldVNData);
    	else 
    		return computeM(e, sol, node, vNData);
    }
    
    
    
    virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
   	 // ASSUME problem.q already contains \rho.q
     VBlockType q(0);

     q = problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
     const FieldVector<DT,n> global_i=this->fvGeom.subContVol[node].global;
     const FMatrix Kharmonic = harmonicMeanKMinc (global_i);
     
     int i = node;
     FieldVector<DT,n> Kharmonicij(0.);
     Kharmonic.umv(this->fvGeom.subContVolFace[node].normal, Kharmonicij);
 
     /***********************************************************************************/     
//// MINC Flux Matrix
     VBlockType IPFlux;
     FieldVector<RT,n> MincPWGrad(0);
     FieldVector<RT,n> MincPNGrad(0);
     MincPWGrad += vNData[i].pWFracture;
     MincPWGrad -= vNData[i].pWMatrix;
     MincPNGrad += vNData[i].pNFracture;
     MincPNGrad -= vNData[i].pNMatrix;
     RT inwardW = MincPWGrad*Kharmonicij;
     RT inwardN = MincPNGrad*Kharmonicij;
//// If the pressure in the fracture is bigger than the one in the matrix the flow occurs from fracture to matrix
       	if (inwardW>0)
      	{IPFlux[2] = vNData[i].density[pWFIdx]*vNData[i].mobilityMatrix[pWFIdx]*inwardW;
     	q[0]-=IPFlux[2];
      	q[2]+=IPFlux[2];
     	}
     else
       	{IPFlux[2] = vNData[i].density[pWMIdx]*vNData[i].mobilityFracture[pWMIdx]*inwardW;
       	q[0]-=IPFlux[2];
       	q[2]+=IPFlux[2];
       	}
     if (inwardN>0)
       	{IPFlux[3] = vNData[i].density[satNFIdx]*vNData[i].mobilityMatrix[satNFIdx]*inwardN;
       	q[1]-=IPFlux[3];
       	q[3]+=IPFlux[3];
       	}
     else
       	{IPFlux[3] = vNData[i].density[satNMIdx]*vNData[i].mobilityFracture[satNMIdx]*inwardN;
       	q[1]-=IPFlux[3];
       	q[3]+=IPFlux[3];
       	}
          
     /************************************************************************************/     
     return q; 
    };
    
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
    	int i = this->fvGeom.subContVolFace[face].i;
    	int j = this->fvGeom.subContVolFace[face].j;
     	 
		  // permeability in edge direction 
		  FieldVector<DT,n> KFractureij(0.);
		  elData.KFracture.umv(this->fvGeom.subContVolFace[face].normal, KFractureij);
		  FieldVector<DT,n> h; 
		  h = this->fvGeom.subContVolFace[face].normal;
		  // permeability of the Interacting Continua
		  FieldVector<DT,n> KMatrixij(0.);
		  elData.KMatrix.umv(this->fvGeom.subContVolFace[face].normal, KMatrixij);
		  const FieldVector<DT,n> global_i = this-> fvGeom.subContVol[i].global;
		  const FieldVector<DT,n> global_j = this-> fvGeom.subContVol[j].global;
		  const FMatrix K= harmonicMeanK(global_i, global_j);
		  const FieldVector<DT,n> normal(this->fvGeom.subContVolFace[face].normal);
		  VBlockType flux(0.0);
		  FieldVector<RT, n> pGrad(0.0);
		  FieldVector<RT, n> pWFGrad(0.0);
		  FieldVector<RT, n> pNFGrad(0.0);
		  FieldVector<RT, n> pWMGrad(0.0);
		  FieldVector<RT, n> pNMGrad(0.0);

		  for (int comp = 0; comp < m; comp++) 
		  {	  	  
			  FieldVector<RT, n> gravity = problem.gravity();
			  switch (comp)	{
		  case WF:
	          for (int k = 0; k < this->fvGeom.nNodes; k++) {
	        	  FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
	        	  grad *= sol[k][pWFIdx];
	        	  pWFGrad += grad;
	          }
	          // adjust by gravity			 
	          gravity *= vNData[i].density[pWFIdx];
			  pWFGrad -= gravity;
			  break;
		  case NF:
			  for (int k = 0; k < this->fvGeom.nNodes; k++) {
			  	  FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
			  	  grad *= vNData[k].pNFracture;
			  	  pNFGrad += grad;
			  }
			  // adjust by gravity			 
			  gravity *= vNData[i].density[satNFIdx];
			  pNFGrad -= gravity;
			  break;
		  case WM:
		  	  break;
		  case NM:
		  	  break;			  
			  
		  }	

			  FieldVector<RT,n>Kij(0.0);
			  K.umv(normal, Kij);
		  
			  RT out_WF = pWFGrad*Kij;
			  RT out_NF = pNFGrad*Kij;
			  
			  int up_WF, down_WF, up_NF, down_NF;
			  if (out_WF<0){
			  up_WF = i; down_WF=j;}
			  else {
			  up_WF = j; down_WF=i;}
			  if (out_NF<0){
			  up_NF = i; down_NF=j;}
			  else {
			  up_NF = j; down_NF=i;} 
		if (out_WF<0)
		flux[pWFIdx]=vNData[i].density[pWFIdx]*vNData[i].mobilityFracture[pWFIdx]*out_WF;
		else
		flux[pWFIdx]=vNData[j].density[pWFIdx]*vNData[j].mobilityFracture[pWFIdx]*out_WF;
		if (out_NF<0)
		flux[satNFIdx]=vNData[i].density[satNFIdx]*vNData[i].mobilityFracture[satNFIdx]*out_NF;
		else
		flux[satNFIdx]=vNData[j].density[satNFIdx]*vNData[j].mobilityFracture[satNFIdx]*out_NF;
		} 

			return flux;
   };
  
    
	  //*********************************************************
	  //*																*
	  //*	Calculation of Data at Elements			 					*
	  //*						 										*
	  //*															 	*
	  //*********************************************************

    virtual void computeElementData (const Entity& e)
    {
  		 // ASSUME element-wise constant parameters for the material law 
 		 elData.parametersFracture = problem.materialLawParametersFracture(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
 		 elData.parametersMatrix = problem.materialLawParametersMatrix(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
		 // ASSUMING element-wise constant permeability, evaluate K at the cell center 
 		 elData.KFracture = problem.KFracture(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);  
 		 elData.KMatrix = problem.KMatrix(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
 		 elData.porosityFracture = problem.porosityFracture(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
 		 elData.porosityMatrix = problem.porosityMatrix(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
 		 
    };

    
	  //*********************************************************
	  //*															*
	  //*	Calculation of Data at Nodes that has to be			 	*
	  //*	determined only once	(statNData)						*
	  //*															*
	  //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, const VBlockType* sol)
    {
	  return;
    }


	  //*********************************************************
	  //*														*
	  //*	Calculation of variable Data at Nodes			 	*
	  //*	(varNData)										 	*
	  //*													 	*
	  //*********************************************************
   

    // analog to EvalPrimaryData in MUFTE, uses members of varNData
	virtual void updateVariableData(const Entity& e, const VBlockType* sol, 
			int i, std::vector<VariableNodeData>& varData) 
    {
   	 varData.resize(this->fvGeom.nNodes);
   	 int size = vNData.size();

   	 for (int i = 0; i < size; i++) {
//   		this->def[i] = 0; // it initialize the deffect to 0
   		varData[i].saturationWFracture = 1.0 - sol[i][satNFIdx];
   		varData[i].saturationWMatrix = 1.0 - sol[i][satNMIdx];
//   		   		

   		// ASSUME element-wise constant parameters for the material law
//   		FieldVector<RT, 4> parameters = problem.materialLawParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
   	        	 
   		FieldVector<RT, 4> parametersFracture = problem.materialLawParametersFracture(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
  		FieldVector<RT, 4> parametersMatrix = problem.materialLawParametersMatrix(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
  	 
//         varNData[i].pC = problem.materialLaw().pC(varNData[i].saturationW, parameters);
         varData[i].pCFracture = problem.materialLaw().pC(varData[i].saturationWFracture, parametersFracture);
         varData[i].pCMatrix = problem.materialLaw().pC(varData[i].saturationWMatrix, parametersMatrix);
         varData[i].pWFracture = sol[i][pWFIdx];
//         varData[i].pWMatrix = sol[i][pWMIdx];
         varData[i].pNFracture = sol[i][pWFIdx] + varData[i].pCFracture;
         varData[i].pNMatrix = sol[i][pWMIdx] + varData[i].pCMatrix;
         
//         varNData[i].mobility[pWFIdx] = problem.materialLaw().mobW(varNData[i].saturationW, parameters);
         varData[i].mobilityFracture[pWFIdx] = problem.materialLaw().mobW(varData[i].saturationWFracture, parametersFracture);
         varData[i].mobilityMatrix[pWMIdx] = problem.materialLaw().mobW(varData[i].saturationWMatrix, parametersMatrix);
         
//         varNData[i].mobility[satNFIdx] = problem.materialLaw().mobN(sol[i][satNIdx], parameters);
         varData[i].mobilityFracture[satNFIdx] = problem.materialLaw().mobN(sol[i][satNFIdx], parametersFracture);
         varData[i].mobilityMatrix[satNMIdx] = problem.materialLaw().mobN(sol[i][satNMIdx], parametersMatrix);
                  
         varData[i].density[pWFIdx] = problem.materialLaw().wettingPhase.density();
         varData[i].density[satNFIdx] = problem.materialLaw().nonwettingPhase.density();
         varData[i].density[pWMIdx] = problem.materialLaw().wettingPhase.density();
         varData[i].density[satNMIdx] = problem.materialLaw().nonwettingPhase.density();
   	 }   

    }

	virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false) 
	{
		if (old)
			updateVariableData(e, sol, i, oldVNData);
		else 
			updateVariableData(e, sol, i, vNData);
	}

    
    
	void updateVariableData(const Entity& e, const VBlockType* sol, bool old = false)
	{
		int size = this->fvGeom.nNodes;
			
		for (int i = 0; i < size; i++) 
				updateVariableData(e, sol, i, old);
	}
    
    struct StaticNodeData 
    {
    	bool visited;
    	FieldMatrix<DT,n,n> K;
    	FieldMatrix<DT,n,n> K1Fracture;
    	FieldMatrix<DT,n,n> K1Matrix;
    };
    
       // the members of the struct are defined here
    struct VariableNodeData  
    {
//       RT saturationW;
       RT saturationWFracture;
       RT saturationWMatrix;
//       RT pC;
       RT pCFracture;
       RT pCMatrix;
//       RT pN;
       RT pNFracture;
       RT pNMatrix;
       RT pWFracture;
       RT pWMatrix;
//       VBlockType mobility;  //Vector with the number of phases
       VBlockType mobilityFracture;  //Vector with the number of phases
       VBlockType mobilityMatrix;  //Vector with the number of phases
       VBlockType density;
    };
        
    struct ElementData {
   	 RT cellVolume;
//     RT porosity;
     RT porosityFracture;
     RT porosityMatrix;
   	 RT gravity;
//   	 FieldVector<RT, 4> parameters;
   	 FieldVector<RT, 4> parametersMatrix;
   	 FieldVector<RT, 4> parametersFracture;

   	 FieldMatrix<DT,n,n> KFracture;
   	 FieldMatrix<DT,n,n> KMatrix;
   	 } elData;
   	 
    
    // parameters given in constructor
    MincProblem<G,RT>& problem;
    std::vector<StaticNodeData> statNData;
    std::vector<VariableNodeData> vNData;
    std::vector<VariableNodeData> oldVNData;
  };

  /** @} */
}
#endif
