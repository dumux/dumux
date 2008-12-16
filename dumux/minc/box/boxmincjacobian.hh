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
//#include<dune/disc/functions/p1function.hh>
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
  template<class G, class RT, int m, class BoxFunction = LeafP1FunctionExtended<G, RT, m> > class BoxMincJacobian
    : public BoxJacobian<BoxMincJacobian<G,RT,m,BoxFunction>,G,RT, m,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxMincJacobian<G,RT,m,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,m>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,m>::MBlockType MBlockType;
 	typedef Dune::FVElementGeometry<G> FVElementGeometry;
// number of phases

 	enum {wPhase = 0, nPhase = 1};
 	enum {F = 0, M = 1};
 	enum {pWFIdx = 0, satNFIdx = 1, pWMIdx = 2, satNMIdx = 3};
 	enum {WF = 0, NF = 1, WM = 2, NM = 3};

  public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
	// m represents the number of equations
    enum {n = G::dimension};
    enum {nPhases = 2}; //number of phases
    // number of interacting continua
    enum {nCont = m/2};

    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<RT,n,n> FMatrix;
    typedef FieldVector<RT,n> FVector;
    //! Constructor
    BoxMincJacobian (MincProblem<G,RT, m>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid,
			      BoxFunction& sol,
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,m,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params),
      statNData(this->vertexMapper.size()),
      vNData(SIZE), oldVNData(SIZE)
    {
      this->analytic = false;
    }
/***********************************************************************************************************/
    ////  Calculating the perimeter of the control volume that will be used as Aj,j+1
    virtual FieldVector<RT, nCont> AreaP (const Entity& e, int face)
//    virtual double AreaP (const Entity& e, int face)
    {
    			  FieldVector<RT,nCont> x(0.0);
    			  FieldVector<RT,nCont> y(0.0);
    			  FieldVector<RT,nCont> L(0.0);
    			  double length_x =0;  	//the length of the subcontrol face on x direction
    			  double length_y =0;  	//the length of the subcontrol face on y direction
    	    	  double ExtLength =0;		//Exterior length of the subcontrol face
    	      for (int face =0; face < this->fvGeom.numVertices; face++){
    	    	  int it =this->fvGeom.subContVolFace[face].i;
    	    	  int jt =this->fvGeom.subContVolFace[face].j;
    	    	  if (it == 0)
    	    	        			{
    	    	        				FieldVector<DT,n> ht;
    	    	        				ht = this->fvGeom.subContVolFace[face].normal;
    	    	        				if (ht[0]<0) {ht[0]*=(-1);}
    	    	        				if (ht[1]<0) {ht[1]*=(-1);}
    	    	        				if (ht[0]==0){ht[0]=ht[1];}
    	    	        				if (ht[1]==0){ht[1]=ht[0];}
    	    	        				length_y = ht[0];
    	    	        			}
    	    	   else if (jt==0)   {
    	    		   					FieldVector<DT,n> ht;
    	    		       	    	    ht = this->fvGeom.subContVolFace[face].normal;
    	    		       	    	    if (ht[0]<0) {ht[0]*=(-1);}
    	    		       	    	    if (ht[1]<0) {ht[1]*=(-1);}
    	    		       	    	    if (ht[0]==0){ht[0]=ht[1];}
    	    		       	    	    if (ht[1]==0){ht[1]=ht[0];}
    	    		       	    	    length_x = ht[0];
    	    	        			}
    	    	   ExtLength = length_x+length_y;

    	      }
//    	      double f = 1.0 / nCont; // volume fraction
//    	      // Calculation of the distance between the nested volume elements delta
//    	      // delta in this case divides the surfaces equidistantly (delta is const)
    	      elData.delta = 0.0;
    	      if (length_x > length_y){
    	    	  elData.delta = 2*length_y / (2.0 * (nCont-1) + 1.0);
    	      }
    	      else {
    	    	  elData.delta = 2*length_x / (2.0 * (nCont-1) + 1.0);
    	      }

    	      double TotalLength = 4*ExtLength;

	    	  L[0]  = TotalLength;
	    	  x[0]  = 2*length_x;	//length of the unitary element (CV) on x direction
	    	  y[0]  = 2*length_y;	//length of the unitary element (CV) on y direction
    	      for (int nC=1; nC<nCont; nC++){
    	      	  x[nC] = x[nC-1] - 2 * elData.delta;
    	      	  y[nC] = y[nC-1] - 2 * elData.delta;
    	      	  L[nC] = 2* (x[nC]+y[nC]);
    	      }
    	      return L;

//    	      return TotalLength;
    }

/***********************************************************************************************************/
/*the harmonicMeanKMinc function computes the harmonic mean of the perameabilities of fractures and rock at
 * the same global coordinate " i". It is used afterwards in the computation of the interporosity flux */
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
/*the harmonicMeanK function computes the harmonic mean of the perameabilities between the two nodes of different
 *  permeabilities in the fracture domain */
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
    		int node, std::vector<VariableNodeData>& vNData)
    {
   	 VBlockType AccumulationTerm;
     MBlockType result;

   	 for (int nC = 0; nC < nCont; nC++){
   		result[wPhase][nC] = - vNData[node].density[wPhase][nC]* elData.porosityFracture * vNData[node].saturation[nPhase][nC];
   		result[nPhase][nC] = vNData[node].density[nPhase][nC]* elData.porosityFracture * vNData[node].saturation[nPhase][nC];
   	 }

     for (int nC = 0; nC < nCont ;nC++ )
     {
     //2 is the number of phases (wetting and non-wetting) and nC - the number of Continua
    	 int count = 2*nC;
    		 AccumulationTerm[count]=result[wPhase][nC];
    		 AccumulationTerm[count+1]=result[nPhase][nC];
     }

   	 return AccumulationTerm;
    };

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
    	if (old)
    		return computeM(e, sol, node, oldVNData);
    	else
    		return computeM(e, sol, node, vNData);
    }

//***********************************************************************************************//
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
		  const FMatrix K = harmonicMeanK(global_i, global_j);
		  const FieldVector<DT,n> normal(this->fvGeom.subContVolFace[face].normal);
		  VBlockType flux(0.0);
		  FieldVector<RT, n> pGrad(0.0);
		  FieldVector<RT, n> pWFGrad(0.0);
		  FieldVector<RT, n> pNFGrad(0.0);
		  FieldVector<RT, n> pWMGrad(0.0);
		  FieldVector<RT, n> pNMGrad(0.0);
		  for (int comp = 0; comp < 4; comp++)
		  {
			  FieldVector<RT, n> gravity = problem.gravity();
			  switch (comp)	{
		  case WF:
	          for (int k = 0; k < this->fvGeom.numVertices; k++) {
	        	  FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
//	        	  grad *= sol[k][pWFIdx];
	        	  grad *= vNData[k].p[wPhase][F];
	        	  pWFGrad += grad;
	          }
	          // adjust by gravity
	          gravity *= vNData[i].density[wPhase][0];
			  pWFGrad -= gravity;
			  break;
		  case NF:
			  for (int k = 0; k < this->fvGeom.numVertices; k++) {
			  	  FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
			  	  grad *= vNData[k].p[nPhase][F];
			  	  pNFGrad += grad;
			  }
			  // adjust by gravity
			  gravity *= vNData[i].density[nPhase][0];
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
		// assigns the fully upwind mobility
	    if (out_WF<0)
		flux[pWFIdx]=vNData[i].density[wPhase][F]*vNData[i].mobility[wPhase][F]*out_WF;
		else
		flux[pWFIdx]=vNData[j].density[wPhase][F]*vNData[j].mobility[wPhase][F]*out_WF;
		if (out_NF<0)
		flux[satNFIdx]=vNData[i].density[nPhase][F]*vNData[i].mobility[nPhase][F]*out_NF;
		else
		flux[satNFIdx]=vNData[j].density[nPhase][F]*vNData[j].mobility[nPhase][F]*out_NF;
		}

			return flux;
   };

//******************************************************************************************//
//******************************************************************************************//

    virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
     {
    	 // ASSUME problem.q already contains \rho.q
      VBlockType q(0);

////      std::ofstream file;
////      file.open("file.txt", std::ios::app);


      q = problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
      const FieldVector<DT,n> global_i=this->fvGeom.subContVol[node].global;
      FieldVector <RT, nCont> A(0.0);
      A = AreaP(e, node);
      const FMatrix Kharmonic = harmonicMeanKMinc (global_i);
      int i = node;
//      double K_used = Kharmonic[0][0];
////  Product of permeability with the nestedvolume area
      FieldVector<RT, nCont> KA_nestVol(0.0);
      for (int nC=0; nC<nCont; nC++)
      {
    	  KA_nestVol[nC] = Kharmonic[0][0] * A[nC];
      }
/***********************************************************************************/
 //// MINC Flux Matrix
      FieldMatrix <RT, n, nCont> IPFlux(0.);
      FieldMatrix <RT, n, nCont> PGrad(0.);
      FieldVector <RT, nCont> inward(0.);
      for (int nC = 0; nC < nCont-1; nC++)
      {
          PGrad[wPhase][nC] += vNData[i].p[wPhase][nC];
          PGrad[wPhase][nC] -= vNData[i].p[wPhase][nC+1];
          PGrad[nPhase][nC] += vNData[i].p[nPhase][nC];
          PGrad[nPhase][nC] -= vNData[i].p[nPhase][nC+1];

////       //if the normal of the face gives negative numbers it multiplies with -1
               inward[wPhase] = PGrad[wPhase][nC] * KA_nestVol[nC] / elData.delta;
                  	if (inward[wPhase] < 0){
               inward[wPhase]*=(-1);
                	}
               inward[nPhase] = PGrad[wPhase][nC] * KA_nestVol[nC] / elData.delta;
                	if (inward[nPhase] < 0){
               inward[nPhase]*=(-1);
                	}
//*******************************************************************//
////      the position of the node
//      for (int node = 0; node < numVertices; node++) {
//                          fvGeom.subContVol[node].local  = referenceElement.position(node, dim);
//                          subContVol[node].global = geometry.global(subContVol[node].local);
//                      }
//*******************************************************************//

 //// If the pressure in the fracture is bigger than the one in the matrix the flow occurs from fracture to matrix
      if (vNData[i].p[wPhase][nC]>vNData[i].p[wPhase][nC+1])
       	{
        	// wetting interporosity flux calculated with the wetting mobility of the Fracture;
    	  	IPFlux[wPhase][nC]=vNData[i].density[wPhase][nC]*vNData[i].mobility[wPhase][nC]*inward[wPhase];
        	q[2*nC]-=IPFlux[wPhase][nC]; //pWFIdx
        	q[2*(nC+1)]+=IPFlux[wPhase][nC];//pWMIdx
      	 }
      else
        	{
    	  	// wetting interporosity flux calculated with the wetting mobility of the Matrix;
    	  	IPFlux[wPhase][nC]=vNData[i].density[wPhase][nC+1]*vNData[i].mobility[wPhase][nC+1]*inward[wPhase];
    	  	q[2*nC]+=IPFlux[wPhase][nC];//pWFIdx
    	  	q[2*(nC+1)]-=IPFlux[wPhase][nC];//pWMIdx
        	}
      if (vNData[i].p[nPhase][nC]>vNData[i].p[nPhase][nC+1])
        	{
    	  	IPFlux[nPhase][nC]=vNData[i].density[nPhase][nC]*vNData[i].mobility[nPhase][nC]*inward[nPhase];
    	  	q[2*nC+1]-=IPFlux[nPhase][0];//satNFIdx
        	q[2*(nC+1)+1]+=IPFlux[nPhase][0];//satNMIdx
        	}
      else
        	{
    	  	IPFlux[nPhase][nC]=vNData[i].density[nPhase][nC+1]*vNData[i].mobility[nPhase][nC+1]*inward[nPhase];
    	  	q[2*nC+1]+=IPFlux[nPhase][nC];//satNFIdx
        	q[2*(nC+1)+1]-=IPFlux[nPhase][nC];//satNMIdx
        	}
      }
/************************************************************************************/
//      VBlockType q_previous(0);
//      q_previous = q;
//      for (int c=0; c<4; c++){
//    	  file << q[c]<<"\t\t";
//      		}
//      file << std::endl;
//      file.close();

      return q;
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
 		 elData.parametersFracture = problem.materialLawParametersFracture(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
 		 elData.parametersMatrix = problem.materialLawParametersMatrix(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
		 // ASSUMING element-wise constant permeability, evaluate K at the cell center
 		 elData.KFracture = problem.KFracture(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
 		 elData.KMatrix = problem.KMatrix(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

 		 elData.porosityFracture = problem.porosityFracture(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
 		 elData.porosityMatrix = problem.porosityMatrix(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

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
	  //*	(vNData)										 	*
	  //*													 	*
	  //*********************************************************


    // analog to EvalPrimaryData in MUFTE, uses members of varNData
	virtual void updateVariableData(const Entity& e, const VBlockType* sol,
			int i, std::vector<VariableNodeData>& vNData)
    {
   	 vNData.resize(this->fvGeom.numVertices);
   	 int size = vNData.size();

   	 for (int i = 0; i < size; i++) {

   		 for (int nEq=0; nEq<m; nEq++)
   		 {
   	         int a = nEq%4;
   	         int pos = nEq/2;
   	         if (a==1){
   	        	vNData[i].saturation[wPhase][pos]  = 1.0 - sol[i][nEq];		//S wPhase F
   	        	vNData[i].saturation[nPhase][pos]  = sol[i][nEq];			//S nPhase F
   	         }
   			 if (a==3){
   				vNData[i].saturation[wPhase][pos]= 1.0 - sol[i][nEq]; 		//S wPhase M
   				vNData[i].saturation[nPhase][pos]= sol[i][nEq];			//S nPhase M
   			 }
   		 }

   		// ASSUME element-wise constant parameters for the material law
//   		FieldVector<RT, 4> parameters = problem.materialLawParameters(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

   		FieldVector<RT, 4> parametersFracture = problem.materialLawParametersFracture(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
  		FieldVector<RT, 4> parametersMatrix = problem.materialLawParametersMatrix(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

         vNData[i].pCFracture = problem.materialLaw().pC(vNData[i].saturation[wPhase][F], parametersFracture);
         vNData[i].pCMatrix   = problem.materialLaw().pC(vNData[i].saturation[wPhase][M], parametersMatrix);

         // it assigns the values of the solution vector which has 4 components to the vNData Matrix (2x2) for 2 continua.
         for (int nEq = 0; nEq < m; nEq++)
              {
         int a = nEq%4;
         int pos = nEq/2;
         if (a==0) {
         vNData[i].p[wPhase][pos]   = sol[i][nEq]; 							//P wPhase F
         vNData[i].p[nPhase][pos]   = sol[i][nEq] + vNData[i].pCFracture; 	//P nPhase F
         }
         if (a==2) {
         vNData[i].p[wPhase][pos] = sol[i][nEq]; 							//P wPhase M
         vNData[i].p[nPhase][pos] = sol[i][nEq] + vNData[i].pCMatrix; 		//P nPhase M
         }
           	 }

         //Mobilities
         for (int nC = 0; nC < nCont ;nC++ )
         {
        	 vNData[i].mobility[wPhase][nC] = problem.materialLaw().mobW(vNData[i].saturation[wPhase][nC], parametersMatrix);
        	 vNData[i].mobility[nPhase][nC] = problem.materialLaw().mobN(vNData[i].saturation[nPhase][nC], parametersMatrix);
         }

         //Densities
         for (int nC = 0; nC < nCont ;nC++ )
                 {
        		 vNData[i].density[wPhase][nC]  = problem.materialLaw().wettingPhase.density();
                 }

         for (int nC = 0; nC < nCont ;nC++ )
        	     {
        	     vNData[i].density[nPhase][nC]  = problem.materialLaw().nonwettingPhase.density();
        	     }
   	 }

    }
//*******************************************************************
	virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false)
	{
		if (old)
			updateVariableData(e, sol, i, oldVNData);
		else
			updateVariableData(e, sol, i, vNData);
	}



	void updateVariableData(const Entity& e, const VBlockType* sol, bool old = false)
	{
		int size = this->fvGeom.numVertices;

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

    	RT pCFracture;
    	RT pCMatrix;

       FieldMatrix<DT, nCont, nCont> p;
       FieldMatrix<DT, nCont, nCont> saturation;
       FieldMatrix<DT, nCont, nCont> mobility;
       FieldMatrix<DT, nCont, nCont> density;
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

   	 double delta; //distance between nested volume faces (kept constant - so all faces are equally distanced)
   	 } elData;


    // parameters given in constructor
    MincProblem<G,RT, m>& problem;
    std::vector<StaticNodeData> statNData;
    std::vector<VariableNodeData> vNData;
    std::vector<VariableNodeData> oldVNData;
  };

  /** @} */
}
#endif
