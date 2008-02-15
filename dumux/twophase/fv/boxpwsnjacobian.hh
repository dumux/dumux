#ifndef DUNE_BOXPWSNJACOBIAN_HH
#define DUNE_BOXPWSNJACOBIAN_HH

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
#include "dumux/operators/localjacobian.hh"
#include "dumux/twophase/twophaseproblem.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
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

int edgeToFaceHex[8][8] = {
		{9, 2, 0, 9, 0, 9, 9, 9}, 
		{4, 9, 9, 4, 9, 2, 9, 9}, 
		{4, 9, 9, 3, 9, 9, 0, 9}, 
		{9, 1, 4, 9, 9, 9, 9, 1}, 
		{2, 9, 9, 9, 9, 5, 0, 9}, 
		{9, 1, 9, 9, 2, 9, 9, 5}, 
		{9, 9, 3, 9, 5, 9, 9, 3}, 
		{9, 9, 9, 3, 9, 1, 5, 9}
};

// ASSUME planar element faces
template <class CoordVec> 
double quadrilateralArea(const CoordVec& node0, const CoordVec& node1, 
							const CoordVec& node2, const CoordVec& node3)
{
	CoordVec distVec(node0);
	distVec -= node1;
	double a = distVec.two_norm();
	
	distVec = node1; 
	distVec -= node2;
	double b = distVec.two_norm();
	
	distVec = node2; 
	distVec -= node3;
	double c = distVec.two_norm();
	
	distVec = node3; 
	distVec -= node0;
	double d = distVec.two_norm();
	
	double s = 0.5*(a + b+ c + d);
	
	return sqrt((s - a)*(s - b)*(s - c)*(s - d));
}

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
  template<class G, class RT, class BoxFunction = LeafP1Function<G, RT, 2> >
  class BoxPwSnLocalJacobian 
    : public LocalJacobian<BoxPwSnLocalJacobian<G,RT>,G,RT,2>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef typename LocalJacobian<BoxPwSnLocalJacobian<G,RT>,G,RT,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<BoxPwSnLocalJacobian<G,RT>,G,RT,2>::MBlockType MBlockType;
 	typedef FVElementGeometry<G> FVElementGeometry;
  	enum {pWIdx = 0, satNIdx = 1};
	
  public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    
    //! Constructor
    BoxPwSnLocalJacobian (TwoPhaseProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : problem(params),levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_), 
    procBoundaryAsDirichlet(procBoundaryAsDirichlet_), 
      currentSolution(sol), oldSolution(grid), dt(1)
    {
      this->analytic = false;
    }
    

    template<class TypeTag>
    void localDefect (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol)
    {
      // extract some important parameters
      const Geometry& geometry = e.geometry();
      GeometryType gt = geometry.type();
      const typename LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
      	sfs=LagrangeShapeFunctions<DT,RT,n>::general(gt, 1);
      int size = sfs.size();
      this->setcurrentsize(size);
      
      // get cell volume
      DT cellVolume = geometry.volume();
      
      // assuming rectangular grids:
      cellVolume /= pow(2.0, n);
      
      // cell center in reference element
      const FieldVector<DT,n> 
		  cellLocal = ReferenceElements<DT,n>::general(gt).position(0,0);
		
	  // get global coordinate of cell center
	  const FieldVector<DT,n> cellGlobal = geometry.global(cellLocal);
      
	  // calculate secondary variables and set defect to zero 
      RT saturationW[size];
      RT pC[size];
      RT pN[size];
      VBlockType mobility[size];
      VBlockType density;
      // ASSUME element-wise constant parameters for the material law 
      FieldVector<RT, 4> parameters = problem.materialLawParameters(cellGlobal, e, cellLocal);
      
      for (int i=0; i < size; i++) {
		  this->def[i] = 0;
		  saturationW[i] = 1.0 - sol[i][satNIdx];
		  pC[i] = problem.materialLaw().pC(saturationW[i], parameters);
		  pN[i] = sol[i][pWIdx] + pC[i];
		  mobility[i][pWIdx] = problem.materialLaw().mobW(saturationW[i], parameters);
		  mobility[i][satNIdx] = problem.materialLaw().mobN(sol[i][satNIdx], parameters);
      }
      density[pWIdx] = problem.materialLaw().wettingPhase.density();
      density[satNIdx] = problem.materialLaw().nonwettingPhase.density();
      
	  // ASSUMING element-wise constant permeability, evaluate K at the cell center 
	  const FieldMatrix<DT,n,n> K = problem.K(cellGlobal, e, cellLocal);  
	  
	  // ASSUMING constant element mapping jacobian
	  FieldMatrix<DT,n,n> jacobianInverseTransposed = geometry.jacobianInverseTransposed(cellLocal);

	  // ASSUMING element-wise constant porosity, evaluate at the cell center 
	  double volumeFactor = problem.porosity(cellGlobal, e, cellLocal)*cellVolume;
	  for (int i=0; i < size; i++) // begin loop over vertices
	  {
		  // current nonwetting phase saturation 
		  RT satNI = sol[i][satNIdx];
		  RT satNOld = uold[i][satNIdx];
		  
		  // time derivative
		  RT diffSatN = volumeFactor*(satNI - satNOld);
		  this->def[i][pWIdx] -= diffSatN;
		  this->def[i][satNIdx] += diffSatN;

		  // local coordinate of vertex 
	      const FieldVector<DT,n> vertexLocal = sfs[i].position();
	      
		  // get global coordinate of vertex
		  const FieldVector<DT,n> vertexGlobal = geometry.global(vertexLocal);
	      
		  // get source term 
		  FieldVector<RT, m> q = problem.q(vertexGlobal, e, vertexLocal);
		  
		  // add source to defect 
		  q *= dt*cellVolume;
		  this->def[i] -= q;
		  
		  for (int j=i+1; j < size; j++) // begin loop over neighboring vertices 
		  {
			  // local coordinate of neighbor 
			  const FieldVector<DT,n> neighborLocal = sfs[j].position();

			  if (!gt.isSimplex()) {
				  // compute the local distance 
				  DT distanceLocal = (vertexLocal - neighborLocal).two_norm();
			  
				  // check whether the two vertices share a cell edge
				  if (distanceLocal > 1.01)
					  continue;
			  }
			  
			  // get the local edge center 
			  FieldVector<DT,n> edgeLocal = vertexLocal + neighborLocal;
			  edgeLocal *= 0.5;
			  
			  // get global coordinate of neighbor
			  const FieldVector<DT,n> neighborGlobal = geometry.global(neighborLocal);
			  
			  // compute the edge vector
			  FieldVector<DT,n>  edgeVector = neighborGlobal - vertexGlobal;
			  
			  // get distance between neighbors 
			  DT oneByDistanceGlobal = 1.0/edgeVector.two_norm(); 
			  
			  // normalize edge vector 
			  edgeVector *= oneByDistanceGlobal; 
			  
			  // permeability in edge direction 
			  FieldVector<DT,n> Kij(0);
			  K.umv(edgeVector, Kij);
			  
			  // calculate pressure differences 
//			  VBlockType pDiff;
//			  pDiff[pWIdx] = sol[j][pWIdx] - sol[i][pWIdx];
//			  pDiff[satNIdx] = pDiff[pWIdx] + pC[j] - pC[i];
			  
			  VBlockType flux;
			  for (int comp = 0; comp < m; comp++) {
				  // calculate pressure component gradient
//				  FieldVector<RT, n> pGrad(edgeVector);
//				  pGrad *= oneByDistanceGlobal*pDiff[comp];
				  
		          // calculate FE gradient
		          FieldVector<RT, n> pGrad(0);
		          for (int k = 0; k < size; k++) {
		        	  FieldVector<DT,n> grad(0),temp;
		        	  for (int l = 0; l < n; l++) 
		        		  temp[l] = sfs[k].evaluateDerivative(0, l, edgeLocal);
		        	  jacobianInverseTransposed.umv(temp, grad);
		        	  grad *= (comp) ? pN[k] : sol[k][pWIdx];
		        	  pGrad += grad;
		          }
				  
				  // adjust by gravity 
				  FieldVector<RT, n> gravity = problem.gravity();
				  gravity *= density[comp];
				  pGrad -= gravity;
				  
				  // calculate the flux using upwind
				  RT outward = pGrad*Kij;
				  if (outward < 0)
					  flux[comp] = mobility[i][comp]*outward;
				  else 
					  flux[comp] = mobility[j][comp]*outward;
			  }
			  
			  // get global coordinate of edge center
			  const FieldVector<DT,n> edgeGlobal = geometry.global(edgeLocal);
			  
			  DT dualFaceArea; 
			  FieldVector<DT,n> leftFaceLocal, leftFaceGlobal, rightFaceLocal, rightFaceGlobal;
			  int leftFaceIdx; 
			  int rightFaceIdx; 
			  
			  switch (n) {
			  case 1:
				  dualFaceArea = 1.0;
				  break;
			  case 2:
				  // distance between cell center and edge center
				  dualFaceArea = (cellGlobal - edgeGlobal).two_norm();
				  break;
			  case 3: 
				  switch (size) {
				  case 8:
					  leftFaceIdx = edgeToFaceHex[i][j];
					  rightFaceIdx = edgeToFaceHex[j][i];
					  break;
				  default:
					  DUNE_THROW(NotImplemented, "BoxPwSnJacobian :: localDefect for dim = " << n 
							  << " and sfs.size = " << size); 
					  break;
				  }

				  // get global coordinates of face centers
			      leftFaceLocal = ReferenceElements<DT,n>::general(gt).position(leftFaceIdx,1);
				  leftFaceGlobal = geometry.global(leftFaceLocal);
			      rightFaceLocal = ReferenceElements<DT,n>::general(gt).position(rightFaceIdx,1);
				  rightFaceGlobal = geometry.global(rightFaceLocal);
				  
				  // calculate area of dual face
				  dualFaceArea = quadrilateralArea(edgeGlobal, rightFaceGlobal, cellGlobal, leftFaceGlobal);
				  break;
			  default:
				  DUNE_THROW(NotImplemented, "BoxPwSnJacobian :: localDefect for dim = " << n); 
				  break;
			  }

			  // obtain integrated Flux 
			  flux *= dualFaceArea; 
			  
			  // add to defect 
			  flux *= dt;
			  this->def[i] -= flux;
			  this->def[j] += flux;
		  } // end loop over neighboring vertices 
	  } // end loop over vertices

	  // adjust by density 
	  for (int i=0; i < size; i++) {
		  this->def[i][pWIdx] *= density[pWIdx];
		  this->def[i][satNIdx] *= density[satNIdx];
	  }
	  
	  // assemble boundary conditions 
	  assembleBC<TypeTag> (e); 
	  
	  // add to defect 
	  for (int i=0; i < size; i++) {
	    this->b[i] *= dt;
		  this->def[i] -= this->b[i];
	  }

//	  for (int i=0; i < 2*size; i++) {
//		  std::cout << "M = " << mass[i] << ", A = " << stiffness[i] << ", Q = " << neumann[i] << std::endl;
//	  }
//	  std::cout << std::endl;
		
      return;
    }
    
    void setLocalSolution (const Entity& e)
    {
        GeometryType gt = e.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
        	sfs=LagrangeShapeFunctions<DT,RT,n>::general(gt, 1);
        int size = sfs.size();
        this->setcurrentsize(size);
        
        for (int i = 0; i < size; i++) 
        	for (int comp = 0; comp < m; comp++) {
        		this->u[i][comp] = currentSolution.evallocal(comp, e, sfs[i].position());
        		uold[i][comp] = oldSolution.evallocal(comp, e, sfs[i].position());
        }
        
    	return;
    }
    
    void setDt (double d) 
    {
    	dt = d;
    }
  
    void setOldSolution (BoxFunction& uOld) 
    {
    	*oldSolution = *uOld;
    }
  
    template<class TypeTag>
    void analyticJacobian (const Entity& e, const VBlockType* sol)
    {
      // extract some important parameters
      const Geometry& geometry = e.geometry();
      GeometryType gt = geometry.type();
      const typename LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
      	sfs=LagrangeShapeFunctions<DT,RT,n>::general(gt, 1);
      int size = sfs.size();
      this->setcurrentsize(size);
      
	  std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
      // get cell volume
      DT cellVolume = geometry.volume();

      // assuming rectangular grids:
      cellVolume /= 4.0;
      
      // cell center in reference element
      const FieldVector<DT,n> 
		  cellLocal = ReferenceElements<DT,n>::general(gt).position(0,0);
		
	  // get global coordinate of cell center
	  const FieldVector<DT,n> cellGlobal = geometry.global(cellLocal);
      
	  // calculate secondary variables and set defect to zero 
      RT saturationW[size];
      RT pC[size], pN[size];
      VBlockType mobility[size];
      VBlockType dmobW[size];
      VBlockType dmobN[size];
      RT dpCSw[size];
      VBlockType density;
      // ASSUME element-wise constant parameters for the material law 
      FieldVector<RT, 4> parameters = problem.materialLawParameters(cellGlobal, e, cellLocal);
      
      for (int i=0; i < size; i++) {
		  saturationW[i] = 1.0 - sol[i][satNIdx];
		  pC[i] = problem.materialLaw().pC(saturationW[i], parameters);
		  pN[i] = sol[i][pWIdx] + pC[i];
		  mobility[i][pWIdx] = problem.materialLaw().mobW(saturationW[i], parameters);
		  mobility[i][satNIdx] = problem.materialLaw().mobN(sol[i][satNIdx], parameters);
		  RT epsSn = 1e-7*(fabs(sol[i][satNIdx]) + 1.0);
		  if (sol[i][satNIdx] > 0.99999) 
		    epsSn *= -1;
		  RT sNPlusEps = sol[i][satNIdx] + epsSn;
		  RT sWPlusEps = 1 - sNPlusEps;
		  RT mobWPlusEps = problem.materialLaw().mobW(sWPlusEps, parameters);
		  RT mobNPlusEps = problem.materialLaw().mobN(sNPlusEps, parameters);
		  RT pCPlusEps = problem.materialLaw().pC(sWPlusEps, parameters);
		  dmobW[i][pWIdx] = 0; // taken from UG's box2p: no explicit dependence on pW ??????
		  dmobW[i][satNIdx] = (mobWPlusEps - mobility[i][pWIdx])/epsSn;
		  dmobN[i][pWIdx] = 0;
		  dmobN[i][satNIdx] = (mobNPlusEps - mobility[i][satNIdx])/epsSn;
		  dpCSw[i] = (pCPlusEps - pC[i])/epsSn;
		  //std::cout << "i = " << i << ", pW = " << sol[i][pWIdx] << ", pN = " << pN[i] << std::endl;
		  for (int j = 0; j < size; j++)
		    this->A[i][j] = 0;
      }
      density[pWIdx] = problem.materialLaw().wettingPhase.density();
      density[satNIdx] = problem.materialLaw().nonwettingPhase.density();
      
          // ASSUMING constant element mapping jacobian
	  FieldMatrix<DT,n,n> jacobianInverseTransposed = geometry.jacobianInverseTransposed(cellLocal);


	  // ASSUMING element-wise constant permeability, evaluate K at the cell center 
	  const FieldMatrix<DT,n,n> K = problem.K(cellGlobal, e, cellLocal);  
	  
	  // ASSUMING element-wise constant porosity, evaluate at the cell center 
	  double volumeFactor = problem.porosity(cellGlobal, e, cellLocal)*cellVolume;
	  for (int i=0; i < size; i++) // begin loop over vertices
	  {
		  // mass part
		  this->A[i][i][pWIdx][satNIdx] = -density[pWIdx]*volumeFactor;
		  this->A[i][i][satNIdx][satNIdx] = density[satNIdx]*volumeFactor;

		  // local coordinate of vertex 
		  const FieldVector<DT,n> vertexLocal = sfs[i].position();
	      
		  // get global coordinate of vertex
		  const FieldVector<DT,n> vertexGlobal = geometry.global(vertexLocal);
		  
		  for (int j=0; j < size; j++) // begin loop over neighboring vertices 
		  {
		    if (i == j) 
		      continue;
			  // local coordinate of neighbor 
			  const FieldVector<DT,n> neighborLocal = sfs[j].position();

			  if (!gt.isSimplex()) {
				  // compute the local distance 
				  DT distanceLocal = (vertexLocal - neighborLocal).two_norm();
			  
				  // check whether the two vertices share a cell edge
				  if (distanceLocal > 1.01)
					  continue;
			  }
			  
			  // get global coordinate of neighbor
			  const FieldVector<DT,n> neighborGlobal = geometry.global(neighborLocal);
			  
			  // compute the edge vector
			  FieldVector<DT,n>  edgeVector = neighborGlobal - vertexGlobal;
			  
			  // get distance between neighbors 
			  DT oneByDistanceGlobal = 1.0/edgeVector.two_norm(); 
			  
			  // normalize edge vector 
			  edgeVector *= oneByDistanceGlobal; 
			  
			  // permeability in edge direction 
			  FieldVector<DT,n> Kij(0);
			  K.umv(edgeVector, Kij);
			  
			  // get the local edge center 
			  FieldVector<DT,n> edgeLocal = vertexLocal + neighborLocal;
			  edgeLocal *= 0.5;

			  // get the local coordinate of the control volume face
			  FieldVector<DT,n> controlVolumeFaceLocal = edgeLocal + cellLocal;
			  controlVolumeFaceLocal *= 0.5;

			  // calculate pressure differences 
			  //VBlockType pDiff;
			  //pDiff[pWIdx] = sol[j][pWIdx] - sol[i][pWIdx];
			  //pDiff[satNIdx] = pDiff[pWIdx] + pC[j] - pC[i];
			  
			  VBlockType flux;
			  int upwindNode[2];
			  VBlockType pGradKn;
			  FieldVector<RT, n> gradients[m];
			  //std::cout << "jacobianInverseTransposed = " << jacobianInverseTransposed;
			  //std::cout << "local = " << controlVolumeFaceLocal << ", grad pw = ";
			  for (int comp = 0; comp < m; comp++) {
				  // calculate pressure component gradient
				  // FieldVector<RT, n> pGrad(edgeVector);
				  // pGrad *= oneByDistanceGlobal*pDiff[comp];

			          // calculate FE gradient
			          FieldVector<RT, n> pGrad(0);
			          for (int k = 0; k < size; k++) {
				    FieldVector<DT,n> grad(0),temp;
				    for (int l = 0; l < n; l++) 
				      temp[l] = sfs[k].evaluateDerivative(0, l, controlVolumeFaceLocal);
				    jacobianInverseTransposed.umv(temp, grad);
				    //std::cout << "grad = " << grad << ", p = " << ((comp) ? pN[k] : sol[k][pWIdx]) << std::endl;
				    grad *= (comp) ? pN[k] : sol[k][pWIdx];
				    pGrad += grad;
				  }
				  //std::cout << pGrad << ", grad pn = ";
				  // adjust by gravity 
				  FieldVector<RT, n> gravity = problem.gravity();
				  gravity *= density[comp];
				  pGrad -= gravity;

				  gradients[comp] = 0;
				  K.umv(pGrad, gradients[comp]);
				  // select the upwind node 
 				  pGradKn[comp] = -(pGrad*Kij);
 				  if (pGradKn[comp] > 0)
				    upwindNode[comp] = i;
 				  else 
				    upwindNode[comp] = j;
				  //std::cout << "pGrad = " << pGrad << ", Kij = " << Kij << ", pGradKn = " << pGradKn[comp] << std::endl;
			  }
			  //std::cout << std::endl;      


			  // get the shapefunction derivative 
			  FieldVector<RT, n> sGrad, shapeGradient(0);
			  for (int k = 0; k < n; k++)
			    sGrad[k] = sfs[i].evaluateDerivative(0, k, controlVolumeFaceLocal);
			  jacobianInverseTransposed.umv(sGrad, shapeGradient);
			  RT shapeGradKnI = shapeGradient*Kij;
			  shapeGradient = 0;
			  for (int k = 0; k < n; k++)
			    sGrad[k] = sfs[j].evaluateDerivative(0, k, controlVolumeFaceLocal);
			  jacobianInverseTransposed.umv(sGrad, shapeGradient);
			  RT shapeGradKnJ = shapeGradient*Kij;
			  
			  
			  // get global coordinate of edge center
			  const FieldVector<DT,n> edgeGlobal = geometry.global(edgeLocal);
			  
			  // distance between cell center and edge center
			  DT distanceEdgeCell = (cellGlobal - edgeGlobal).two_norm();
			  
			  ////////////////////////////////////////////////////////////
			  // CAREFUL: only valid in 2D 
			  ////////////////////////////////////////////////////////////
			  // obtain integrated Flux 
			  //flux *= distanceEdgeCell; 
			  
			  RT entry = dt*distanceEdgeCell*density[pWIdx]*(dmobW[upwindNode[pWIdx]][pWIdx]*pGradKn[pWIdx]);
			  //std::cout << dmobW[upwindNode[pWIdx]][pWIdx] << ", " << pGradKn[pWIdx] << std::endl;;
			  RT entryIJ = 0;
			  RT entryJI = 0;
			  if (upwindNode[pWIdx] == i)
			    entryIJ = entry;
			  else 
			    entryJI = entry;
			  entryIJ += dt*distanceEdgeCell*density[pWIdx]*mobility[upwindNode[pWIdx]][pWIdx]*shapeGradKnI;
			  entryJI += dt*distanceEdgeCell*density[pWIdx]*mobility[upwindNode[pWIdx]][pWIdx]*shapeGradKnJ;
			  this->A[i][i][pWIdx][pWIdx] -= entryIJ;
			  this->A[j][i][pWIdx][pWIdx] += entryIJ;
			  //std::cout << dt << ", " << distanceEdgeCell << ", " << density[pWIdx] << ", " << 
			  //  mobility[upwindNode[pWIdx]][pWIdx] << ", " << shapeGradKnI << std::endl;;
			  //this->A[j][j][pWIdx][pWIdx] += entryJI;
			  //this->A[i][j][pWIdx][pWIdx] -= entryJI;
			  //RT g0 = vertexGlobal[0];
			  //RT g1 = vertexGlobal[1];
			  //RT n0 = neighborGlobal[0];
			  //RT n1 = neighborGlobal[1];
			  //			  if ((g0 == 2.0) && (g1 == 4.0) && (n0 == 2.5) && (n1 == 4.0)) {
			  //std::cout << "i = " << i << ", j = " << j 
			  //	    << ", pwI = " << sol[i][pWIdx] << ", sNI = " << sol[i][satNIdx] 
			  //	      << "pwJ = " << sol[j][pWIdx] << ", sNJ = " << sol[j][satNIdx] << std::endl; 
			  //std::cout << "lambda = " << mobility[upwindNode[pWIdx]][pWIdx] << ", sGrad = " << shapeGradient  
			  //	      << ", upwind = " << upwindNode[pWIdx] << ", entryIJ = " << entryIJ << ", entryJI = " << entryJI << std::endl;
			    //  }
// 			  std::cout << "node = " << vertexGlobal << ", neighbor = " << neighborGlobal << ", upwind = " << upwindNode[pWIdx];
// 			  std::cout << ", entryIJ = " << entryIJ << ", entryJI = " << entryJI << std::endl; 
			  entryIJ = entryJI = 0;
			  if (upwindNode[pWIdx] == i) {
			    entryIJ = dt*distanceEdgeCell*density[pWIdx]*(dmobW[upwindNode[pWIdx]][satNIdx]*pGradKn[pWIdx]);
			    this->A[i][i][pWIdx][satNIdx] += entryIJ;
			    this->A[j][i][pWIdx][satNIdx] -= entryIJ;
			  }
			  else {
			    //entryJI = dt*distanceEdgeCell*density[pWIdx]*(dmobW[upwindNode[pWIdx]][satNIdx]*pGradKn[pWIdx]);
			    //this->A[j][j][pWIdx][satNIdx] -= entryJI;
			    //this->A[i][j][pWIdx][satNIdx] += entryJI;
			  }
			  // if ((g0 == 2.0) && (g1 == 4.0) && (n0 == 2.5) && (n1 == 4.0)) {
			  //  std::cout << "lambda = " << dmobW[upwindNode[pWIdx]][satNIdx] << ", pGradKn = " << pGradKn[pWIdx] 
			  //	      << ", upwind = " << upwindNode[pWIdx] << ", entryIJ = " << entryIJ << ", entryJI = " << entryJI << std::endl;
			    //   }
			  //std::cout << ", entryIJ = " << entryIJ << ", entryJI = " << entryJI << std::endl; 
			  

			  entry = dt*distanceEdgeCell*density[satNIdx]*(dmobN[upwindNode[satNIdx]][pWIdx]*pGradKn[satNIdx]);
			  entryIJ = entryJI = 0;
			  if (upwindNode[satNIdx] == i)
			    entryIJ = entry;
			  else 
			    entryJI = entry;
			  entryIJ += dt*distanceEdgeCell*density[satNIdx]*mobility[upwindNode[satNIdx]][satNIdx]*shapeGradKnI;
			  entryJI += dt*distanceEdgeCell*density[satNIdx]*mobility[upwindNode[satNIdx]][satNIdx]*shapeGradKnJ;
			  this->A[i][i][satNIdx][pWIdx] -= entryIJ;
			  this->A[j][i][satNIdx][pWIdx] += entryIJ;
			  //this->A[j][j][satNIdx][pWIdx] += entryJI;
			  //this->A[i][j][satNIdx][pWIdx] -= entryJI;
// 	  std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
			  //if ((g0 == 2.0) && (g1 == 4.0) && (n0 == 2.5) && (n1 == 4.0)) {
			  //std::cout << "lambda = " << mobility[upwindNode[satNIdx]][satNIdx] 
			  //	      << ", upwind = " << upwindNode[satNIdx] << ", entryIJ = " << entryIJ << ", entryJI = " << entryJI << std::endl;
			    //  }
			  
			  entryIJ = entryJI = 0;
			  if (upwindNode[satNIdx] == i) {
			    entryIJ = dt*distanceEdgeCell*density[satNIdx]*(dmobN[upwindNode[satNIdx]][satNIdx]*pGradKn[satNIdx]
									    + mobility[upwindNode[satNIdx]][satNIdx]*shapeGradKnJ*dpCSw[i]);
			    this->A[i][i][satNIdx][satNIdx] += entryIJ;
			    this->A[j][i][satNIdx][satNIdx] -= entryIJ;
			  }
			  else {
			    entryJI = dt*distanceEdgeCell*density[satNIdx]*(dmobN[upwindNode[satNIdx]][satNIdx]*pGradKn[satNIdx]);
			    //this->A[j][j][satNIdx][satNIdx] -= entryJI;
			    //this->A[i][j][satNIdx][satNIdx] += entryJI;
			  }
			  // if ((g0 == 2.0) && (g1 == 4.0) && (n0 == 2.5) && (n1 == 4.0)) {
			  //std::cout << "lambda_ = " << dmobN[upwindNode[satNIdx]][satNIdx] << ", spn = " << distanceEdgeCell*pGradKn[satNIdx]
			  //	      << ", lambda = " << mobility[upwindNode[satNIdx]][satNIdx] << ", sp_ = " << distanceEdgeCell*shapeGradKnJ*dpCSw[i]
			  //	      << ", upwind = " << upwindNode[satNIdx] << ", entryIJ = " << entryIJ << ", entryJI = " << entryJI << std::endl;
			    //  }
			  //std::cout << ", entryIJ = " << entryIJ << ", entryJI = " << entryJI << std::endl; 

		  } // end loop over neighboring vertices 
	  } // end loop over vertices

      return;

    }

    template<class TypeTag>
    void assembleBC (const Entity& e)
		{
		  // extract some important parameters
		  GeometryType gt = e.geometry().type();
		  const typename LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
			sfs=LagrangeShapeFunctions<DT,RT,n>::general(gt,1);
		  setcurrentsize(sfs.size());

		  for (int i = 0; i < sfs.size(); i++) {
			  this->b[i] = 0;
			  this->bctype[i].assign(BoundaryConditions::neumann);
		  }
		  
		  // determine quadrature order
		  int p=0;
		  // evaluate boundary conditions via intersection iterator
		  typedef typename IntersectionIteratorGetter<G,TypeTag>::IntersectionIterator
		    IntersectionIterator;
		  
		  IntersectionIterator endit = IntersectionIteratorGetter<G,TypeTag>::end(e);
		  for (IntersectionIterator it = IntersectionIteratorGetter<G,TypeTag>::begin(e); 
		       it!=endit; ++it)
			{
			  // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
			  // in level assemble treat non-level neighbors as boundary
			  if (it.neighbor())
				{
				  if (levelBoundaryAsDirichlet && it.outside()->level()==e.level()) 
					continue;
				  if (!levelBoundaryAsDirichlet)
					continue;
				}

			  // determine boundary condition type for this face, initialize with processor boundary
			  FieldVector<typename BoundaryConditions::Flags, m> bctypeface(BoundaryConditions::process);

			  // handle face on exterior boundary, this assumes there are no interior boundaries
			  if (it.boundary())
				{
				  GeometryType gtface = it.intersectionSelfLocal().type();
				  for (size_t g=0; g<QuadratureRules<DT,n-1>::rule(gtface,p).size(); ++g)
					{
					  const FieldVector<DT,n-1>& facelocal = QuadratureRules<DT,n-1>::rule(gtface,p)[g].position();
					  FieldVector<DT,n> local = it.intersectionSelfLocal().global(facelocal);
					  FieldVector<DT,n> global = it.intersectionGlobal().global(facelocal);
					  bctypeface = problem.bctype(global,e,it,local); // eval bctype


					  if (bctypeface[0]!=BoundaryConditions::neumann) break;

					  VBlockType J = problem.J(global,e,it,local);
					  if (J.two_norm() < 1e-10) 
						  continue;
					  double weightface = QuadratureRules<DT,n-1>::rule(gtface,p)[g].weight();
					  DT detjacface = it.intersectionGlobal().integrationElement(facelocal);
					  J *= 1.0/(pow(2.0, n-1))*weightface*detjacface;
					  for (int i=0; i<sfs.size(); i++) // loop over test function number
						if (this->bctype[i][0]==BoundaryConditions::neumann)
						  {
							//////////////////////////////////////////////////////////////////////////
							// HACK: piecewise constants with respect to dual grid not implemented yet 
							// works only if exactly one quadrature point is located within each dual 
							// cell boundary (which should be the case for p = 2)
							//////////////////////////////////////////////////////////////////////////
							//if (sfs[i].evaluateFunction(0,local) > 0.5) {
								this->b[i] -= J;
							//}
						  }
					}
				  if (bctypeface[0]==BoundaryConditions::neumann) continue; // was a neumann face, go to next face
				}

			  // If we are here, then it is 
			  // (i)   an exterior boundary face with Dirichlet condition, or
			  // (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
			  // (iii) a level boundary in case of level-wise assemble
			  // How processor boundaries are handled depends on the processor boundary mode
			  if (bctypeface[0]==BoundaryConditions::process && procBoundaryAsDirichlet==false 
				  && levelBoundaryAsDirichlet==false) 
				continue; // then it acts like homogeneous Neumann

			  // now handle exterior or interior Dirichlet boundary
			  for (int i=0; i<sfs.size(); i++) // loop over test function number
				{
				  if (sfs[i].codim()==0) continue; // skip interior dof
				  if (sfs[i].codim()==1) // handle face dofs
					{
					  if (sfs[i].entity()==it.numberInSelf())
						{
						  if (this->bctype[i][0]<bctypeface[0])
							{
							  this->bctype[i].assign(bctypeface[0]);
							  if (bctypeface[0]==BoundaryConditions::process)
								this->b[i] = 0;
							  if (bctypeface[0]==BoundaryConditions::dirichlet)
								{
								  this->b[i] = 0;
								}
							}
						}
					  continue;
					}
				  // handle subentities of this face
				  for (int j=0; j<ReferenceElements<DT,n>::general(gt).size(it.numberInSelf(),1,sfs[i].codim()); j++)
					if (sfs[i].entity()==ReferenceElements<DT,n>::general(gt).subEntity(it.numberInSelf(),1,j,sfs[i].codim()))
					  {
						if (this->bctype[i][0]<bctypeface[0])
						  {
						    this->bctype[i].assign(bctypeface[0]);
							if (bctypeface[0]==BoundaryConditions::process)
							  this->b[i] = 0;
							if (bctypeface[0]==BoundaryConditions::dirichlet)
							  {
								this->b[i] = 0;
							  }
						  }
					  }
				}
			}
		}


  private:
    // parameters given in constructor
    TwoPhaseProblem<G,RT>& problem;
    bool levelBoundaryAsDirichlet;
    bool procBoundaryAsDirichlet;
    const BoxFunction& currentSolution;
    BoxFunction oldSolution;
  public:
    double dt;
    VBlockType uold[SIZE];
  };

  /** @} */
}
#endif
