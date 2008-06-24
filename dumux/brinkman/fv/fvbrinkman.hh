#ifndef DUNE_FVBRINKMAN_HH
#define DUNE_FVBRINKMAN_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/brinkman/brinkman.hh"

/**
 * @file
 * @brief  Finite Volume Brinkman Model
 * @author Bernd Flemisch, Jochen Fritz
 */

namespace Dune
{
  //! \ingroup diffusion
  //! Finite Volume Brinkman Model
  /*! Provides a Finite Volume implementation for the evaluation 
   * of equations of the form 
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$, 
   * \f$p = g\f$ on \f$\Gamma_1\f$, and 
   * \f$\lambda K \text{grad}\, p \cdot \mathbf{n} = J\f$ 
   * on \f$\Gamma_2\f$. Here, 
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability, 
   * and \f$\lambda\f$ the total mobility, possibly depending on the 
   * saturation, \f$q\f$ the source term.
	Template parameters are:

	- G         a DUNE grid type
	- RT        type used for return values 
   */
  template<class G, class RT>
  class FVBrinkman 
  : public Brinkman< G, RT, BlockVector< FieldVector<RT,1> >, BlockVector< FieldVector<RT,G::dimension> > >
  {
	  template<int dim>
	  struct ElementLayout
	  {
		  bool contains (GeometryType gt)
	      {
			  return gt.dim() == dim;
	      }
	  }; 
	  
	  enum{dim = G::dimension};	
	  enum{dimworld = G::dimensionworld};
	  
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename G::Traits::LeafIndexSet IS;
	  typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
	  typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	  typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	  typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
	  typedef typename G::ctype ct; 
	  typedef FieldMatrix<double,1,1> MB;
	  typedef BCRSMatrix<MB> PressureMatrixType;
	  typedef FieldMatrix<double,dim,dim> MBV;
	  typedef BCRSMatrix<MBV> VelocityMatrixType;
	  typedef FieldVector<double, 1> VB;
	  typedef BlockVector<VB> Vector;
	  typedef BlockVector< FieldVector<double, dim> > VelType;
	  
  public:
	typedef BlockVector< FieldVector<RT,1> > RepresentationType;

	void updateVelocityRHS()
	{
		fV = fVBound; 
		
		Iterator eendit = indexset.template end<0,All_Partition>();
		  for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		  {		
			  // cell geometry type
			  GeometryType gt = it->geometry().type();

			  // cell center in reference element
			  const FieldVector<ct,dim>& 
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0);

			  // get global coordinate of cell center
			  FieldVector<ct,dim> global = it->geometry().global(local);

			  // cell index
			  int indexi = elementmapper.map(*it);

			  RT pI = this->pressure[indexi];
			  
			  IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			  for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
			  is!=endit; ++is)
			  {

				  // get geometry type of face
				  GeometryType gtf = is->intersectionSelfLocal().type();

				  // center in face's reference element
				  const FieldVector<ct,dim-1>& 
				  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

				  // get normal vector scaled with volume
				  FieldVector<ct,dimworld> integrationOuterNormal 
				  = is->integrationOuterNormal(facelocal);

				  FieldVector<ct,dimworld> entry(integrationOuterNormal); 

				  // handle interior face
				  if (is->neighbor()) 
				  {
					  // access neighbor
					  EntityPointer outside = is->outside();
					  int indexj = elementmapper.map(*outside);
					  
					  RT pJ = this->pressure[indexj];
					  
					  entry *= 0.5*(pI + pJ);
				  }
				  else 
				  {
					  // get geometry type of face
					  GeometryType gtf = is->intersectionSelfLocal().type();

					  // center in face's reference element
					  const FieldVector<ct,dim-1>& 
					  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

					  // center of face inside volume reference element
					  const FieldVector<ct,dim>& 
					  facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is->numberInSelf(),1);

					  // center of face in global coordinates
					  FieldVector<ct,dimworld> 
					  faceglobal = is->intersectionGlobal().global(facelocal);

					  //get boundary condition for boundary face center
					  BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);

					  // Dirichlet for velocity is Neumann for pressure and vice versa 
					  if (bctype == BoundaryConditions::dirichlet) 
					  { 
						  FieldVector<ct,dimworld> distVec(global - faceglobal);
						  double dist = distVec.two_norm();

						  RT dPdn = this->problem.JPressure(faceglobal, *it, facelocalDim);
						  std::cout << "normal = " << entry << ", pI = " << pI << 
						  ", dist = " << dist << ", dPdN = " << dPdn << std::endl;
						  entry *= (pI + dist*dPdn);
					  }
					  else
					  {
						  RT pB = this->problem.gPressure(faceglobal, *it, facelocalDim);
						  
						  entry *= pB;
					  }
					  
				  }
			  
				  fV[indexi] -= entry;
			  }
		  }
	}

	template<class IFPressure>
	void calculateInterfacePressures(IFPressure& ifPressure)
	{
		Iterator eendit = indexset.template end<0,All_Partition>();
		  for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		  {		
			  // cell geometry type
			  GeometryType gt = it->geometry().type();

			  // cell center in reference element
			  const FieldVector<ct,dim>& 
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0);

			  // get global coordinate of cell center
			  FieldVector<ct,dim> global = it->geometry().global(local);

			  // cell index
			  int indexi = elementmapper.map(*it);

			  RT pI = this->pressure[indexi];

			  IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			  for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
			  is!=endit; ++is)
			  {
				  int faceIdx = is->numberInSelf();
				  
				  // handle interior face
				  if (is->neighbor()) 
				  {
					  // access neighbor
					  EntityPointer outside = is->outside();
					  int indexj = elementmapper.map(*outside);
					  
					  RT pJ = this->pressure[indexj];
					  
					  ifPressure[indexi][faceIdx] = 0.5*(pI + pJ);
				  }
				  else 
				  {
					  // get geometry type of face
					  GeometryType gtf = is->intersectionSelfLocal().type();

					  // center in face's reference element
					  const FieldVector<ct,dim-1>& 
					  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

					  // center of face inside volume reference element
					  const FieldVector<ct,dim>& 
					  facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is->numberInSelf(),1);

					  // center of face in global coordinates
					  FieldVector<ct,dimworld> 
					  faceglobal = is->intersectionGlobal().global(facelocal);

					  //get boundary condition for boundary face center
					  BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);

					  // Dirichlet for velocity is Neumann for pressure and vice versa 
					  if (bctype == BoundaryConditions::dirichlet) 
					  { 
						  FieldVector<ct,dimworld> distVec(global - faceglobal);
						  double dist = distVec.two_norm();

						  RT dPdn = this->problem.JPressure(faceglobal, *it, facelocalDim);

						  ifPressure[indexi][faceIdx] = (pI + dist*dPdn);
					  }
					  else
					  {
						  RT pB = this->problem.gPressure(faceglobal, *it, facelocalDim);
						  
						  ifPressure[indexi][faceIdx] = pB;
					  }
					  
				  }
			  }
		  }
	}
	
	void updatePressureRHS()
	{
	    const int blocksize = 2*dim;
	    typedef BlockVector< FieldVector<RT, blocksize> > IFPressure;
	    IFPressure ifPressure(this->grid.size(0));
	    calculateInterfacePressures(ifPressure);
	    
		fP = 0; 
		
		Iterator eendit = indexset.template end<0,All_Partition>();
		  for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		  {		
			  // cell geometry type
			  GeometryType gt = it->geometry().type();

			  // cell center in reference element
			  const FieldVector<ct,dim>& 
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0);

			  // cell volume 
			  double volume = it->geometry().integrationElement(local)
			  *ReferenceElements<ct,dim>::general(gt).volume();

			  // get global coordinate of cell center
			  FieldVector<ct,dim> global = it->geometry().global(local);

			  // cell index
			  int indexi = elementmapper.map(*it);

			  RT pI = this->pressure[indexi];
			  FieldVector<RT,dim> velI = this->velocity[indexi];
			  
			  IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			  for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
			  is!=endit; ++is)
			  {
				  // eastIndex for 0,1: 0, for 2,3: 2, for 4,5: 4
				  int faceIdx = is->numberInSelf(); 
				  int eastIndex = faceIdx/2; 
				  eastIndex *= 2; 
				  int westIndex = eastIndex + 1;

				  // get geometry type of face
				  GeometryType gtf = is->intersectionSelfLocal().type();

				  // center in face's reference element
				  const FieldVector<ct,dim-1>& 
				  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

				  // get normal vector scaled with volume
				  FieldVector<ct,dimworld> integrationOuterNormal 
				  = is->integrationOuterNormal(facelocal);

				  FieldVector<ct,dimworld> velFace; 

				  // handle interior face
				  if (is->neighbor()) 
				  {
					  // access neighbor
					  EntityPointer outside = is->outside();
					  int indexj = elementmapper.map(*outside);
					  
					  // compute factor in neighbor
					  GeometryType nbgt = outside->geometry().type();
					  const FieldVector<ct,dim>& 
					  nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);

					  // neighbor cell center in global coordinates
					  FieldVector<ct,dimworld> 
					  nbglobal = outside->geometry().global(nblocal);

					  // distance vector between barycenters
					  FieldVector<ct,dimworld> 
					  distVec = global - nbglobal;

					  // compute distance between cell centers
					  double dist = distVec.two_norm();

					  RT pJ = this->pressure[indexj];
					  FieldVector<RT,dim> velJ = this->velocity[indexj];
					  
					  RT pDiff;
					  if (faceIdx%2) 
						  pDiff = pJ - pI;
					  else 
						  pDiff = pI - pJ;						  
					  
					  for (int comp = 0; comp < dim; comp++)
						  velFace[comp] = 0.5*(velI[comp] + volume/(dist*AV[indexi][indexi][comp][comp])*(ifPressure[indexi][eastIndex] - ifPressure[indexi][westIndex]))
							  + 0.5*(velJ[comp] + volume/(dist*AV[indexj][indexj][comp][comp])*(ifPressure[indexj][eastIndex] - ifPressure[indexj][westIndex]))
							  - 0.5*volume/dist*pDiff*(1.0/AV[indexi][indexi][comp][comp] + 1.0/AV[indexj][indexj][comp][comp]);
				  }
				  else 
				  {
					  // get geometry type of face
					  GeometryType gtf = is->intersectionSelfLocal().type();

					  // center in face's reference element
					  const FieldVector<ct,dim-1>& 
					  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

					  // center of face inside volume reference element
					  const FieldVector<ct,dim>& 
					  facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is->numberInSelf(),1);

					  // center of face in global coordinates
					  FieldVector<ct,dimworld> 
					  faceglobal = is->intersectionGlobal().global(facelocal);

					  //get boundary condition for boundary face center
					  BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);

					  // Dirichlet for velocity is Neumann for pressure and vice versa 
					  if (bctype == BoundaryConditions::neumann) 
					  { 
						  velFace = velI;
					  }
					  else
					  {
						  velFace = this->problem.gVelocity(faceglobal, *it, facelocalDim);
					  }
				  }
			  
				  fP[indexi] -= velFace*integrationOuterNormal;
			  }
		  }
	}

	void computeVelocityCorrection()
	{
		this->velocityCorrection = 0;
		
		Iterator eendit = indexset.template end<0,All_Partition>();
		  for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		  {		
			  // cell geometry type
			  GeometryType gt = it->geometry().type();

			  // cell center in reference element
			  const FieldVector<ct,dim>& 
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0);

			  // cell volume 
			  double volume = it->geometry().integrationElement(local)
			  *ReferenceElements<ct,dim>::general(gt).volume();

			  // get global coordinate of cell center
			  FieldVector<ct,dim> global = it->geometry().global(local);

			  // cell index
			  int indexi = elementmapper.map(*it);

			  RT pCI = this->pressureCorrection[indexi];

			  IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			  for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
			  is!=endit; ++is)
			  {
				  int faceIdx = is->numberInSelf(); 

				  // get geometry type of face
				  GeometryType gtf = is->intersectionSelfLocal().type();

				  // center in face's reference element
				  const FieldVector<ct,dim-1>& 
				  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

				  // handle interior face
				  if (is->neighbor()) 
				  {
					  // access neighbor
					  EntityPointer outside = is->outside();
					  int indexj = elementmapper.map(*outside);
					  
					  // compute factor in neighbor
					  GeometryType nbgt = outside->geometry().type();
					  const FieldVector<ct,dim>& 
					  nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);

					  // neighbor cell center in global coordinates
					  FieldVector<ct,dimworld> 
					  nbglobal = outside->geometry().global(nblocal);

					  // distance vector between barycenters
					  FieldVector<ct,dimworld> 
					  distVec = global - nbglobal;

					  // compute distance between cell centers
					  double dist = distVec.two_norm();

					  RT pCJ = this->pressureCorrection[indexj];
					  
					  RT pCAvg = 0.5*(pCJ + pCI);

					  switch (faceIdx) {
					  case 0: 
						  this->velocityCorrection[indexi][0] += volume/(dist*AV[indexi][indexi][0][0])*pCAvg;
						  break;
					  case 1: 
						  this->velocityCorrection[indexi][0] -= volume/(dist*AV[indexi][indexi][0][0])*pCAvg;
						  break;
					  case 2: 
						  this->velocityCorrection[indexi][1] += volume/(dist*AV[indexi][indexi][1][1])*pCAvg;
						  break;
					  case 3: 
						  this->velocityCorrection[indexi][1] -= volume/(dist*AV[indexi][indexi][1][1])*pCAvg;
						  break;
					  case 4: 
						  this->velocityCorrection[indexi][2] += volume/(dist*AV[indexi][indexi][2][2])*pCAvg;
						  break;
					  case 5: 
						  this->velocityCorrection[indexi][2] -= volume/(dist*AV[indexi][indexi][2][2])*pCAvg;
						  break;
					  }
				  }
				  else 
				  {
					  // get geometry type of face
					  GeometryType gtf = is->intersectionSelfLocal().type();

					  // center in face's reference element
					  const FieldVector<ct,dim-1>& 
					  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

					  // center of face inside volume reference element
					  const FieldVector<ct,dim>& 
					  facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is->numberInSelf(),1);

					  // center of face in global coordinates
					  FieldVector<ct,dimworld> 
					  faceglobal = is->intersectionGlobal().global(facelocal);

					  //get boundary condition for boundary face center
					  BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);

					  // Dirichlet for velocity is Neumann for pressure and vice versa 
					  if (bctype == BoundaryConditions::dirichlet) 
					  { 
						  FieldVector<ct,dimworld> distVec(global - faceglobal);
						  double dist = distVec.two_norm();

						  switch (faceIdx) {
						  case 0: 
							  this->velocityCorrection[indexi][0] += volume/(2.0*dist*AV[indexi][indexi][0][0])*pCI;
							  break;
						  case 1: 
							  this->velocityCorrection[indexi][0] -= volume/(2.0*dist*AV[indexi][indexi][0][0])*pCI;
							  break;
						  case 2: 
							  this->velocityCorrection[indexi][1] += volume/(2.0*dist*AV[indexi][indexi][1][1])*pCI;
							  break;
						  case 3: 
							  this->velocityCorrection[indexi][1] -= volume/(2.0*dist*AV[indexi][indexi][1][1])*pCI;
							  break;
						  case 4: 
							  this->velocityCorrection[indexi][2] += volume/(2.0*dist*AV[indexi][indexi][2][2])*pCI;
							  break;
						  case 5: 
							  this->velocityCorrection[indexi][2] -= volume/(2.0*dist*AV[indexi][indexi][2][2])*pCI;
							  break;
						  }
					  }
				  }
			  }
		  }
	}

	void solveVelocitySystem(); 

	void solvePressureSystem(); 

	void computeVelocity()
	{
		updateVelocityRHS();
		solveVelocitySystem();
		return;
	}

	void computePressureCorrection()
	{
		updatePressureRHS();
		solvePressureSystem();
		return;
	}

	void vtkout (const char* name, int k) const 
	{
		VTKWriter<G> vtkwriter(this->grid);
		char fname[128];	
		sprintf(fname,"%s-%05d",name,k);
		vtkwriter.addCellData(this->pressure,"total pressure p~");
		vtkwriter.write(fname,VTKOptions::ascii);		
	}
	
	void initializeMatrices();
	
	void assembleMatrices();
	
	FVBrinkman(G& g, BrinkmanProblem<G, RT>& prob)
	            : Brinkman<G, RT, RepresentationType, VelType>(g, prob), 
	              elementmapper(g, g.leafIndexSet()), 
	              indexset(g.leafIndexSet()), 
	              AV(g.size(0), g.size(0), (2*dim+1)*g.size(0), BCRSMatrix<MBV>::random), 
	              AP(g.size(0), g.size(0), (2*dim+1)*g.size(0), BCRSMatrix<MB>::random), 
	              fP(g.size(0)), fV(g.size(0)), fVBound(g.size(0))
	{
		this->pressure.resize(g.size(0));
		this->pressure = 0;
		this->pressureCorrection.resize(g.size(0));
		this->pressureCorrection = 0;
		this->velocity.resize(g.size(0));
		this->velocity = 0;
		this->velocityCorrection.resize(g.size(0));
		this->velocityCorrection = 0;
		fV = 0;
		fVBound = 0;
		initializeMatrices();
		assembleMatrices();
	}
	
  //private:
	  EM elementmapper;
	  const IS& indexset;
	  VelocityMatrixType AV;
	  PressureMatrixType AP;
	  RepresentationType fP;
	  VelType fV;
	  VelType fVBound;
  };

  
  
  template<class G, class RT>
  void FVBrinkman<G, RT>::initializeMatrices()
  {
	    // determine matrix row sizes 
	    Iterator eendit = indexset.template end<0,All_Partition>();
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	      {
			// cell index
			int indexi = elementmapper.map(*it);
	
			// initialize row size
			int rowSize = 1;
	
			// run through all intersections with neighbors 
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
				  is!=endit; ++is)
			    if (is->neighbor()) 
			      rowSize++;
			AV.setrowsize(indexi, rowSize);
			AP.setrowsize(indexi, rowSize);
	      }
	    AV.endrowsizes();
	    AP.endrowsizes();

	    // determine position of matrix entries 
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	      {
			// cell index
			int indexi = elementmapper.map(*it);
	
			// add diagonal index
			AV.addindex(indexi, indexi);
			AP.addindex(indexi, indexi);
	
			// run through all intersections with neighbors 
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
			  	  is!=endit; ++is)
			    if (is->neighbor()) 
			      {
					// access neighbor
					EntityPointer outside = is->outside();
					int indexj = elementmapper.map(*outside);
		
					// add off diagonal index
					AV.addindex(indexi, indexj);
					AP.addindex(indexi, indexj);
			      }
	      }
	    AV.endindices();		
	    AP.endindices();		

	    return;
  }
  
  template<class G, class RT>
  void FVBrinkman<G, RT>::assembleMatrices()
  {
	  // initialization: set matrix A to zero	   
	  AV = 0;
	  AP = 0;

	  Iterator eendit = indexset.template end<0,All_Partition>();
	  for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	  {		
		  // cell geometry type
		  GeometryType gt = it->geometry().type();

		  // cell center in reference element
		  const FieldVector<ct,dim>& 
		  local = ReferenceElements<ct,dim>::general(gt).position(0,0);

		  // get global coordinate of cell center
		  FieldVector<ct,dim> global = it->geometry().global(local);

		  // cell index
		  int indexi = elementmapper.map(*it);

		  // cell volume 
		  double volume = it->geometry().integrationElement(local)
		  *ReferenceElements<ct,dim>::general(gt).volume();

		  // get absolute permeability 
		  FieldMatrix<ct,dim,dim> KinvI(this->problem.Kinv(global,*it,local));

		  // get effective viscosity
		  RT muEffI = this->problem.muEff(global,*it,local);

		  // get viscosity
		  RT muI = this->problem.mu(global,*it,local);

		  KinvI *= volume*muI;

		  AV[indexi][indexi] = KinvI;

		  IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
		  for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
		  is!=endit; ++is)
		  {

			  // get geometry type of face
			  GeometryType gtf = is->intersectionSelfLocal().type();

			  // center in face's reference element
			  const FieldVector<ct,dim-1>& 
			  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

			  // center of face inside volume reference element
			  const FieldVector<ct,dim>& 
			  facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is->numberInSelf(),1);

			  // get normal vector 
			  FieldVector<ct,dimworld> unitOuterNormal 
			  = is->unitOuterNormal(facelocal);

			  // get normal vector scaled with volume
			  FieldVector<ct,dimworld> integrationOuterNormal 
			  = is->integrationOuterNormal(facelocal);

			  // get face volume 
			  double faceVol = is->intersectionGlobal().volume();

			  // handle interior face
			  if (is->neighbor()) 
			  {
				  // access neighbor
				  EntityPointer outside = is->outside();
				  int indexj = elementmapper.map(*outside);

				  // compute factor in neighbor
				  GeometryType nbgt = outside->geometry().type();
				  const FieldVector<ct,dim>& 
				  nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);

				  // neighbor cell center in global coordinates
				  FieldVector<ct,dimworld> 
				  nbglobal = outside->geometry().global(nblocal);

				  // distance vector between barycenters
				  FieldVector<ct,dimworld> 
				  distVec = global - nbglobal;

				  // compute distance between cell centers
				  double dist = distVec.two_norm();

				  // get the effective viscosity 
				  RT muEffJ = this->problem.muEff(nbglobal, *outside, nblocal);

				  // average the effective viscosity 
				  RT muEff;
				  if (muEffI || muEffJ)
					  muEff = 2.0*muEffI*muEffJ/(muEffI + muEffJ);
				  else 
					  muEff = 0;

				  double entry = -muEff*(distVec*integrationOuterNormal)/(dist*dist);

				  for (int k = 0; k < dim; k++) {
					  // update diagonal entry 
					  AV[indexi][indexi][k][k] += entry;

					  // set off-diagonal entry 
					  AV[indexi][indexj][k][k] = -entry;
				  }
			  }
			  // boundary face 
			  else 
			  { 
				  // center of face in global coordinates
				  FieldVector<ct,dimworld> 
				  faceglobal = is->intersectionGlobal().global(facelocal);

				  //get boundary condition for boundary face center
				  BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);
				  if (bctype == BoundaryConditions::dirichlet) 
				  { 
					  FieldVector<ct,dimworld> distVec(global - faceglobal);
					  double dist = distVec.two_norm();

					  double entry = -muEffI*(distVec*integrationOuterNormal)/(dist*dist);

					  for (int k = 0; k < dim; k++) {
						  // update diagonal entry 
						  AV[indexi][indexi][k][k] += entry;
					  }

					  FieldVector<RT, dim> g = this->problem.gVelocity(faceglobal, *it, facelocalDim);
					  g *= entry;
					  fVBound[indexi] += g;
				  } 
				  else
				  {
					  FieldVector<RT, dim> J = this->problem.JVelocity(faceglobal, *it, facelocalDim);
					  J *= faceVol;
					  fVBound[indexi] += J;
				  }

			  }
		  } // end all intersections         
	  } // end grid traversal 

	  for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	  {		
		  int indexi = elementmapper.map(*it);
		  
		  double AVI = AV[indexi][indexi][0][0];

		  IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
		  for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
		  is!=endit; ++is)
		  {
			  // cell geometry type
			  GeometryType gt = it->geometry().type();

			  
			  // cell center in reference element
			  const FieldVector<ct,dim>& 
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0);

			  // get global coordinate of cell center
			  FieldVector<ct,dim> global = it->geometry().global(local);				  // compute factor in neighbor

			  // get geometry type of face
			  GeometryType gtf = is->intersectionSelfLocal().type();

			  // center in face's reference element
			  const FieldVector<ct,dim-1>& 
			  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

			  // center of face inside volume reference element
			  const FieldVector<ct,dim>& 
			  facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is->numberInSelf(),1);

			  if (is->neighbor()) 
			  {
				  EntityPointer outside = is->outside();
				  int indexj = elementmapper.map(*outside);

				  double AVJ = AV[indexj][indexj][0][0];

				  // changes


				  GeometryType nbgt = outside->geometry().type();
				  
				  const FieldVector<ct,dim>& 
				  nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);

				  // neighbor cell center in global coordinates
				  FieldVector<ct,dimworld> 
				  nbglobal = outside->geometry().global(nblocal);
	
				  // distance vector between barycenters
				  FieldVector<ct,dimworld> 
				  distVec = global - nbglobal;
	
				  // compute distance between cell centers
				  double dist = distVec.two_norm();
				  //CHANGES UNTIL HERE
				  
				  double averageEntry = 2.0/(AVI + AVJ)*dist*dist;
				  
				  AP[indexi][indexi] += averageEntry;
				  AP[indexi][indexj] -= averageEntry;
			  }
			  else 
			  {
				  // center of face in global coordinates
				  FieldVector<ct,dimworld> 
				  faceglobal = is->intersectionGlobal().global(facelocal);

				  //get boundary condition for boundary face center
				  BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);

				  // Neumann for velocity means Dirichlet for pressure correction and vice versa
				  if (bctype == BoundaryConditions::neumann) 
				  { 
					  FieldVector<ct,dimworld> distVec(global - faceglobal);
					  double dist = 2.0*distVec.two_norm();
					  AP[indexi][indexi] += 1.0/AVI*dist*dist;
				  }
			  }
		  }
	  }

	  return;
	}
	
	
  template<class G, class RT>
  void FVBrinkman<G, RT>::solveVelocitySystem()
  {
	  MatrixAdapter<VelocityMatrixType,VelType,VelType> op(AV); 
	  InverseOperatorResult r;
	  
	  SeqILU0<VelocityMatrixType,VelType,VelType> preconditioner(AV, 1.0);
	  BiCGSTABSolver<VelType> solver(op, preconditioner, 1E-14, 10000, 1);
	  solver.apply(this->velocity, fV, r);

	  return;
  }
	
  template<class G, class RT>
  void FVBrinkman<G, RT>::solvePressureSystem()
  {
	  MatrixAdapter<PressureMatrixType, RepresentationType,RepresentationType> op(AP); 
	  InverseOperatorResult r;
	  
	  SeqILU0<PressureMatrixType,RepresentationType,RepresentationType> preconditioner(AP, 1.0);
	  BiCGSTABSolver<RepresentationType> solver(op, preconditioner, 1E-14, 10000, 1);
	  solver.apply(this->pressureCorrection, fP, r);

	  return;
  }
	
  
}
#endif
