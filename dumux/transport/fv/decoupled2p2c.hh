#ifndef DECOUPLED2P2C_HH
#define DECOUPLED2P2C_HH
  
// commons:
#include <float.h>
#include "dumux/diffusion/diffusionproblem.hh"
#include <cmath>
#include <string>
#include <fstream>
#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/bvector.hh>

// transport:
#include "dumux/transport/transport.hh"
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem2p2c.hh"

// pressure:
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/diffusion/diffusion.hh"
#include "dumux/pardiso/pardiso.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/transport/problems/simpleproblem.hh"


// author: Jochen Fritz
// last change: 13.08.08

namespace Dune
{

/*####################################################*
 * 																										*
 *     CLASS DECLARATION															*
 * 																										*
 *####################################################*/

	template<class G, class RT>
	class Decoupled2p2c
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
	  typedef typename G::Traits::LevelIndexSet IS;
	  typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
	  typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	  typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	  typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
	  typedef typename G::ctype ct; 
	  typedef FieldMatrix<double,1,1> MB;
	  typedef BCRSMatrix<MB> MatrixType;
	  typedef FieldVector<double, 1> VB;
	  typedef BlockVector<VB> Vector;
	  typedef FieldVector<double,dim> R2;
	  typedef BlockVector<R2> R3;
	  typedef BlockVector<R3> LocVelType;
	
	public:  
		typedef BlockVector< Dune::FieldVector<RT,1> > RepresentationType;
		typedef BlockVector< Dune::FieldVector<Dune::FieldVector<RT, G::dimension>, 2*G::dimension> > VelocityType;
		
		void initial()
		{
			upd = 0;
			timestep = 0;
			initializeMatrix();
			initialguess();
			pressure(true, 0);
			transportInitial();
			pressure(false,0);
			transportInitial();
			totalVelocity(0);
		}
		
		void update(double t, double& dt, RepresentationType& updateVec)
		{
			upd = 0;
			concentrationUpdate(t, dt, upd);
			upd *= dt;
			timestep = dt;
			
			pressure(false, t);
			totalVelocity(t);
			concentrationUpdate(t, dt, updateVec);
		}
		
		int level()
		{
			return level_;
		}
		
		RepresentationType& operator*()
		{
			return problem.variables.totalConcentration;
		}
		
		// pressure equation functions:
	  void initializeMatrix();
	  
		void assemble(bool first, const RT t); 

		void solve(); 

		void pressure(bool first, const RT t=0)
		{
			assemble(first, t);
			solve();
			return;
		}

		void totalVelocity(const RT t);
	  
	  // transport equation functions:
	  void initialguess();
	  
	  void transportInitial(); 
	  
	  int concentrationUpdate(const RT t, RT& dt, RepresentationType& updateVec);
	  
	  void flashCalculation(double Z1, double p, double temp, double poro, double& sat, double& C1, double& C2, double& Cw1, double& Cn1, double& Cw2, double& Cn2);
	   
	  void satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Cw1, double& Cn1, double& Cw2, double& Cn2);
	   
	  void postupdate(double t, double dt);
	  
	  // graphical output
	  void vtkout(const char* name, int k)
	  {  
		  problem.variables.volErr *= timestep;
		  problem.variables.vtkout(name, k);
		  problem.variables.volErr /= timestep;
	  }
	  
		// constructor
		Decoupled2p2c(
				G& g, 
				TransportProblem2p2c<G, RT>& prob, 
				int lev = 0,
	   	  DiffusivePart<G, RT>& diffPart = *(new DiffusivePart<G, RT>),
	   	  bool rec = false, 
	   	  double amax = 0.8, 
	   	  const NumericalFlux<RT>& numFl = *(new Upwind<RT>),
	   	  const std::string solverName = "BiCGSTAB", 
			  const std::string preconditionerName = "SeqILU0" )
	   	  :	grid(g), level_(lev), indexset(g.levelIndexSet(lev)), reconstruct(rec), 
	   	  	numFlux(numFl), diffusivePart(diffPart), alphamax(amax), 
	   	  	problem(prob), 
	   	  	elementmapper(g, g.levelIndexSet(lev)),
	   	  	A(g.size(lev, 0),g.size(lev, 0), (2*dim+1)*g.size(lev, 0), BCRSMatrix<MB>::random), f(g.size(lev, 0)),
	   	  	solverName_(solverName), preconditionerName_(preconditionerName)
 	  {
			problem.variables.volErr = 0;
			upd.resize(elementmapper.size());
 	  };
		
	private:
		// common variables
		G& grid;
		int level_;
  	const IS& indexset;
  	EM elementmapper;
  	double timestep;
		
		// variables for transport equation:
		TransportProblem2p2c<G, RT>& problem;
		RepresentationType upd;
		
  	bool reconstruct;
  	const NumericalFlux<RT>& numFlux;
  	const DiffusivePart<G, RT>& diffusivePart;
  	double alphamax;
		
		// variables for pressure equation:
	  MatrixType A;
	  RepresentationType f;
	  std::string solverName_;
	  std::string preconditionerName_;
		
	}; //end class declaration
	
	
	
	/*####################################################*
	 * 																										*
	 *     FUNCTION DEFINITIONS 1: PRESSURE EQUATION			*
	 * 																										*
	 *####################################################*/
	
  template<class G, class RT>
  void Decoupled2p2c<G, RT>::initializeMatrix()
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
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); is!=endit; ++is)
		    if (is.neighbor()) rowSize++;
			A.setrowsize(indexi, rowSize);
    }
	  A.endrowsizes();

    // determine position of matrix entries 
    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
    {
			// cell index
			int indexi = elementmapper.map(*it);
	
			// add diagonal index
			A.addindex(indexi, indexi);
	
			// run through all intersections with neighbors 
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); is!=endit; ++is)
		    if (is.neighbor()) 
	      {
						// access neighbor
						EntityPointer outside = is.outside();
						int indexj = elementmapper.map(*outside);
			
						// add off diagonal index
						A.addindex(indexi, indexj);
	      }
    }
    A.endindices();		

	  return;
  } // end function initialize matrix
  
  
  /** for first == true, this function assembles the matrix and right hand side for 
   * the solution of the pressure field in the same way as in the class FVDiffusion. 
   * for first == false, the approach is changed to \f$ \[-\frac{\partial V}{\partial p}
   * \frac{\partial p}{\partial t}+\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}\nabla\cdot
   * \left(\sum_{\alpha}C_{\alpha}^{\kappa}\mathbf{v}_{\alpha}\right)
   * =\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}q^{\kappa}\] \f$. See Paper SPE 99619.
   * This is done to account for the volume effects which appear when gas and liquid are dissolved iin each other.
   */
  template<class G, class RT>
  void Decoupled2p2c<G, RT>::assemble(bool first, const RT t=0)
	{
	    // initialization: set matrix A to zero	   
        A = 0;
        
    // iterate over all cells
    Iterator eendit = indexset.template end<0,All_Partition>();
    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
    {	
	    // get geometry infos about the cell...
			GeometryType gt = it->geometry().type(); // cell geometry type
			const FieldVector<ct,dim>& 
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
			FieldVector<ct,dim> global = it->geometry().global(local); 			// global coordinate of cell center
			double volume = it->geometry().integrationElement(local)
				*ReferenceElements<ct,dim>::general(gt).volume(); // cell volume
				
			// cell index
			int indexi = elementmapper.map(*it);
	
			// get absolute permeability 
			FieldMatrix<ct,dim,dim> Ki(this->problem.K(global,*it,local));
	
			// get the cell's saturation
			double sati = problem.variables.saturation[indexi];

			// total mobility and fractional flow factors
			double lambdaI = problem.materialLaw.mobTotal(sati);
			double fw_I = problem.materialLaw.fractionalW(sati);
			double fn_I = problem.materialLaw.fractionalN(1-sati);
			
			// derivatives of the fluid volume with respect to mass of compnents and pressure
			double dV_dm1, dV_dm2, dV_dp;
			// specific volume of the phases
			double Vg, Vw;
			
			if (first)
			{
				// specific volume of the phases
				double Vg = 1. / problem.materialLaw.nonwettingPhase.density(283.15, 1e5);
				double Vw = 1. / problem.materialLaw.wettingPhase.density();
				
				FieldVector<RT,2> q = problem.q(global,*it,local);
				f[indexi] = volume * (q[0] * Vw + q[1] * Vg);
			}
			else
			{
				// specific volume of the phases
				double Vg = 1. / problem.materialLaw.nonwettingPhase.density(283.15, problem.variables.pressure[indexi]);
				double Vw = 1. / problem.materialLaw.wettingPhase.density();
				
				// mass of components inside the cell
				double m1 = problem.variables.totalConcentration[indexi]*volume; 
				double m2 = problem.variables.totalConcentration[indexi+elementmapper.size()]*volume;
				// mass fraction of wetting phase
				double nuw1 = sati / Vw / (sati/Vw + (1-sati)/Vg);
				// actual fluid volume
				double volalt = (m1+m2) * (nuw1 * Vw + (1-nuw1) * Vg);
				
				// increments for numerical derivatives
				double inc1 = (fabs(upd[indexi][0]) > 1e-8 /Vw) ?  upd[indexi][0] : 1e-8/Vw;
				double inc2 =(fabs(upd[indexi+elementmapper.size()][0]) > 1e-8 / Vg) ?  upd[indexi+elementmapper.size()][0] : 1e-8 / Vg;
				inc1 *= volume; 
				inc2 *= volume;
				
				// numerical derivative of fluid volume with respect to mass of component 1 
				m1 +=  inc1;
				double Z1 = m1 / (m1 + m2);
			  double dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, satt;
				flashCalculation(Z1, problem.variables.pressure[indexi], 283.15, problem.porosity(), satt, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6);
				double nuw = satt / Vw / (satt/Vw + (1-satt)/Vg);
				dV_dm1 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt) /inc1;
				m1 -= inc1;
				
				// numerical derivative of fluid volume with respect to mass of component 2 
				m2 += inc2;
				Z1 = m1 / (m1 + m2);
				flashCalculation(Z1, problem.variables.pressure[indexi], 283.15, problem.porosity(), satt, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6);
				nuw = satt / Vw / (satt/Vw + (1-satt)/Vg);
				dV_dm2 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt)/ inc2;
				m2 -= inc2;
				
				// numerical derivative of fluid volume with respect to pressure 
				double incp = 1e-5;
				double p_ = problem.variables.pressure[indexi] + incp;
				Vg = 1. / problem.materialLaw.nonwettingPhase.density(283.15, p_);
				Vw = 1. / problem.materialLaw.wettingPhase.density(283.15, p_);
				dV_dp = ((m1+m2) * (nuw1 * Vw + (1-nuw1) * Vg) - volalt) /incp;
				
				// right hand side entry: sources
				FieldVector<RT,2> q = problem.q(global,*it,local);
				f[indexi] = volume * (dV_dm1 * q[0] + dV_dm2 * q[1]);
			}

			// iterate over all faces of the cell
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
				  is!=endit; ++is)
			{
	    	// some geometry informations of the face
		    GeometryType gtf = is.intersectionSelfLocal().type(); // get geometry type of face
		    const FieldVector<ct,dim-1>& 
		      facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0); // center in face's reference element
		    const FieldVector<ct,dim>& 
		      facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1); // center of face inside volume reference element
		    FieldVector<ct,dimworld> unitOuterNormal = is.unitOuterNormal(facelocal);// get normal vector 
		    FieldVector<ct,dimworld> integrationOuterNormal = is.integrationOuterNormal(facelocal); 
		    integrationOuterNormal *= ReferenceElements<ct,dim-1>::general(gtf).volume(); //normal vector scaled with volume
		    double faceVol = is.intersectionGlobal().volume(); // get face volume 
			    
				// compute directed permeability vector Ki.n
				FieldVector<ct,dim> Kni(0);
				Ki.umv(unitOuterNormal, Kni);

				// handle interior face
			  if (is.neighbor()) 
			  {
					// acces neighbor
					EntityPointer outside = is.outside();
					int indexj = elementmapper.map(*outside);
					
					// some geometry infos of the neighbor
					GeometryType nbgt = outside->geometry().type();
					const FieldVector<ct,dim>& nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
					FieldVector<ct,dimworld> 
					  nbglobal = outside->geometry().global(nblocal); // neighbor cell center in global coordinates
		
					// distance vector between barycenters
					FieldVector<ct,dimworld> distVec = global - nbglobal;
		
					// compute distance between cell centers
					double dist = distVec.two_norm();
		
					// get absolute permeability 
					FieldMatrix<ct,dim,dim> Kj(problem.K(nbglobal, *outside, nblocal));
						
					// compute vectorized permeabilities
          FieldVector<ct,dim> Knj(0);
          Kj.umv(unitOuterNormal, Knj);
          double K_n_i = Kni * unitOuterNormal;
          double K_n_j = Knj * unitOuterNormal;
          double Kn    = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
          // compute permeability tangential to intersection and take arithmetic mean
          FieldVector<ct,dim> uON = unitOuterNormal;
          FieldVector<ct,dim> K_t_i = Kni - (uON *= K_n_i);
          uON = unitOuterNormal;
          FieldVector<ct,dim> K_t_j = Knj - (uON *= K_n_j);
          FieldVector<ct,dim> Kt = (K_t_i += K_t_j);
          Kt *= 0.5;
          // Build vectorized averaged permeability
          uON = unitOuterNormal;
          FieldVector<ct,dim> K = (Kt += (uON *=Kn));
	
					//compute mobilities
					double lambdaJ;
					double fw_J, fn_J;
					double satj = problem.variables.saturation[indexj];
					lambdaJ = problem.materialLaw.mobTotal(satj);
					
					if (!first)
					{
						fw_J = problem.materialLaw.fractionalW(satj);
						fn_J = problem.materialLaw.fractionalN(1-satj);
					}
				
					// compute averaged total mobility
					// CAREFUL: Harmonic weightig can generate zero matrix entries, 
			    // use arithmetic weighting instead: 
					double lambda;
						lambda = 0.5*(lambdaI + lambdaJ);
					
					// update diagonal entry 
					double entry;
					if (first)
						entry = fabs(lambda*faceVol*(K*distVec)/(dist*dist));
					else
					{
						if (problem.variables.pressure[indexi] > problem.variables.pressure[indexj])
						{
							entry = fabs( 
									     dV_dm1 * ( problem.variables.wet_c1[indexi] * fw_I + problem.variables.nonwet_c1[indexi] * fn_I) 
										 + dV_dm2 * ( problem.variables.wet_c2[indexi] * fw_I + problem.variables.nonwet_c2[indexi] * fn_I)
											 );
						}
						else
						{
							entry = fabs( 
									     dV_dm1 * ( problem.variables.wet_c1[indexj] * fw_J + problem.variables.nonwet_c1[indexj] * fn_J) 
										 + dV_dm2 * ( problem.variables.wet_c2[indexj] * fw_J + problem.variables.nonwet_c2[indexj] * fn_J)
											 );
						}	
						entry *= lambda * fabs(faceVol*(K*distVec)/(dist*dist));
					}
					// set diagonal entry
					A[indexi][indexi] += entry;
		
					// set off-diagonal entry 
					A[indexi][indexj] = -entry;
			  }
			  
			  // boundary face 
			  else 
			  { 
					// center of face in global coordinates
					FieldVector<ct,dimworld> 
					  faceglobal = is.intersectionGlobal().global(facelocal);
					  
					// compute total mobility
					double lambda = lambdaI;
		
					//get boundary condition for boundary face center
					BoundaryConditions::Flags bctype = problem.pbctype(faceglobal, *it, facelocalDim);
					
					//dirichlet boundary
					if (bctype == BoundaryConditions::dirichlet) 
					{ 
						// distance vector to boundary face
						FieldVector<ct,dimworld> distVec(global - faceglobal);
						double dist = distVec.two_norm();
						if (first)
						{	
							A[indexi][indexi] -= lambda * faceVol * (Kni * distVec) / (dist * dist);
							double pressBC = problem.gPress(faceglobal, *it, facelocalDim);
							f[indexi] -= lambda * faceVol * pressBC * (Kni * distVec) / (dist * dist);
						}
						else
						{
							double pressBC = problem.gPress(faceglobal, *it, facelocalDim);
							
							double satBound, C1Bound, C2Bound, cw1Bound, cn1Bound, cw2Bound, cn2Bound;
							
							//get boundary condition type for compositional transport
							BoundaryConditions2p2c::Flags bctype = problem.cbctype(faceglobal, *it, facelocalDim);
							if (bctype == BoundaryConditions2p2c::saturation) // saturation given
							{
								satBound = problem.gS(faceglobal, *it, facelocalDim);
								satFlash(satBound, pressBC, 283.15, problem.porosity(), C1Bound, C2Bound, cw1Bound, cn1Bound, cw2Bound, cn2Bound);
							}
							if (bctype == BoundaryConditions2p2c::concentration) // mass fraction given
							{
								double Z1Bound = problem.gZ(faceglobal, *it, facelocalDim);
								double dummy;
								flashCalculation(Z1Bound, pressBC, 283.15, problem.porosity(), satBound, C1Bound, C2Bound, cw1Bound, cn1Bound, cw2Bound, cn2Bound);
							}
							
							double fw_J = problem.materialLaw.fractionalW(satBound);
							double fn_J = problem.materialLaw.fractionalN(1-satBound);
							
							double entry;
							if (problem.variables.pressure[indexi] > pressBC)
								entry = fabs( 
							     dV_dm1 * ( problem.variables.wet_c1[indexi] * fw_I + problem.variables.nonwet_c1[indexi] * fn_I) 
									 + dV_dm2 * ( problem.variables.wet_c2[indexi] * fw_I + problem.variables.nonwet_c2[indexi] * fn_I)
									 );
							else
								entry = fabs( 
							     dV_dm1 * ( cw1Bound * fw_I + cn1Bound * fn_I) 
									 + dV_dm2 * ( cw2Bound * fw_I + cn2Bound * fn_I)
									 );
								
							entry *= - lambda * faceVol*(Kni*distVec)/(dist*dist);
							
							// set diagonal entry and right hand side entry
							A[indexi][indexi] += entry;
							f[indexi] += entry * pressBC;
						}
					} 
					else
					{
						FieldVector<RT,2> J = problem.J(faceglobal, *it, facelocalDim);
						if (first)
							f[indexi] += faceVol * (J[0] * Vw + J[1] * Vg);
						else
						{
							f[indexi] += faceVol* (J[0] * dV_dm1 + J[1] * dV_dm2);
						}
					}
			  }
			} // end all intersections
			
			// compressibility term
			if (!first && timestep != 0) 
			{
				A[indexi][indexi] -= dV_dp / timestep;
				f[indexi] -= problem.variables.pressure[indexi] *dV_dp / timestep;
			}
			
			// error reduction routine: volumetric error is inserted to right hand side
			double maxErr = fabs(problem.variables.volErr.infinity_norm());
			double erri = fabs(problem.variables.volErr[indexi]);
			double x_lo = 0.6;
			double x_mi = 0.9;
			double fac  = 0.5;
			double lofac = 0.;
			double hifac = 0.; 
			hifac /= fac;
			
			if (erri*timestep > 5e-4)
			if (erri > x_lo * maxErr) 
			{
				if (erri <= x_mi * maxErr)
					f[indexi] += fac* (1-x_mi*(lofac/fac-1)/(x_lo-x_mi) + (lofac/fac-1)/(x_lo-x_mi)*erri/maxErr)* problem.variables.volErr[indexi]*volume;
				else 
					f[indexi] += fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxErr) * problem.variables.volErr[indexi]*volume;
			}
	  } // end grid traversal 
	    
//	    printmatrix(std::cout,A,"stiffnesmatrix","row");
//	    printvector(std::cout,f,"right hand side","row");
	    
	    return;
	} // end function assemble
	
	
  template<class G, class RT>
  void Decoupled2p2c<G, RT>::solve()
  {
	  MatrixAdapter<MatrixType,Vector,Vector> op(A); 
	  InverseOperatorResult r;
	  
	  if (preconditionerName_ == "SeqILU0") {
	      SeqILU0<MatrixType,Vector,Vector> preconditioner(A, 1.0);
	      if (solverName_ == "CG") {
	    	  CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
	    	  solver.apply(problem.variables.pressure, f, r);
	      }
	      else if (solverName_ == "BiCGSTAB") {
	    	  BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
	    	  solver.apply(problem.variables.pressure, f, r);
	      }
	      else 
			  DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination " << preconditionerName_ 
					  << " and " << solverName_ << ".");
	  }
	  else if (preconditionerName_ == "SeqPardiso") {
	      SeqPardiso<MatrixType,Vector,Vector> preconditioner(A);
	      if (solverName_ == "Loop") {
	    	  LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
	    	  solver.apply(problem.variables.pressure, f, r);
	      }
	      else 
	    	  DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination " << preconditionerName_ 
	    			  << " and " << solverName_ << ".");
	  }
	  else 
		  DUNE_THROW(NotImplemented, "FVDiffusion :: solve : preconditioner " << preconditionerName_ << ".");

	  return;
  }
	
  template<class G, class RT>
  void Decoupled2p2c<G, RT>::totalVelocity(const RT t=0)
	{
      // find out whether gravity effects are relevant
      bool hasGravity = false;
      const FieldVector<ct,dim>& gravity = problem.gravity();
      for (int k = 0; k < dim; k++)
      	if (gravity[k] != 0) 
      		hasGravity = true;
      
	    Iterator eendit = indexset.template end<0,All_Partition>();
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	    {
	    	// cell geometry type
	    	GeometryType gt = it->geometry().type();
		      
	    	// cell center in reference element
	    	const FieldVector<ct,dim>& 
				local = ReferenceElements<ct,dim>::general(gt).position(0,0);
		      
			// cell center in global coordinates
			FieldVector<ct,dimworld> 
				global = it->geometry().global(local);
		      
		    // cell index
			int indexi = elementmapper.map(*it);
		    
			// get pressure and permeability in element
			double pressi = this->problem.variables.pressure[indexi];
		      
			// get absolute permeability 
			FieldMatrix<ct,dim,dim> Ki(problem.K(global,*it,local));
	
			//compute total mobility
			double lambdaI, fractionalWI;
			double sati = problem.variables.saturation[indexi];
			lambdaI = problem.materialLaw.mobTotal(sati);
			if (hasGravity) 
				fractionalWI = problem.materialLaw.fractionalW(sati);
			
			double faceVol[2*dim];

			// run through all intersections with neighbors and boundary
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
				  is!=endit; ++is)
			{
			  // get geometry type of face
			  GeometryType gtf = is.intersectionSelfLocal().type();

			  //Geometry dg = is.intersectionSelfLocal();
			  // local number of facet 
			  int numberInSelf = is.numberInSelf();

			  faceVol[numberInSelf] = is.intersectionGlobal().volume();
			    
			  // center in face's reference element
			  const FieldVector<ct,dim-1>& 
			    facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);
			  
			  // center of face inside volume reference element
			  const FieldVector<ct,dim>& 
			  	facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
				    
			  // get normal vector
			  FieldVector<ct,dimworld> unitOuterNormal 
			    = is.unitOuterNormal(facelocal);
			  
			  // center of face in global coordinates
			  FieldVector<ct,dimworld> 
			    faceglobal = is.intersectionGlobal().global(facelocal);
			  
			  // handle interior face
			  if (is.neighbor()) 
			    {
			      // access neighbor
			      EntityPointer outside = is.outside();
			      int indexj = elementmapper.map(*outside);
			      
			      // get neighbor pressure and permeability
			      double pressj = this->problem.variables.pressure[indexj];
			      
			      // compute factor in neighbor
			      GeometryType nbgt = outside->geometry().type();
			      const FieldVector<ct,dim>& 
			      	nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
			      
			      // neighbor cell center in global coordinates
			      FieldVector<ct,dimworld> 
			      	nbglobal = outside->geometry().global(nblocal);
			      
			      // distance vector between barycenters
			      FieldVector<ct,dimworld> distVec = global - nbglobal;
			      
			      // compute distance between cell centers
			      double dist = distVec.two_norm();
			      
			      // get absolute permeability 
			      FieldMatrix<ct,dim,dim> Kj(problem.K(nbglobal, *outside, nblocal));

			      // compute vectorized permeabilities
                  FieldVector<ct,dim> Kni(0);
                  FieldVector<ct,dim> Knj(0);
                  Ki.umv(unitOuterNormal, Kni);
                  Kj.umv(unitOuterNormal, Knj);
                  // compute permeability normal to intersection and take harmonic mean
                  double K_n_i = Kni * unitOuterNormal;
                  double K_n_j = Knj * unitOuterNormal;
                  double Kn    = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
                  // compute permeability tangential to intersection and take arithmetic mean
                  FieldVector<ct,dim> uON = unitOuterNormal;
                  FieldVector<ct,dim> K_t_i = Kni - (uON *= K_n_i);
                  uON = unitOuterNormal;
                  FieldVector<ct,dim> K_t_j = Knj - (uON *= K_n_j);
                  FieldVector<ct,dim> Kt = (K_t_i += K_t_j);
                  Kt *= 0.5;
                  // Build vectorized averaged permeability
                  uON = unitOuterNormal;
                  FieldVector<ct,dim> K = (Kt += (uON *=Kn));

			      //compute total mobility
			      double lambdaJ, fractionalWJ;
			      double satj = problem.variables.saturation[indexj];
			    	  lambdaJ = problem.materialLaw.mobTotal(satj);
			    	  if (hasGravity) 	
			    		  fractionalWJ = problem.materialLaw.fractionalW(satj);
			      
			      // compute averaged total mobility
			      // CAREFUL: Harmonic weightig can generate zero matrix entries, 
			      // use arithmetic weighting instead: 
			      double lambda = 1;
			      double fractionalW;
			      lambda = 0.5*(lambdaI + lambdaJ); 
			      if (hasGravity) 
			    	  fractionalW = 0.5*(fractionalWI + fractionalWJ);
			      
			      FieldVector<ct,dimworld> vTotal(K);
			      vTotal *= lambda*(pressi - pressj)/dist;
			      problem.variables.velocity[indexi][numberInSelf] = vTotal;	
			    }
			  // boundary face 
			  else 
		    { 
		      //get boundary condition for boundary face center
			  BoundaryConditions::Flags bctype = problem.pbctype(faceglobal, *it, facelocalDim);
			  if (bctype == BoundaryConditions::dirichlet) 
			  {
			      // distance vector between barycenters
			      FieldVector<ct,dimworld> distVec = global - faceglobal;
			      
				  double dist = distVec.two_norm();
				  distVec /= dist;
				  
			      // compute directed permeability vector Ki.n
				  FieldVector<ct,dim> Kni(0);
				  Ki.umv(distVec, Kni);

				  // compute averaged total mobility
				  double lambda = 1.;
				  double fractionalW = 1.;
					  lambda = lambdaI;
					  if (hasGravity) fractionalW = fractionalWI;
					  
				  double pressBC = problem.gPress(faceglobal, *it, facelocalDim);
				  
				  FieldVector<ct,dim> vTotal(Kni);
				  vTotal *= lambda * (pressBC - pressi) / dist;
			      if (hasGravity) {
			    	  FieldVector<ct,dimworld> gEffect(0);
			    	  Ki.umv(gravity, gEffect);
			    	  double factor = fractionalW*(problem.materialLaw.wettingPhase.density()) 
							+ (1 - fractionalW)*(problem.materialLaw.nonwettingPhase.density());
			    	  gEffect *= lambda*factor;
			    	  vTotal -= gEffect;
			      }
			      problem.variables.velocity[indexi][numberInSelf] = vTotal;
		      } 
		      else
		      {		   
		      	// for Neumann boundary, only mass inflow but no volumetric flow is known.
		      	// To avoid the use of wrong values, set velocity to nonsens value. 
		      	problem.variables.velocity[indexi][numberInSelf] = DBL_MAX; // highest number in double precission.
		      }   
		    }
			} // end all intersections  
    } // end grid traversal          
	    
	  return;
	} // end function totalVelocity
  
	
	/*####################################################*
	 * 																										*
	 *     FUNCTION DEFINITIONS 2: TRANSPORT EQUATION			*
	 * 																										*
	 *####################################################*/
  
  template<class G, class RT>
  void Decoupled2p2c<G,RT>::initialguess()
  {
		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = indexset.template end<0,All_Partition>();
		for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		{
			int indexi = elementmapper.map(*it);
			
			// get geometry information of cell
			GeometryType gt = it->geometry().type();
			const FieldVector<ct,dim>& 
				local = ReferenceElements<ct,dim>::general(gt).position(0,0);
			FieldVector<ct,dimworld> global = it->geometry().global(local);
			
			double sat_0, Z1_0;
			// get type of initial condition
			BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *it, local);
			// saturation initial condition
			if (ictype == BoundaryConditions2p2c::saturation)
				sat_0 = problem.S0(global, *it, local);
			// concentration initial condition
			else if(ictype == BoundaryConditions2p2c::concentration)
			{
				Z1_0 = problem.Z1_0(global, *it, local);
				sat_0 = Z1_0 / problem.materialLaw.wettingPhase.density(283.15, 1e5);
				sat_0 /= Z1_0 / problem.materialLaw.wettingPhase.density(283.15, 1e5) + (1 - Z1_0) * problem.materialLaw.nonwettingPhase.density(283.15, 1e5);
			}
			// initialize cell saturation
			this->problem.variables.saturation[indexi][0] = sat_0;
		}	
	return;
  }//end function initialguess
  
  
  
  template<class G, class RT>
  void Decoupled2p2c<G,RT>::transportInitial() 
	{ 
		// iterate through leaf grid an evaluate c0 at cell center
	    Iterator eendit = indexset.template end<0,All_Partition>();
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		{
	    	int indexi = elementmapper.map(*it);
	    	
			// get geometry information of cell
			Dune::GeometryType gt = it->geometry().type();
			const Dune::FieldVector<ct,dim>& 
				local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
			Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);
			
			double sat_0, C1_0, C2_0;
			// get type of initial condition
			Dune::BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *it, local);
			// saturation initial condition
			if (ictype == Dune::BoundaryConditions2p2c::saturation)
			{
				sat_0 = problem.S0(global, *it, local);
				satFlash(sat_0, problem.variables.pressure[indexi], 283.15, problem.porosity(), C1_0, C2_0, problem.variables.wet_c1[indexi][0], problem.variables.nonwet_c1[indexi][0], problem.variables.wet_c2[indexi][0], problem.variables.nonwet_c2[indexi][0]);
			}
			// concentration initial condition
			if (ictype == Dune::BoundaryConditions2p2c::concentration)
			{
				double Z1_0 = problem.Z1_0(global, *it, local);
				double dummy;
				flashCalculation(Z1_0, problem.variables.pressure[indexi], 283.15, problem.porosity(), sat_0, C1_0, C2_0, problem.variables.wet_c1[indexi][0], problem.variables.nonwet_c1[indexi][0], problem.variables.wet_c2[indexi][0], problem.variables.nonwet_c2[indexi][0]);
			}
			// initialize cell concentration
			problem.variables.totalConcentration[indexi] = C1_0;
			problem.variables.totalConcentration[indexi + elementmapper.size()] = C2_0;
			this->problem.variables.saturation[indexi][0] = sat_0;
		}
		return;
	} //end function transportInitial
  
  
  template<class G, class RT>
  int Decoupled2p2c<G,RT>::concentrationUpdate(const RT t, RT& dt, RepresentationType& updateVec)
  {
		// initialize dt very large
		dt = 1E100;

		// set update vector to zero
		updateVec = 0;

		// compute update vector 
		Iterator eendit = indexset.template end<0,All_Partition>();
		for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		{
			// cell geometry type
			Dune::GeometryType gt = it->geometry().type();

			// cell center in reference element
			const Dune::FieldVector<ct,dim>& 
				local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
			
			// cell center in global coordinates
			FieldVector<ct,dimworld> global = it->geometry().global(local);

			// cell volume, assume linear map here
			double volume = it->geometry().integrationElement(local)
				*Dune::ReferenceElements<ct,dim>::general(gt).volume();

			// cell index
			int indexi = elementmapper.map(*it);
			
			// compute source term
			updateVec[indexi] += problem.q(global, *it, local)[0] * volume;
			updateVec[indexi+elementmapper.size()] += problem.q(global, *it, local)[1] * volume;
			
	    // get saturation and concentration value at cell center
	    double satI = this->problem.variables.saturation[indexi];
	    double cw1_I = problem.variables.wet_c1[indexi];
	    double cn1_I = problem.variables.nonwet_c1[indexi];
	    double cw2_I = problem.variables.wet_c2[indexi];
	    double cn2_I = problem.variables.nonwet_c2[indexi];

			// for time step calculation
			double sumfactor = 0;
			double sumfactor2 = 0;
			double sumDiff = 0;
			double sumDiff2 = 0;

			// run through all intersections with neighbors and boundary
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
				  is!=endit; ++is)
			{
				// local number of facet 
				int numberInSelf = is.numberInSelf();

				// get geometry type of face
				Dune::GeometryType gtf = is.intersectionSelfLocal().type();
				  
				// center in face's reference element
				const Dune::FieldVector<ct,dim-1>& 
					facelocal = Dune::ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

				// center of face inside volume reference element
				const Dune::FieldVector<ct,dim>& 
					facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1);
					    
				// get normal vector scaled with volume of face
				Dune::FieldVector<ct,dimworld> integrationOuterNormal 
					= is.integrationOuterNormal(facelocal);
				integrationOuterNormal 
					*= Dune::ReferenceElements<ct,dim-1>::general(gtf).volume();

				// compute factor occuring in flux formula	
				double velocityIJ = std::max(problem.variables.velocity[indexi][numberInSelf]*integrationOuterNormal/(volume), 0.0);
				double factor, diffFactor, factorC1, factorC2;

				// handle interior face
				if (is.neighbor())
				{
					// access neighbor
					EntityPointer outside = is.outside();
					int indexj = elementmapper.map(*outside);
					  
					if ( true )
					{
						// compute factor in neighbor
						Dune::GeometryType nbgt = outside->geometry().type();
						const Dune::FieldVector<ct,dim>& 
							nblocal = Dune::ReferenceElements<ct,dim>::general(nbgt).position(0,0);

					    double velocityJI = std::max(-(problem.variables.velocity[indexi][numberInSelf]*integrationOuterNormal/volume), 0.0);
					      
					    // cell center in global coordinates
					    Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);
					      
					    // neighbor cell center in global coordinates
					    Dune::FieldVector<ct,dimworld> nbglobal = outside->geometry().global(nblocal);
					      
					    // distance vector between barycenters
					    Dune::FieldVector<ct,dimworld> distVec = global - nbglobal;
					      
					    // compute distance between cell centers
					    double dist = distVec.two_norm();
					    
					    // get saturation and concentrtion value at neighbor cell center
					    double satJ = this->problem.variables.saturation[indexj];
					    double cw1_J = problem.variables.wet_c1[indexj];
					    double cn1_J = problem.variables.nonwet_c1[indexj];
					    double cw2_J = problem.variables.wet_c2[indexj];
					    double cn2_J = problem.variables.nonwet_c2[indexj];
						
					    // calculate the cocentration gradients 
					    Dune::FieldVector<ct,dim> cn1Gradient = distVec;		
					    cn1Gradient *= (cn1_J - cn1_I)/(dist*dist);
					    Dune::FieldVector<ct,dim> cw1Gradient = distVec;		
					    cw1Gradient *= (cw1_J - cw1_I)/(dist*dist);
					    
					    // the arithmetic average 
					    double satAvg = 0.5*(satI + satJ);
					    
					    // get the diffusive part
					    double n1_diffPart = 0;//this->diffusivePart(*it, numberInSelf, satAvg, cn1Gradient, t)*integrationOuterNormal;
					    double w1_diffPart = 0;//this->diffusivePart(*it, numberInSelf, satAvg, cw1Gradient, t)*integrationOuterNormal;

//					    // CAREFUL: works only for axisymmetric grids 
//					    if (reconstruct) {
//					    	for (int k = 0; k < dim; k++) 
//					    		if (fabs(distVec[k]) > 0.5*dist) 
//					    		{
//					    			satI -= fabs(distVec[k])/distVec[k]*0.5*dist;//*slope[indexi][k];
//					    			satJ += fabs(distVec[k])/distVec[k]*0.5*dist;//*slope[indexj][k];
//					    		}
//					    }
					    
			 		    double fwI = problem.materialLaw.fractionalW(satI);
					    double fwJ = problem.materialLaw.fractionalW(satJ);
					    double fnI = problem.materialLaw.fractionalN(1.0-satI);
					    double fnJ = problem.materialLaw.fractionalN(1.0-satJ);
					    					    
					    // for timestep control
					    {
					    	double coeffW = isnan(fwI / satI) ? 0 : fwI / satI;
					    	double coeffN = isnan(fnI / (1-satI)) ? 0 : fnI / (1-satI);
					    	factor = (velocityJI - velocityIJ) * std::max(std::max(coeffW, coeffN),1.);
					    }
					    
					  //  diffFactor =diffPart / volume; TODO include diffusion into timestep control
					    factorC1 = (n1_diffPart + w1_diffPart) / volume
					    		+ velocityJI * cw1_J * fwJ /*numFlux(satJ, satI, fwJ, fwI)*/
					    		- velocityIJ * cw1_I * fwI /*numFlux(satI, satJ, fwI, fwJ)*/
					    		+ velocityJI * cn1_J * fnJ /*numFlux(1.0-satJ, 1.0-satI, fnJ, fnI)*/
					    		- velocityIJ * cn1_I * fnI /*numFlux(1.0-satI, 1.0-satJ, fnI, fnJ)*/; 
					    factorC2 = 
					    		 velocityJI * cw2_J * fwJ /*numFlux(satJ, satI, fwJ, fwI)*/
					    		- velocityIJ * cw2_I * fwI /*numFlux(satI, satJ, fwI, fwJ)*/
					    		+ velocityJI * cn2_J * fnJ /*numFlux(1.0-satJ, 1.0-satI, fnJ, fnI)*/
					    		- velocityIJ * cn2_I * fnI /*numFlux(1.0-satI, 1.0-satJ, fnI, fnJ)*/; 
					    
					    Dune::FieldVector<ct,dimworld> faceglobal = is.intersectionGlobal().global(facelocal);
					}
				}
				  
				// handle boundary face
				else
				{
					// cell center in global coordinates
					Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);
				    
					Dune::FieldVector<ct,dim> faceglobal = is.intersectionGlobal().global(facelocal);
					
					// distance vector between barycenters
					Dune::FieldVector<ct,dimworld> distVec = global - faceglobal;
					
					// get saturation value at cell center
					double satI = this->problem.variables.saturation[indexi];

					double velocityJI = std::max(-(problem.variables.velocity[indexi][numberInSelf]*integrationOuterNormal/volume), 0.0);
					
					//get boundary conditions
					BoundaryConditions::Flags pressBCtype = problem.pbctype(faceglobal, *it, facelocalDim);
					if (pressBCtype == BoundaryConditions::dirichlet)
					{
						double satBound, C1Bound, C2Bound, cw1Bound, cn1Bound, cw2Bound, cn2Bound;
						BoundaryConditions2p2c::Flags bctype = problem.cbctype(faceglobal, *it, facelocalDim);
						if (bctype == BoundaryConditions2p2c::saturation)
						{
							satBound = problem.gS(faceglobal, *it, facelocalDim);
							satFlash(satBound, problem.gPress(faceglobal, *it, facelocalDim), 283.15, problem.porosity(), C1Bound, C2Bound, cw1Bound, cn1Bound, cw2Bound, cn2Bound);
						}
						if (bctype == BoundaryConditions2p2c::concentration)
						{
							double Z1Bound = problem.gZ(faceglobal, *it, facelocalDim);
							double dummy;
							flashCalculation(Z1Bound, problem.gPress(faceglobal, *it, facelocalDim), 283.15, problem.porosity(), satBound, C1Bound, C2Bound, cw1Bound, cn1Bound, cw2Bound, cn2Bound);
						}
					
						double dist = distVec.two_norm();
					    
						// CAREFUL: works only for axisymmetric grids 
						if (reconstruct) {
							for (int k = 0; k < dim; k++) 
								if (fabs(distVec[k]) > 0.5*dist) 
								{
									satI -= fabs(distVec[k])/distVec[k]*dist;//*slope[indexi][k];
								}
						}
					
			    	double cw1_I = problem.variables.wet_c1[indexi];
			    	double cn1_I = problem.variables.nonwet_c1[indexi];
			    	double cw2_I = problem.variables.wet_c2[indexi];
			    	double cn2_I = problem.variables.nonwet_c2[indexi];
			    				    	
			    	double fwI     = problem.materialLaw.fractionalW(satI);
			    	double fwBound = problem.materialLaw.fractionalW(satBound);
			    	double fnI     = problem.materialLaw.fractionalN(1.-satI);
			    	double fnBound = problem.materialLaw.fractionalN(1.-satBound);
			    	
				    // for timestep control
				    {
				    	double coeffW = isnan(fwI / satI) ? 0 : fwI / satI;
				    	double coeffN = isnan(fnI / (1-satI)) ? 0 : fnI / (1-satI);
				    	factor = (velocityJI - velocityIJ) * std::max(std::max(coeffW, coeffN),1.);
				    }
			    	
			    	factorC1 =
							+ velocityJI * cw1Bound * numFlux(satBound, satI, fwBound, fwI)
							- velocityIJ * cw1_I * numFlux(satI, satBound, fwI, fwBound)
							+ velocityJI * cn1Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
							- velocityIJ * cn1_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound); 
			    	factorC2 = 
							 velocityJI * cw2Bound * numFlux(satBound, satI, fwBound, fwI)
							- velocityIJ * cw2_I * numFlux(satI, satBound, fwI, fwBound)
							+ velocityJI * cn2Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
							- velocityIJ * cn2_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound);
					}
					else if (pressBCtype == BoundaryConditions::neumann)
					{
						FieldVector<RT,2> J = problem.J(faceglobal, *it, facelocalDim);
						double faceVol = integrationOuterNormal.two_norm();
						factorC1 = J[0] * faceVol;
						factorC2 = J[1] * faceVol;
						
						double fnI     = problem.materialLaw.fractionalN(1.-satI);
						double fwI     = problem.materialLaw.fractionalW(satI);
						
				    // for timestep control
				    {
				    	double coeffW = satI==0 ? 0 : fwI / satI;
				    	double coeffN = satI==1 ? 0 : fnI / (1-satI);
				    	factor = fabs(J[0] * problem.materialLaw.wettingPhase.density(283.15, problem.variables.pressure[indexi]) + J[1] * problem.materialLaw.nonwettingPhase.density(283.15, problem.variables.pressure[indexi]));
				    	factor *= std::max(std::max(coeffW, coeffN),1.);
				    }
					}
					else DUNE_THROW(Dune::NotImplemented, "there is no process boundary condition implemented");
				}

				// add to update vector 
				updateVec[indexi] += factorC1;
				updateVec[elementmapper.size()+indexi] += factorC2;

				// for time step calculation
				if (factor>=0) 
					sumfactor += factor;
				else 
					sumfactor2 += (-factor);
				if (diffFactor>=0)
					sumDiff += diffFactor;
				else
					sumDiff += (-diffFactor);
			} // end all intersections    
			// compute dt restriction
//			volInc[indexi] = sumfactor - sumfactor2;
			
			sumfactor = std::max(sumfactor,sumfactor2);
			sumDiff = std::max(sumDiff,sumDiff2);
			sumfactor = std::max(sumfactor,100*sumDiff);
			dt = std::min(dt,1.0/sumfactor);     
						
		} // end grid traversal                 
		dt *= problem.porosity();
		updateVec /= problem.porosity();
		return 0;
  } // end function update
  
  template<class G, class RT>
  void Decoupled2p2c<G,RT>::flashCalculation(double Z1, double p, double temp, double poro, double& sat, double& C1, double& C2, double& Cw1, double& Cn1, double& Cw2, double& Cn2)
  {
	  double K1 = problem.materialLaw.wettingPhase.vaporPressure(temp) / p;
    double K2 = 1. / (p * problem.henry(temp));
    double xw1 = (1. - K2) / (K1 -K2);
    double xn1 = xw1 * K1;
    double Xw1 = xw1 * problem.materialLaw.wettingPhase.molarMass() 
                / ( xw1 * problem.materialLaw.wettingPhase.molarMass() + (1.-xw1) * problem.materialLaw.nonwettingPhase.molarMass() );
    double Xn1 = xn1 * problem.materialLaw.wettingPhase.molarMass() 
                / ( xn1 * problem.materialLaw.wettingPhase.molarMass() + (1.-xn1) * problem.materialLaw.nonwettingPhase.molarMass() );
    K1 = Xn1 / Xw1;
    K2 = (1.-Xn1) / (1.-Xw1);
    double Xw2 = 1- Xw1;
    double Xn2 = 1- Xn1;
    double rho_w = problem.materialLaw.wettingPhase.density(temp, p);
    double rho_n = problem.materialLaw.nonwettingPhase.density(temp, p);
    
    double nu2;
    
    if (Z1 > Xn1 && Z1 < Xw1) 
    	nu2 = -((K1-1)*Z1 + (K2-1)*(1-Z1)) / (K1-1) / (K2-1);
    
    else if (Z1 < Xn1)
    {
    	nu2 = 1;
    	Xn1 = Z1;
    	Xn2 = 1 - Z1;
    }
    else if (Z1 > Xw1)
    {
    	nu2 = 0;
    	Xw1 = Z1;
    	Xw2 = 1 - Z1;
		}
		         
    sat = (1-nu2) / rho_w;
    sat /= ((1-nu2)/rho_w + nu2/rho_n);
    
    Cw1 = rho_w * Xw1;
    Cw2 = rho_w * Xw2;
    Cn1 = rho_n * Xn1;
    Cn2 = rho_n * Xn2;
    C1 = Cw1 * sat + Cn1 * (1-sat);
    C2 = Cw2 * sat + Cn2 * (1-sat);
    
  } // end function flashCalculation
  
  template<class G, class RT>
  void Decoupled2p2c<G,RT>::satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Cw1, double& Cn1, double& Cw2, double& Cn2)
  {
  	if (sat <= 0 || sat >= 1)
	  DUNE_THROW(RangeError,
	  	"FVTransport2p2c :: saturation initial and boundary conditions may not equal zero or one!");
    double K1 = problem.materialLaw.wettingPhase.vaporPressure(temp) / p;
    double K2 = 1. / (p * problem.henry(temp));
    double xw1 = (1. - K2) / (K1 -K2);
    double xn1 = xw1 * K1;
    double Xw1 = xw1 * problem.materialLaw.wettingPhase.molarMass() 
                / ( xw1 * problem.materialLaw.wettingPhase.molarMass() + (1.-xw1) * problem.materialLaw.nonwettingPhase.molarMass() );
    double Xn1 = xn1 * problem.materialLaw.nonwettingPhase.molarMass() 
                / ( xn1 * problem.materialLaw.nonwettingPhase.molarMass() + (1.-xn1) * problem.materialLaw.wettingPhase.molarMass() );
    K1 = Xn1 / Xw1;
    K2 = (1.-Xn1) / (1.-Xw1);
    double rho_w = problem.materialLaw.wettingPhase.density(temp, p);
    double rho_n = problem.materialLaw.nonwettingPhase.density(temp, p);
    
    Cw1 = rho_w * Xw1;
    Cn1 = rho_n * Xn1;
    Cw2 = rho_w - Cw1;
    Cn2 = rho_n - Cn1;
    C1  = poro* (sat * Cw1 + (1-sat) * Cn1);
    C2  = poro* (sat * Cw2 + (1-sat) * Cn2);
  }
  
  template<class G, class RT>
  void Decoupled2p2c<G,RT>::postupdate(double t, double dt)
  {
  	int size = elementmapper.size();
  	// iterate through leaf grid an evaluate c0 at cell center
    Iterator eendit = indexset.template end<0,All_Partition>();
    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		{
    	int indexi = elementmapper.map(*it);
    	double Z1 = problem.variables.totalConcentration[indexi] / (problem.variables.totalConcentration[indexi] + problem.variables.totalConcentration[elementmapper.size()+indexi]);
    	double C1 = problem.variables.totalConcentration[indexi][0];
    	double C2 = problem.variables.totalConcentration[size+indexi][0];
			flashCalculation(Z1, problem.variables.pressure[indexi], 283.15, problem.porosity(), problem.variables.saturation[indexi][0], C1, C2, problem.variables.wet_c1[indexi][0], problem.variables.nonwet_c1[indexi][0], problem.variables.wet_c2[indexi][0], problem.variables.nonwet_c2[indexi][0]);
			
			double nuw = problem.variables.saturation[indexi] * problem.materialLaw.wettingPhase.density(283.15, problem.variables.pressure[indexi]) / (problem.variables.saturation[indexi] * problem.materialLaw.wettingPhase.density(283.15, problem.variables.pressure[indexi]) + (1-problem.variables.saturation[indexi])*problem.materialLaw.nonwettingPhase.density(283.15, problem.variables.pressure[indexi]));
			double massw = (problem.variables.totalConcentration[indexi][0] + problem.variables.totalConcentration[size+indexi][0]) * nuw;
			double massn = (problem.variables.totalConcentration[indexi][0] + problem.variables.totalConcentration[size+indexi][0]) * (1-nuw);
			double vol = massw / problem.materialLaw.wettingPhase.density(283.15, problem.variables.pressure[indexi]) + massn / problem.materialLaw.nonwettingPhase.density(283.15, problem.variables.pressure[indexi]);
			problem.variables.volErr[indexi] = (vol - 1) / dt;
		}
		timestep = dt;
  }
  
  

}//end namespace Dune

#endif /*DECOUPLED2P2C_HH_*/
