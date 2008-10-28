// $Id$

#ifndef DECOUPLED2P2C_HH
#define DECOUPLED2P2C_HH

// commons:
#include <float.h>
#include <cmath>
#include <string>
#include <fstream>
#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

// transport:
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem2p2c.hh"

// pressure:
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/pardiso/pardiso.hh"

// author: Jochen Fritz
// last change: 27.08.08

namespace Dune
{

/*####################################################*
 * 																										*
 *     CLASS DECLARATION															*
 * 																										*
 *####################################################*/

	//! Implementation of a decoupled formulation of a two phase two component process
	/**
	 * 	The pressure equation is given as \f$ -\frac{\partial V}{\partial p}\frac{\partial p}{\partial t}+\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}\nabla\cdot\left(\sum_{\alpha}C_{\alpha}^{\kappa}\mathbf{v}_{\alpha}\right)=\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}q^{\kappa}\f$
	 *  See paper SPE 99619 for derivation.
	 *  The transport equation is \f$ \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \sum{C_\alpha^\kappa f_\alpha {\bf v}} + q^\kappa \f$
	 */
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

	  typedef typename G::LevelGridView GV;
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename GV::IndexSet IS;
	  typedef typename GV::template Codim<0>::Iterator Iterator;
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
		typedef BlockVector< FieldVector<RT,1> > RepresentationType;
		typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> > VelocityType;

		void initial()
		{
			upd = 0;
			timestep = 0;
			problem.variables.pressure = 1e5;
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

	  void flashCalculation(double Z1, double p, double temp, double poro, double& sat, double& C1, double& C2, double& Xw1, double& Xn1, double& Xw2, double& Xn2);

	  void satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Xw1, double& Xn1, double& Xw2, double& Xn2);

	  void postupdate(double t, double dt);

	  // graphical output
	  void vtkout(const char* name, int k)
	  {
		  problem.variables.volErr *= timestep;
		  problem.variables.vtkout(name, k);
		  problem.variables.volErr /= timestep;
	  }

	private:
		// common variables
		G& grid;
		int level_;
  	const IS& indexset;
  	EM elementmapper;
  	double timestep;
  	const double T; //Temperature

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

	public:
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
	   	  :	grid(g), level_(lev), indexset(g.levelView(lev).indexSet()), reconstruct(rec),
	   	  	numFlux(numFl), diffusivePart(diffPart), alphamax(amax),
	   	  	problem(prob),
	   	  	elementmapper(g, g.levelView(lev).indexSet()),
	   	  	A(g.size(lev, 0),g.size(lev, 0), (2*dim+1)*g.size(lev, 0), BCRSMatrix<MB>::random), f(g.size(lev, 0)),
	   	  	solverName_(solverName), preconditionerName_(preconditionerName),
	   	  	T(283.15)
 	  {
			problem.variables.volErr = 0;
			upd.resize(2*elementmapper.size());
 	  };

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
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
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
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
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
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
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
			FieldMatrix<ct,dim,dim> Ki(this->problem.soil.K(global,*it,local));

			// get porosity
			double poroI = problem.soil.porosity(global, *it, local);

			// get the cell's saturation
			double sati = problem.variables.saturation[indexi];

			// relative permeabilities
			std::vector<double> kr(problem.materialLaw.kr(sati, global, *it, local, T));

			// phase viscosities
			double viscosityL, viscosityG;

			// total mobility and fractional flow factors
			double lambdaI, fw_I, fn_I;

			// derivatives of the fluid volume with respect to mass of compnents and pressure
			double dV_dm1 = 0;
			double dV_dm2 = 0;
			double dV_dp = 0;

			// specific volume of the phases
			double Vg, Vw;

			if (first)
			{
				// total mobility and fractional flow factors
				viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], 0.);
				viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], 0.);
				lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
				fw_I = kr[0] / viscosityL / lambdaI;
				fn_I = kr[1] / viscosityG / lambdaI;

				// specific volume of the phases
				double Vg = 1. / problem.gasPhase.density(T, 1e5, 0.);
				double Vw = 1. / problem.liquidPhase.density(T, 1e5, 0.);

				FieldVector<RT,2> q = problem.q(global,*it,local);
				f[indexi] = volume * (q[0] * Vw + q[1] * Vg);
			}
			else
			{
				// total mobility and fractional flow factors
				viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.wet_X2[indexi]);
				viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.nonwet_X1[indexi]);
				lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
				fw_I = kr[0] / viscosityL / lambdaI;
				fn_I = kr[1] / viscosityG / lambdaI;

				// specific volume of the phases
				Vg = 1. / problem.gasPhase.density(T, problem.variables.pressure[indexi], problem.variables.nonwet_X1[indexi]);
				Vw = 1. / problem.liquidPhase.density(T, problem.variables.pressure[indexi], problem.variables.wet_X2[indexi]);

				// mass of components inside the cell
				double m1 = problem.variables.totalConcentration[indexi] * volume * poroI;
				double m2 = problem.variables.totalConcentration[indexi+elementmapper.size()] * volume * poroI;
				// mass fraction of wetting phase
				double nuw1 = sati / Vw / (sati/Vw + (1-sati)/Vg);
				// actual fluid volume
				double volalt = (m1+m2) * (nuw1 * Vw + (1-nuw1) * Vg);

				// increments for numerical derivatives
				double inc1 = (fabs(upd[indexi][0]) * poroI > 1e-8 /Vw) ?  upd[indexi][0] * poroI : 1e-8/Vw;
				double inc2 =(fabs(upd[indexi+elementmapper.size()][0]) * poroI > 1e-8 / Vg) ?  upd[indexi+elementmapper.size()][0] * poroI : 1e-8 / Vg;
				inc1 *= volume;
				inc2 *= volume;

				// numerical derivative of fluid volume with respect to mass of component 1
				m1 +=  inc1;
				double Z1 = m1 / (m1 + m2);
			  double dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, satt;
				flashCalculation(Z1, problem.variables.pressure[indexi], T, problem.soil.porosity(global, *it, local), satt, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6);
				double nuw = satt / Vw / (satt/Vw + (1-satt)/Vg);
				dV_dm1 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt) /inc1;
				m1 -= inc1;

				// numerical derivative of fluid volume with respect to mass of component 2
				m2 += inc2;
				Z1 = m1 / (m1 + m2);
				flashCalculation(Z1, problem.variables.pressure[indexi], 283.15, problem.porosity(global, *it, local), satt, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6);
				nuw = satt / Vw / (satt/Vw + (1-satt)/Vg);
				dV_dm2 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt)/ inc2;
				m2 -= inc2;

				// numerical derivative of fluid volume with respect to pressure
				double incp = 1e-5;
				double p_ = problem.variables.pressure[indexi] + incp;
				double Vg_ = 1. / problem.gasPhase.density(T, p_, problem.variables.nonwet_X1[indexi]);
				double Vw_ = 1. / problem.liquidPhase.density(T, p_, problem.variables.wet_X2[indexi]);
				dV_dp = ((m1+m2) * (nuw1 * Vw_ + (1-nuw1) * Vg_) - volalt) /incp;

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
					double fw_J, fn_J, lambdaJ;
					double satj = problem.variables.saturation[indexj];
					kr = problem.materialLaw.kr(satj, nbglobal, *outside, nblocal, T);
					if (first)
					{
						viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexj], 0.);
						viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexj], 0.);
						lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
					}
					else
					{
						viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexj], problem.variables.wet_X2[indexj]);
						viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexj], problem.variables.nonwet_X1[indexj]);
						lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
						fw_J = kr[0] / viscosityL / lambdaJ;
						fn_J = kr[1] / viscosityG / lambdaJ;
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
						// phase densities in cell in neighbor
						double rho_w_I = 1 / Vw;
						double rho_n_I = 1 / Vg;
						double rho_w_J = problem.liquidPhase.density(T, problem.variables.pressure[indexj], problem.variables.wet_X2[indexj]);
						double rho_n_J = problem.gasPhase.density(T, problem.variables.pressure[indexj], problem.variables.nonwet_X1[indexj]);
						if (problem.variables.pressure[indexi] > problem.variables.pressure[indexj])
						{
							entry = fabs(
									     dV_dm1 * ( rho_w_I * problem.variables.wet_X1[indexi] * fw_I + rho_n_I * problem.variables.nonwet_X1[indexi] * fn_I)
										 + dV_dm2 * ( rho_w_I * problem.variables.wet_X2[indexi] * fw_I + rho_n_I * problem.variables.nonwet_X2[indexi] * fn_I)
											 );
						}
						else
						{
							entry = fabs(
									     dV_dm1 * ( rho_w_J * problem.variables.wet_X1[indexj] * fw_J + rho_n_J * problem.variables.nonwet_X1[indexj] * fn_J)
										 + dV_dm2 * ( rho_w_J * problem.variables.wet_X2[indexj] * fw_J + rho_n_J * problem.variables.nonwet_X2[indexj] * fn_J)
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

							double satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound, Xw2Bound, Xn2Bound;

							//get boundary condition type for compositional transport
							BoundaryConditions2p2c::Flags bctype = problem.cbctype(faceglobal, *it, facelocalDim);
							if (bctype == BoundaryConditions2p2c::saturation) // saturation given
							{
								satBound = problem.gS(faceglobal, *it, facelocalDim);
								satFlash(satBound, pressBC, T, problem.soil.porosity(global, *it, local), C1Bound, C2Bound, Xw1Bound, Xn1Bound, Xw2Bound, Xn2Bound);
							}
							if (bctype == BoundaryConditions2p2c::concentration) // mass fraction given
							{
								double Z1Bound = problem.gZ(faceglobal, *it, facelocalDim);
								flashCalculation(Z1Bound, pressBC, T, problem.porosity(global, *it, local), satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound, Xw2Bound, Xn2Bound);
							}

							// phase densities in cell and on boundary
							double rho_w_I = 1 / Vw;
							double rho_n_I = 1 / Vg;
							double rho_w_J = problem.liquidPhase.density(T, pressBC, Xw2Bound);
							double rho_n_J = problem.gasPhase.density(T, pressBC, Xn1Bound);

							double entry;
							if (problem.variables.pressure[indexi] > pressBC)
								entry = fabs(
							     dV_dm1 * ( rho_w_I * problem.variables.wet_X1[indexi] * fw_I + rho_n_I * problem.variables.nonwet_X1[indexi] * fn_I)
									 + dV_dm2 * ( rho_w_I * problem.variables.wet_X2[indexi] * fw_I + rho_n_I * problem.variables.nonwet_X2[indexi] * fn_I)
									 );
							else
								entry = fabs(
							     dV_dm1 * ( rho_w_J * Xw1Bound * fw_I + rho_n_J * Xn1Bound * fn_I)
									 + dV_dm2 * ( rho_w_J * Xw2Bound * fw_I + rho_n_J * Xn2Bound * fn_I)
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

			// error reduction routine: volumetric error is damped and inserted to right hand side
			// if damping is not done, the solution method gets unstable!
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
	    	  CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
	    	  solver.apply(problem.variables.pressure, f, r);
	      }
	      else if (solverName_ == "BiCGSTAB") {
	    	  BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
	    	  solver.apply(problem.variables.pressure, f, r);
	      }
	      else
			  DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination " << preconditionerName_
					  << " and " << solverName_ << ".");
	  }
	  else if (preconditionerName_ == "SeqPardiso") {
	      SeqPardiso<MatrixType,Vector,Vector> preconditioner(A);
	      if (solverName_ == "Loop") {
	    	  LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
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

	    Iterator eendit = grid.template lend<0>(level_);
	    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
	    {
	    	// some geometry infos about the cell
	    	GeometryType gt = it->geometry().type();  // cell geometry type
	    	const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
	    	FieldVector<ct,dimworld> global = it->geometry().global(local);  // cell center in global coordinates

		    // cell index
	    	int indexi = elementmapper.map(*it);

	    	// get pressure  in element
	    	double pressi = this->problem.variables.pressure[indexi];

	    	// get absolute permeability
	    	FieldMatrix<ct,dim,dim> Ki(problem.soil.K(global,*it,local));


				// total mobility and fractional flow factors
	    	double sati = problem.variables.saturation[indexi];
				std::vector<double> kr(problem.materialLaw.kr(sati, global, *it, local, T));
				double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.wet_X2[indexi]);
				double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.nonwet_X1[indexi]);
				double lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
	    	double fractionalWI = kr[0] / viscosityL / lambdaI;

	    	double faceVol[2*dim];

	    	// run through all intersections with neighbors and boundary
	    	IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
	    	for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); is!=endit; ++is)
				{
				  // get geometry type of face
				  GeometryType gtf = is.intersectionSelfLocal().type();

				  //Geometry dg = is.intersectionSelfLocal();
				  // local number of facet
				  int numberInSelf = is.numberInSelf();

				  faceVol[numberInSelf] = is.intersectionGlobal().volume();

				  // center in face's reference element
				  const FieldVector<ct,dim-1>& facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

				  // center of face inside volume reference element
				  const FieldVector<ct,dim>& facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);

				  // get normal vector
				  FieldVector<ct,dimworld> unitOuterNormal = is.unitOuterNormal(facelocal);

				  // center of face in global coordinates
				  FieldVector<ct,dimworld> faceglobal = is.intersectionGlobal().global(facelocal);

				  // handle interior face
				  if (is.neighbor())
				    {
				      // access neighbor
				      EntityPointer outside = is.outside();
				      int indexj = elementmapper.map(*outside);

				      // get neighbor pressure and permeability
				      double pressj = this->problem.variables.pressure[indexj];

				      // geometry informations of neighbor
				      GeometryType nbgt = outside->geometry().type(); //geometry type
				      const FieldVector<ct,dim>& nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0); // cell center in local coordinates
				      FieldVector<ct,dimworld> nbglobal = outside->geometry().global(nblocal); // neighbor cell center in global coordinates

				      // distance vector between barycenters
				      FieldVector<ct,dimworld> distVec = global - nbglobal;

				      // compute distance between cell centers
				      double dist = distVec.two_norm();

				      // get absolute permeability
				      FieldMatrix<ct,dim,dim> Kj(problem.soil.K(nbglobal, *outside, nblocal));

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


      				// total mobility and fractional flow factors
      	    	double satj = problem.variables.saturation[indexj];
      				kr = problem.materialLaw.kr(satj, nbglobal, *outside, nblocal, T);
      				viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.wet_X2[indexj]);
      				viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.nonwet_X1[indexj]);
      				double lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
				    	double fractionalWJ = kr[0] / viscosityL / lambdaJ;

				      // compute averaged total mobility
				      // CAREFUL: Harmonic weightig can generate zero matrix entries,
				      // use arithmetic weighting instead:
				      double lambda = 1;
				      double fractionalW;
				      lambda = 0.5 * (lambdaI + lambdaJ);
				      if (hasGravity)
				    	  fractionalW = 0.5 * (fractionalWI + fractionalWJ);

				      FieldVector<ct,dimworld> vTotal(K);
				      vTotal *= lambda * (pressi - pressj) / dist;
				      problem.variables.velocity[indexi][numberInSelf] = vTotal;
				    }
				  // boundary face
				  else
			    {
				      //get boundary condition for boundary face center
					  BoundaryConditions::Flags bctype = problem.pbctype(faceglobal, *it, facelocalDim);
					  if (bctype == BoundaryConditions::dirichlet)
					  {
				      // uniform direction vector between barycenters
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
		Iterator eendit = grid.template lend<0>(level_);
		for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
		{
			int indexi = elementmapper.map(*it);

			// get geometry information of cell
			GeometryType gt = it->geometry().type();
			const FieldVector<ct,dim>&
				local = ReferenceElements<ct,dim>::general(gt).position(0,0);
			FieldVector<ct,dimworld> global = it->geometry().global(local);

			// initial conditions
			double sat_0 = 0;
			double Z1_0 = 0;
			BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *it, local); // get type of initial condition

			if (ictype == BoundaryConditions2p2c::saturation)// saturation initial condition
				sat_0 = problem.S0(global, *it, local);
			else if(ictype == BoundaryConditions2p2c::concentration)			// saturation initial condition
			{
				Z1_0 = problem.Z1_0(global, *it, local);
				double rho_l = problem.liquidPhase.density(T, 1e5, 0.);
				sat_0 = Z1_0 / rho_l;
				sat_0 /= Z1_0 / rho_l + (1 - Z1_0) * problem.materialLaw.nonwettingPhase.density(T, 1e5, 0.);
			}
			else
			{
				DUNE_THROW(Dune::NotImplemented, "Boundary condition " << ictype);
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
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
		{
    	int indexi = elementmapper.map(*it);

			// get geometry information of cell
			GeometryType gt = it->geometry().type();
			const FieldVector<ct,dim>&
				local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
			FieldVector<ct,dimworld> global = it->geometry().global(local);

			// initial conditions
			double sat_0, C1_0, C2_0;
			Dune::BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *it, local);			// get type of initial condition

			if (ictype == Dune::BoundaryConditions2p2c::saturation)  // saturation initial condition
			{
				sat_0 = problem.S0(global, *it, local);
				satFlash(sat_0, problem.variables.pressure[indexi], T, problem.soil.porosity(global, *it, local), C1_0, C2_0, problem.variables.wet_X1[indexi][0], problem.variables.nonwet_X1[indexi][0], problem.variables.wet_X2[indexi][0], problem.variables.nonwet_X2[indexi][0]);
			}
			else if (ictype == Dune::BoundaryConditions2p2c::concentration) // concentration initial condition
			{
				double Z1_0 = problem.Z1_0(global, *it, local);
				flashCalculation(Z1_0, problem.variables.pressure[indexi], T, problem.soil.porosity(global, *it, local), sat_0, C1_0, C2_0, problem.variables.wet_X1[indexi][0], problem.variables.nonwet_X1[indexi][0], problem.variables.wet_X2[indexi][0], problem.variables.nonwet_X2[indexi][0]);
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
		// initialize timestep dt very large
		dt = 1E100;

		// set update vector to zero
		updateVec = 0;

		int which;

		// compute update vector
		Iterator eendit = grid.template lend<0>(level_);
		for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
		{
			// get cell geometry informations
			GeometryType gt = it->geometry().type(); //geometry type
			const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
			FieldVector<ct,dimworld> global = it->geometry().global(local); // cell center in global coordinates
			double volume = it->geometry().integrationElement(local) * ReferenceElements<ct,dim>::general(gt).volume(); // cell volume, assume linear map here

			// cell index
			int indexi = elementmapper.map(*it);

			// get source term
			updateVec[indexi] += problem.q(global, *it, local)[0] * volume;
			updateVec[indexi+elementmapper.size()] += problem.q(global, *it, local)[1] * volume;

	    // get saturation and concentration value at cell center
	    double satI = problem.variables.saturation[indexi];
	    double Xw1_I = problem.variables.wet_X1[indexi];
	    double Xn1_I = problem.variables.nonwet_X1[indexi];
	    double Xw2_I = problem.variables.wet_X2[indexi];
	    double Xn2_I = problem.variables.nonwet_X2[indexi];

			// phase densities in cell
			double rho_w_I = problem.liquidPhase.density(T, problem.variables.pressure[indexi], Xw2_I);
			double rho_n_I = problem.gasPhase.density(T, problem.variables.pressure[indexi], Xn1_I);

			// some variables for time step calculation
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

				// get geometry informations of face
				GeometryType gtf = is.intersectionSelfLocal().type(); //geometry type
				const FieldVector<ct,dim-1>& facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0); // center in face's reference element
				const FieldVector<ct,dim>& facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1); // center of face inside volume reference element

				// get normal vector scaled with volume of face
				FieldVector<ct,dimworld> integrationOuterNormal = is.integrationOuterNormal(facelocal);
				integrationOuterNormal *= ReferenceElements<ct,dim-1>::general(gtf).volume();

				// standardizes velocity
				double velocityIJ = std::max(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / (volume), 0.0);

				// variables for timestep calculation
				double factor, factorC1, factorC2;


				if (is.neighbor()) // handle interior face
				{
					// access neighbor
					EntityPointer outside = is.outside();
					int indexj = elementmapper.map(*outside);

					// neighbor geometry informations
					GeometryType nbgt = outside->geometry().type();
					const FieldVector<ct,dim>& nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
					FieldVector<ct,dimworld> nbglobal = outside->geometry().global(nblocal); // neighbor cell center in global coordinates

					// standardized velocity
					double velocityJI = std::max(-(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / volume), 0.0);

					// distance vector between barycenters
					FieldVector<ct,dimworld> distVec = global - nbglobal;

					// distance between barycenters
					double dist = distVec.two_norm();

					// get saturation and concentration value at neighbor cell center
					double satJ = this->problem.variables.saturation[indexj];
					double Xw1_J = problem.variables.wet_X1[indexj];
					double Xn1_J = problem.variables.nonwet_X1[indexj];
					double Xw2_J = problem.variables.wet_X2[indexj];
					double Xn2_J = problem.variables.nonwet_X2[indexj];

					// phase densities in neighbor
					double rho_w_J = problem.liquidPhase.density(T, problem.variables.pressure[indexj], Xw2_J);
					double rho_n_J = problem.gasPhase.density(T, problem.variables.pressure[indexj], Xn1_J);

					// total mobility and fractional flow factors
					std::vector<double> kr(problem.materialLaw.kr(satI, global, *it, local, T));
					double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], Xw2_I);
					double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], Xn1_I);
					double lambda = kr[0] / viscosityL + kr[1] / viscosityG;
					double fwI = kr[0] / viscosityL / lambda;
					double fnI = kr[1] / viscosityG / lambda;
					// neighbor cell
					viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexj], Xw2_J);
					viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexj], Xn1_J);
					kr = problem.materialLaw.kr(satJ, nbglobal, *outside, nblocal, T);
					lambda = kr[0] / viscosityL + kr[1] / viscosityG;
					double fwJ = kr[0] / viscosityL / lambda;
					double fnJ = kr[1] / viscosityG / lambda;

					// for timestep control
					{
						double coeffW = fwI / satI;
						if(isinf(coeffW) || isnan(coeffW)) coeffW = 0;
						double coeffN = fnI / (1-satI);
						if(isinf(coeffN) || isnan(coeffN)) coeffN = 0;
						factor = (velocityJI - velocityIJ) * std::max(std::max(coeffW, coeffN),1.);
					}

					//  diffFactor =diffPart / volume; TODO include diffusion into timestep control
					factorC1 =
							  velocityJI * Xw1_J * rho_w_J * fwJ /*numFlux(satJ, satI, fwJ, fwI)*/
							- velocityIJ * Xw1_I * rho_w_I * fwI /*numFlux(satI, satJ, fwI, fwJ)*/
							+ velocityJI * Xn1_J * rho_n_J * fnJ /*numFlux(1.0-satJ, 1.0-satI, fnJ, fnI)*/
							- velocityIJ * Xn1_I * rho_n_I * fnI /*numFlux(1.0-satI, 1.0-satJ, fnI, fnJ)*/;
					factorC2 =
							  velocityJI * Xw2_J * rho_w_J * fwJ /*numFlux(satJ, satI, fwJ, fwI)*/
							- velocityIJ * Xw2_I * rho_w_I * fwI /*numFlux(satI, satJ, fwI, fwJ)*/
							+ velocityJI * Xn2_J * rho_n_J * fnJ /*numFlux(1.0-satJ, 1.0-satI, fnJ, fnI)*/
							- velocityIJ * Xn2_I * rho_n_I * fnI /*numFlux(1.0-satI, 1.0-satJ, fnI, fnJ)*/;
				}

				else // handle boundary face
				{
					// face center in globel coordinates
					FieldVector<ct,dim> faceglobal = is.intersectionGlobal().global(facelocal);

					// distance vector between cell and face center
					FieldVector<ct,dimworld> distVec = global - faceglobal;

					// get saturation value at cell center
					double satI = this->problem.variables.saturation[indexi];

					// standardized velocity
					double velocityJI = std::max(-(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / volume), 0.0);

					//get boundary conditions
					BoundaryConditions::Flags pressBCtype = problem.pbctype(faceglobal, *it, facelocalDim);
					if (pressBCtype == BoundaryConditions::dirichlet)
					{
						double satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound, Xw2Bound, Xn2Bound;
						BoundaryConditions2p2c::Flags bctype = problem.cbctype(faceglobal, *it, facelocalDim);
						if (bctype == BoundaryConditions2p2c::saturation)
						{
							satBound = problem.gS(faceglobal, *it, facelocalDim);
							satFlash(satBound, problem.gPress(faceglobal, *it, facelocalDim), T, problem.soil.porosity(global, *it, local), C1Bound, C2Bound, Xw1Bound, Xn1Bound, Xw2Bound, Xn2Bound);
						}
						if (bctype == BoundaryConditions2p2c::concentration)
						{
							double Z1Bound = problem.gZ(faceglobal, *it, facelocalDim);
							flashCalculation(Z1Bound, problem.gPress(faceglobal, *it, facelocalDim), T, problem.porosity(global, *it, local), satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound, Xw2Bound, Xn2Bound);
						}

						double dist = distVec.two_norm();

						// CAREFUL: works only for axisymmetric grids
						if (reconstruct) {
							for (int k = 0; k < dim; k++)
								if (fabs(distVec[k]) > 0.5 * dist)
								{
									satI -= fabs(distVec[k]) / distVec[k] * dist;//*slope[indexi][k];
								}
						}

			    	double Xw1_I = problem.variables.wet_X1[indexi];
			    	double Xn1_I = problem.variables.nonwet_X1[indexi];
			    	double Xw2_I = problem.variables.wet_X2[indexi];
			    	double Xn2_I = problem.variables.nonwet_X2[indexi];

						// phase densities on boundary
						double rho_w_Bound = problem.liquidPhase.density(T, problem.gPress(faceglobal, *it, facelocalDim), Xw2Bound);
						double rho_n_Bound = problem.gasPhase.density(T, problem.gPress(faceglobal, *it, facelocalDim), Xn1Bound);

						// total mobility and fractional flow factors
						std::vector<double> kr(problem.materialLaw.kr(satI, global, *it, local, T));
						double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], Xw2_I);
						double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], Xn1_I);
						double lambda = kr[0] / viscosityL + kr[1] / viscosityG;
			    	double fwI = kr[0] / viscosityL / lambda;
			    	double fnI = kr[1] / viscosityG / lambda;
			    	// neighbor cell
			    	kr = problem.materialLaw.kr(satBound, global, *it, local, T);
						lambda = kr[0] / viscosityL + kr[1] / viscosityG;
			    	double fwBound = kr[0] / viscosityL / lambda;
			    	double fnBound = kr[1] / viscosityG / lambda;

				    // for timestep control
				    {
				    	double coeffW = fwI / satI;
				    	if(isinf(coeffW) || isnan(coeffW)) coeffW = 0;
				    	double coeffN = fnI / (1-satI);
				    	if(isinf(coeffN) || isnan(coeffN)) coeffN = 0;
				    	factor = (velocityJI - velocityIJ) * std::max(std::max(coeffW, coeffN),1.);
				    }

			    	factorC1 =
							+ velocityJI * Xw1Bound * rho_w_Bound * numFlux(satBound, satI, fwBound, fwI)
							- velocityIJ * Xw1_I * rho_w_I * numFlux(satI, satBound, fwI, fwBound)
							+ velocityJI * Xn1Bound * rho_n_Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
							- velocityIJ * Xn1_I * rho_n_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound);
			    	factorC2 =
							 velocityJI * Xw2Bound * rho_w_Bound * numFlux(satBound, satI, fwBound, fwI)
							- velocityIJ * Xw2_I * rho_w_I * numFlux(satI, satBound, fwI, fwBound)
							+ velocityJI * Xn2Bound * rho_n_Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
							- velocityIJ * Xn2_I * rho_n_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound);
					}
					else if (pressBCtype == BoundaryConditions::neumann)
					{
						FieldVector<RT,2> J = problem.J(faceglobal, *it, facelocalDim);
						double faceVol = integrationOuterNormal.two_norm();
						factorC1 = J[0] * faceVol / volume;
						factorC2 = J[1] * faceVol / volume;

						// total mobility and fractional flow factors
						std::vector<double> kr(problem.materialLaw.kr(satI, global, *it, local, T));
						double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.wet_X2[indexi]);
						double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.nonwet_X1[indexi]);
						double lambda = kr[0] / viscosityL + kr[1] / viscosityG;
			    	double fwI = kr[0] / viscosityL / lambda;
			    	double fnI = kr[1] / viscosityG / lambda;

				    // for timestep control
				    {
				    	double coeffW = fwI / satI;
				    	if(isinf(coeffW) || isnan(coeffW)) coeffW = 0;
				    	double coeffN = fnI / (1-satI);
				    	if(isinf(coeffN) || isnan(coeffN)) coeffN = 0;
				    	factor = fabs(J[0] * problem.liquidPhase.density(T, problem.variables.pressure[indexi][0], 0.) + J[1] * problem.gasPhase.density(T, problem.variables.pressure[indexi][0], 0.));
				    	factor *= 0;//std::max(std::max(coeffW, coeffN),1.);
				    }
					}

					else DUNE_THROW(NotImplemented, "there is no process boundary condition implemented");
				}

				// add to update vector
				updateVec[indexi] += factorC1;
				updateVec[elementmapper.size()+indexi] += factorC2;

				// for time step calculation
				if (factor>=0)
					sumfactor += factor;
				else
					sumfactor2 += (-factor);
//				if (diffFactor>=0)
//					sumDiff += diffFactor;
//				else
//					sumDiff += (-diffFactor);
			} // end all intersections
			// compute dt restriction
//			volInc[indexi] = sumfactor - sumfactor2;

			// account for porosity
			double poro = problem.porosity(global, *it, local);

			sumfactor = std::max(sumfactor,sumfactor2) / poro;
			sumDiff = std::max(sumDiff,sumDiff2);
			sumfactor = std::max(sumfactor,100*sumDiff);
			if ( 1./sumfactor < dt)
			{
				dt = 1./sumfactor;
				which= indexi;
			}

		} // end grid traversal
//		printvector(std::cout,updateVec,"update","row");
		return 0;
  } // end function "update"

  template<class G, class RT>
  void Decoupled2p2c<G,RT>::flashCalculation(double Z1, double p, double temp, double poro, double& sat, double& C1, double& C2, double& Xw1, double& Xn1, double& Xw2, double& Xn2)
  {
	  double K1 = problem.liquidPhase.p_vap(temp) / p;
    double K2 = 1. / (p * problem.liquidPhase.henry(temp));
    double xw1 = (1. - K2) / (K1 -K2);
    double xn1 = xw1 * K1;
    Xw1 = xw1 * problem.liquidPhase.molarMass_w()
                / ( xw1 * problem.liquidPhase.molarMass_w() + (1.-xw1) * problem.liquidPhase.molarMass_a() );
    Xn1 = xn1 * problem.liquidPhase.molarMass_w()
                / ( xn1 * problem.liquidPhase.molarMass_w() + (1.-xn1) * problem.liquidPhase.molarMass_a() );
    K1 = Xn1 / Xw1;
    K2 = (1.-Xn1) / (1.-Xw1);
    Xw2 = 1- Xw1;
    Xn2 = 1- Xn1;

    double nu2 = 0;

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

    double rho_w = problem.liquidPhase.density(temp, p, Xw2);
    double rho_n = problem.gasPhase.density(temp, p, Xn1);

    Xw2 = 1 - Xw1;
    Xn2 = 1 - Xn1;

    sat = (1-nu2) / rho_w;
    sat /= ((1-nu2)/rho_w + nu2/rho_n);

    C1 = poro * (Xw1 * sat * rho_w + Xn1 * (1-sat) * rho_n);
    C2 = poro * (Xw2 * sat * rho_w + Xn2 * (1-sat) * rho_n);

  } // end function flashCalculation

  template<class G, class RT>
  void Decoupled2p2c<G,RT>::satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Xw1, double& Xn1, double& Xw2, double& Xn2)
  {
  	if (sat <= 0 || sat >= 1)
	  DUNE_THROW(RangeError,
	  	"Decoupled2p2c :: saturation initial and boundary conditions may not equal zero or one!");
    double K1 = problem.liquidPhase.p_vap(temp) / p;
    double K2 = 1. / (p * problem.liquidPhase.henry(temp));
    double xw1 = (1. - K2) / (K1 -K2);
    double xn1 = xw1 * K1;
    Xw1 = xw1 * problem.liquidPhase.molarMass_w()
                    / ( xw1 * problem.liquidPhase.molarMass_w() + (1.-xw1) * problem.liquidPhase.molarMass_a() );
    Xn1 = xn1 * problem.liquidPhase.molarMass_w()
                    / ( xn1 * problem.liquidPhase.molarMass_w() + (1.-xn1) * problem.liquidPhase.molarMass_a() );
    K1 = Xn1 / Xw1;
    K2 = (1.-Xn1) / (1.-Xw1);

    double rho_w = problem.liquidPhase.density(temp, p, 1.-Xw1);
    double rho_n = problem.gasPhase.density(temp, p, Xn1);

    Xw2 = 1 - Xw1;
    Xn2 = 1 - Xn1;

    C1  = poro* (sat * Xw1 * rho_w + (1-sat) * Xn1 * rho_n);
    C2  = poro* (sat * Xw2 * rho_w + (1-sat) * Xn2 * rho_n);
  }

  template<class G, class RT>
  void Decoupled2p2c<G,RT>::postupdate(double t, double dt)
  {
  	int size = elementmapper.size();
  	// iterate through leaf grid an evaluate c0 at cell center
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
		{
    	int indexi = elementmapper.map(*it);
    	// get cell geometry informations
			GeometryType gt = it->geometry().type(); //geometry type
			const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
			FieldVector<ct,dimworld> global = it->geometry().global(local); // cell center in global coordinates

			double poro = problem.porosity(global, *it, local);
			double rho_l = problem.liquidPhase.density(T, problem.variables.pressure[indexi][0],0.);
			double rho_g = problem.gasPhase.density(T, problem.variables.pressure[indexi][0],0.);

    	double Z1 = problem.variables.totalConcentration[indexi] / (problem.variables.totalConcentration[indexi] + problem.variables.totalConcentration[elementmapper.size()+indexi]);
    	double C1 = problem.variables.totalConcentration[indexi][0];
    	double C2 = problem.variables.totalConcentration[size+indexi][0];
			flashCalculation(Z1, problem.variables.pressure[indexi], T, poro, problem.variables.saturation[indexi][0], C1, C2, problem.variables.wet_X1[indexi][0], problem.variables.nonwet_X1[indexi][0], problem.variables.wet_X2[indexi][0], problem.variables.nonwet_X2[indexi][0]);

			double nuw = problem.variables.saturation[indexi] * rho_l / (problem.variables.saturation[indexi] * rho_l + (1-problem.variables.saturation[indexi]) * rho_g);
			double massw = (problem.variables.totalConcentration[indexi][0] + problem.variables.totalConcentration[size+indexi][0]) * nuw;
			double massn = (problem.variables.totalConcentration[indexi][0] + problem.variables.totalConcentration[size+indexi][0]) * (1-nuw);
			double vol = massw / rho_l + massn / rho_g;
			problem.variables.volErr[indexi] = (vol - poro) / dt;
		}
		timestep = dt;
  }



}//end namespace Dune

#endif /*DECOUPLED2P2C_HH_*/
