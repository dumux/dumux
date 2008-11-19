#ifndef MULTIPHYSICS2P2C_HH
#define MULTIPHYSICS2P2C_HH

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
#include <dune/subgrid/subgrid.hh>

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
// last change: 19.11.08

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
	class Multiphysics2p2c
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

	  // grid typedefs
	  typedef typename G::LevelGridView GV;
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename GV::IndexSet IS;
	  typedef typename GV::template Codim<0>::Iterator Iterator;
	  typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	  typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	  typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
	  typedef typename G::ctype ct;

	  // subgrid typedefs
	  typedef SubGrid<dim,G> SG;
	  typedef typename SG::Traits::LevelIndexSet SIS;
	  typedef typename SG::Traits::template Codim<0>::Entity SubEntity;
	  typedef typename SG::template Codim<0>::LevelIterator SubIterator;
		typedef typename SG::template Codim<0>::LevelIntersectionIterator SubIntersectionIterator;
		typedef typename SG::template Codim<0>::EntityPointer SubEntityPointer;

	  // data typedefs
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
			initializeMatrix();
			initialguess();
			pressure(true, 0);
			transportInitial();
			pressure(false,0);
			transportInitial();
			totalVelocity(0.);
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

		void pressure(bool first, const RT t)
		{
			assemble(first, t);
			solve();
			return;
		}

		void totalVelocity(const RT t);

	  // transport equation functions:
	  void initialguess();

	  void transportInitial();

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

	  int concentrationUpdate(const RT t, RT& dt, RepresentationType& updateVec);

	  void flashCalculation(double Z1, double p, double temp, double poro, double& sat, double& C1, double& C2, double& Xw1, double& Xn1);

	  void satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Xw1, double& Xn1);

	  void postupdate(double t, double dt);

	  // graphical output
	  void vtkout(const char* name, int k)
	  {
		  problem.variables.volErr *= timestep;

			char fname1[128];
			char fname2[128];
			sprintf(fname1, "%s2p2c-%05d", name, k);
			sprintf(fname2, "%s1p2c-%05d", name, k);

		  int globalsize = indexset.size(0);
		  int subsize = subindexset.size(0);
		  RepresentationType C1(subsize);
		  RepresentationType C2(globalsize);
		  for (int i = 0; i < globalsize; i++)
		  {
		  	C2[i] = problem.variables.totalConcentration[i + subsize];
		  	if (i < subsize)
		  		C1[i] = problem.variables.totalConcentration[i];
		  }

	    typedef Dune::VTKWriter<SG, typename SG::LevelGridView> VTKWrite;
	    VTKWrite vtkwriter(subgrid.levelView(level_));
			vtkwriter.addCellData(problem.variables.saturation, "Saturation [-]");
			vtkwriter.addCellData(problem.variables.wet_X1, "Mass fraction 1 in wetting phase [-]");
			vtkwriter.addCellData(problem.variables.nonwet_X1, "Mass fraction 1 in nonwetting phase [-]");
			vtkwriter.addCellData(C1, "Total concentration 1 [kg/m^3]");
			vtkwriter.addCellData(problem.variables.volErr, "volumetric error [-]");
			vtkwriter.write(fname1, Dune::VTKOptions::ascii);

			Dune::VTKWriter<G, typename G::LevelGridView> vtkwriter2(grid.levelView(level_));
			vtkwriter2.addCellData(problem.variables.pressure, "Pressure[Pa]");
			vtkwriter2.addCellData(C2, "Total concentration 2 [kg/m^3]");
			vtkwriter2.write(fname2, Dune::VTKOptions::ascii);

		  problem.variables.volErr /= timestep;
	  }

	  // functions for the interrelations of grid and subgrid
	  // internal method to save the subgrid indices belonging to the hostgrid indices and vice versa
	  void initializeIndexMaps()
	  {

	  	SubIterator endit = subgrid.template lend<0>(level_);
			for (SubIterator it = subgrid.template lbegin<0>(level_); it != endit; ++it)
			{
				// cell index in subgrid
				int subindex = subindexset.index(*it);
				// cell indexi in hostgrid
				int index = indexset.index(*(subgrid.template getHostEntity<0>(*it)));
				// save indices to the mapping vectors
				mH2S[index] = subindex;
				mS2H[subindex] = index;
			}
	  }

	  // this function prevents accidental changes of the mapping vectors
	  // returns the index on the subgrid for a given index on the hostgrid
	  inline int mapHost2Sub(int index) const
	  {
	  	return mH2S[index];
	  }

	  // this function prevents accidental changes of the mapping vectors
	  // returns the index on the hostgrid for a given index on the subgrid
	  inline int mapSub2Host(int subindex) const
	  {
	  	if (subindex > mS2H.size() - 1)
	  		DUNE_THROW(RangeError, "given subindex exceeds size of subgrid!");
	  	return mS2H[subindex];
	  }

	private:
		// common variables
		G& grid;
		SG& subgrid;
		int level_;
  	const IS& indexset;
  	const SIS& subindexset;
  	EM elementmapper;
  	double timestep;
  	const double T; //Temperature
  	// index maps
  	std::vector<int> mH2S;
  	std::vector<int> mS2H;

		// variables for transport equation:
		TransportProblem2p2c<G, RT>& problem;
		RepresentationType upd;

		typename VariableClass2p2c<G,RT>::VelType factorcheck;

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
		Multiphysics2p2c(
				G& g,
				SG& sg,
				TransportProblem2p2c<G, RT>& prob,
				int lev = 0,
	   	  DiffusivePart<G, RT>& diffPart = *(new DiffusivePart<G, RT>),
	   	  bool rec = false,
	   	  double amax = 0.8,
	   	  const NumericalFlux<RT>& numFl = *(new Upwind<RT>),
	   	  const std::string solverName = "BiCGSTAB",
			  const std::string preconditionerName = "SeqILU0" )
	   	  :	grid(g), subgrid(sg), level_(lev), indexset(g.levelView(lev).indexSet()), subindexset(sg.levelIndexSet(lev)), reconstruct(rec),
	   	  	numFlux(numFl), diffusivePart(diffPart), alphamax(amax),
	   	  	problem(prob),
	   	  	elementmapper(g, g.levelView(lev).indexSet()),
	   	  	A(g.size(lev, 0),g.size(lev, 0), (2*dim+1)*g.size(lev, 0), BCRSMatrix<MB>::random), f(g.size(lev, 0)),
	   	  	solverName_(solverName), preconditionerName_(preconditionerName),
	   	  	T(283.15), mH2S(grid.size(lev,0), -1), mS2H(subgrid.size(lev,0))
 	  {
			problem.variables.volErr = 0;
			int subsize = subgrid.size(lev,0);
			int hostsize = grid.size(lev,0);
			upd.resize(subsize + hostsize);
			problem.variables.saturation.resize(subsize);
			problem.variables.totalConcentration.resize(subsize + hostsize);
			problem.variables.wet_X1.resize(subsize);
			problem.variables.nonwet_X1.resize(subsize);
			problem.variables.wet_X2.resize(subsize);
			problem.variables.nonwet_X2.resize(subsize);
			problem.variables.volErr.resize(subsize);
			factorcheck.resize(hostsize);
			initializeIndexMaps();
 	  };

	}; //end class declaration



	/*####################################################*
	 * 																										*
	 *     FUNCTION DEFINITIONS 1: PRESSURE EQUATION			*
	 * 																										*
	 *####################################################*/

  template<class G, class RT>
  void Multiphysics2p2c<G, RT>::initializeMatrix()
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
  void Multiphysics2p2c<G, RT>::assemble(bool first, const RT t=0)
	{
		// initialization: set matrix A to zero
			A = 0;

		if (first) problem.variables.pressure = 1e5;

    // iterate over all cells in the grid
    Iterator eendit = grid.template lend<0>(level_);
    Iterator it = grid.template lbegin<0>(level_);
    for (; it != eendit; ++it)
    {
	    // get geometry infos about the cell...
			GeometryType gt = it->geometry().type(); // cell geometry type
			const FieldVector<ct,dim>&
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
			FieldVector<ct,dim> global = it->geometry().global(local); 			// global coordinate of cell center
			double volume = it->geometry().integrationElement(local)
				*ReferenceElements<ct,dim>::general(gt).volume(); // cell volume

			// cell index in host grid and in suubgrid
			int indexi = indexset.index(*it);
			int subindexi = mapHost2Sub(indexi);

			// get absolute permeability
			FieldMatrix<ct,dim,dim> Ki(this->problem.soil.K(global,*it,local));

			// get porosity
			double poroI = problem.soil.porosity(global, *it, local);

			// phase viscosities
			double viscosityL, viscosityG;

			// total mobility and fractional flow factors
			double satI, lambdaI, fw_I, fn_I;

			// get the cell's saturation
			// if subindexi == 1, the subgrid does not contain this entity! (see function initializeIndexMaps)
			if (subindexi != -1) // entity is contained in subgrid
				satI = problem.variables.saturation[subindexi];
			else //entity lies outside of subgrid
				satI = 1.;

			// relative permeabilities
			std::vector<double> kr(problem.materialLaw.kr(satI, global, *it, local, T));

			// derivatives of the fluid volume with respect to mass of compnents and pressure
			double dV_dm1 = 0;
			double dV_dm2 = 0;
			double dV_dp = 0;

			// specific volume of the phases
			double Vg, Vw;

			if (first || subindexi == -1)
			{
				// total mobility and fractional flow factors
				viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], 0.);
				viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], 0.);
				lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
				fw_I = kr[0] / viscosityL / lambdaI;
				fn_I = kr[1] / viscosityG / lambdaI;

				// specific volume of the phases
				Vg = 1. / problem.gasPhase.density(T, problem.variables.pressure[indexi], 0.);
				Vw = 1. / problem.liquidPhase.density(T, problem.variables.pressure[indexi], 0.);

				FieldVector<RT,2> q = problem.q(global,*it,local);
				f[indexi] = volume * (q[0] * Vw + q[1] * Vg);
			}
			else
			{
				// total mobility and fractional flow factors
				viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], (1. - problem.variables.wet_X1[subindexi]));
				viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.nonwet_X1[subindexi]);
				lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
				fw_I = kr[0] / viscosityL / lambdaI;
				fn_I = kr[1] / viscosityG / lambdaI;

				// specific volume of the phases
				Vg = 1. / problem.gasPhase.density(T, problem.variables.pressure[indexi], problem.variables.nonwet_X1[subindexi]);
				Vw = 1. / problem.liquidPhase.density(T, problem.variables.pressure[indexi], (1. - problem.variables.wet_X1[subindexi]));

				// mass of components inside the cell
				double m1 = problem.variables.totalConcentration[subindexi] * volume * poroI;
				double m2 = problem.variables.totalConcentration[indexi+subgrid.size(0)] * volume * poroI;
				// mass fraction of wetting phase
				double nuw1 = satI / Vw / (satI/Vw + (1-satI)/Vg);
				// actual fluid volume
				double volalt = (m1+m2) * (nuw1 * Vw + (1-nuw1) * Vg);

				// increments for numerical derivatives
				double inc1 = (fabs(upd[subindexi][0]) * poroI > 1e-8 /Vw) ?  upd[subindexi][0] * poroI : 1e-8/Vw;
				double inc2 =(fabs(upd[indexi+subgrid.size(0)][0]) * poroI > 1e-8 / Vg) ?  upd[indexi+subgrid.size(0)][0] * poroI : 1e-8 / Vg;
				inc1 *= volume;
				inc2 *= volume;

				// numerical derivative of fluid volume with respect to mass of component 1
				m1 +=  inc1;
				double Z1 = m1 / (m1 + m2);
			  double dummy1, dummy2, dummy3, dummy4, satt;
				flashCalculation(Z1, problem.variables.pressure[indexi], T, poroI, satt, dummy1, dummy2, dummy3, dummy4);
				double nuw = satt / Vw / (satt/Vw + (1-satt)/Vg);
				dV_dm1 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt) /inc1;
				m1 -= inc1;

				// numerical derivative of fluid volume with respect to mass of component 2
				m2 += inc2;
				Z1 = m1 / (m1 + m2);
				flashCalculation(Z1, problem.variables.pressure[indexi], 283.15, poroI, satt, dummy1, dummy2, dummy3, dummy4);
				nuw = satt / Vw / (satt/Vw + (1-satt)/Vg);
				dV_dm2 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt)/ inc2;
				m2 -= inc2;

				// numerical derivative of fluid volume with respect to pressure
				double incp = 1e-5;
				double p_ = problem.variables.pressure[indexi] + incp;
				double Vg_ = 1. / problem.gasPhase.density(T, p_, problem.variables.nonwet_X1[subindexi]);
				double Vw_ = 1. / problem.liquidPhase.density(T, p_, (1. - problem.variables.wet_X1[subindexi]));
				dV_dp = ((m1+m2) * (nuw1 * Vw_ + (1-nuw1) * Vg_) - volalt) /incp;

				// right hand side entry: sources
				FieldVector<RT,2> q = problem.q(global,*it,local);
				f[indexi] = volume * (dV_dm1 * q[0] + dV_dm2 * q[1]);
			}

			// iterate over all faces of the cell
			IntersectionIterator endis = it->ilevelend();
			for (IntersectionIterator is = it->ilevelbegin(); is!=endis; ++is)
			{
	    	// some geometry informations of the face
		    GeometryType gtf = is.intersectionSelfLocal().type(); // get geometry type of face
		    const FieldVector<ct,dim-1>& facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0); // center in face's reference element
		    const FieldVector<ct,dim>& facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1); // center of face inside volume reference element
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

					int indexj = indexset.index(*outside);
					int subindexj = mapHost2Sub(indexj);

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
					double satJ, fw_J, fn_J, lambdaJ;

					if (subindexj != -1)
						satJ = problem.variables.saturation[subindexj];
					else
						satJ = 1.;

					kr = problem.materialLaw.kr(satJ, nbglobal, *outside, nblocal, T);
					if (first || subindexj == -1)
					{
						viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexj], 0.);
						viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexj], 0.);
						lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
					}
					else
					{
						viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexj], (1. - problem.variables.wet_X1[subindexj]));
						viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexj], problem.variables.nonwet_X1[subindexj]);
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
					if (first || subindexi == -1) // first guess pressure calculation or the current cell is not contained in subgrid
						entry = fabs(lambda*faceVol*(K*distVec)/(dist*dist));
					else if (subindexj != -1) // current cell and neighbor are contained in subgrid
					{
						// phase densities in cell in neighbor
						double rho_w_I = 1 / Vw;
						double rho_n_I = 1 / Vg;
						double rho_w_J = problem.liquidPhase.density(T, problem.variables.pressure[indexj], (1. - problem.variables.wet_X1[subindexj]));
						double rho_n_J = problem.gasPhase.density(T, problem.variables.pressure[indexj], problem.variables.nonwet_X1[subindexj]);
						if (problem.variables.pressure[indexi] > problem.variables.pressure[indexj])
						{
							entry = fabs(
									     dV_dm1 * ( rho_w_I * problem.variables.wet_X1[subindexi] * fw_I + rho_n_I * problem.variables.nonwet_X1[subindexi] * fn_I)
										 + dV_dm2 * ( rho_w_I * (1. - problem.variables.wet_X1[subindexi]) * fw_I + rho_n_I * (1. - problem.variables.nonwet_X1[subindexi]) * fn_I)
											 );
						}
						else
						{
							entry = fabs(
									     dV_dm1 * ( rho_w_J * problem.variables.wet_X1[subindexj] * fw_J + rho_n_J * problem.variables.nonwet_X1[subindexj] * fn_J)
										 + dV_dm2 * ( rho_w_J * (1. - problem.variables.wet_X1[subindexj]) * fw_J + rho_n_J * (1. - problem.variables.nonwet_X1[subindexj]) * fn_J)
											 );
						}
						entry *= lambda * fabs(faceVol*(K*distVec)/(dist*dist));
					}
					else // current cell is contained in subgrid, neighbor is not
					{
						// phase densities in cell in neighbor
						double rho_w_I = 1 / Vw;
						double rho_n_I = 1 / Vg;
						double rho_w_J = problem.liquidPhase.density(T, problem.variables.pressure[indexj], 0.);
						double rho_n_J = problem.gasPhase.density(T, problem.variables.pressure[indexj], 0.);
						if (problem.variables.pressure[indexi] > problem.variables.pressure[indexj])
						{
							entry = fabs(
											 dV_dm1 * ( rho_w_I * problem.variables.wet_X1[subindexi] * fw_I + rho_n_I * problem.variables.nonwet_X1[subindexi] * fn_I)
										 + dV_dm2 * ( rho_w_I * (1. - problem.variables.wet_X1[subindexi]) * fw_I + rho_n_I * (1. - problem.variables.nonwet_X1[subindexi]) * fn_I)
											 );
						}
						else
						{
							entry = fabs(
											 dV_dm1 * ( rho_w_J * (1. - problem.variables.totalConcentration[indexj + subindexset.size(0)]) * fw_J + rho_n_J * 0. * fn_J)
										 + dV_dm2 * ( rho_w_J * problem.variables.totalConcentration[indexj + subindexset.size(0)] * fw_J + rho_n_J * 0. * fn_J)
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
					FieldVector<ct,dimworld> faceglobal = is.intersectionGlobal().global(facelocal);

					// compute total mobility
					double lambda = lambdaI;

					//get boundary condition for boundary face center
					BoundaryConditions::Flags bctype = problem.pbctype(faceglobal, *it, facelocalDim);


					if (bctype == BoundaryConditions::dirichlet) //dirichlet boundary
					{
						// distance vector to boundary face
						FieldVector<ct,dimworld> distVec(global - faceglobal);
						double dist = distVec.two_norm();
						if (first || subindexi == -1)
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
								satFlash(satBound, pressBC, T, problem.soil.porosity(global, *it, local), C1Bound, C2Bound, Xw1Bound, Xn1Bound);
							}
							if (bctype == BoundaryConditions2p2c::concentration) // mass fraction given
							{
								double Z1Bound = problem.gZ(faceglobal, *it, facelocalDim);
								flashCalculation(Z1Bound, pressBC, T, problem.porosity(global, *it, local), satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound);
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
									 + dV_dm2 * ( rho_w_I * (1. - problem.variables.wet_X1[indexi]) * fw_I + rho_n_I * (1. - problem.variables.nonwet_X1[indexi]) * fn_I)
									 );
							else
								entry = fabs(
									 dV_dm1 * ( rho_w_J * Xw1Bound * fw_I + rho_n_J * Xn1Bound * fn_I)
									 + dV_dm2 * ( rho_w_J * (1. - Xw2Bound) * fw_I + rho_n_J * (1.- Xn1Bound) * fn_I)
									 );

							entry *= - lambda * faceVol*(Kni*distVec)/(dist*dist);

							// set diagonal entry and right hand side entry
							A[indexi][indexi] += entry;
							f[indexi] += entry * pressBC;
						}
					}
					else //neumann boundary
					{
						FieldVector<RT,2> J = problem.J(faceglobal, *it, facelocalDim);
						if (first || subindexi == -1)
							f[indexi] += faceVol * (J[0] * Vw + J[1] * Vg);
						else
						{
							f[indexi] += faceVol* (J[0] * dV_dm1 + J[1] * dV_dm2);
						}
					}
			  }
			} // end all intersections

			// compressibility term
			if (!first && subindexi != -1)
			{
				if (timestep != 0.)
				{
					A[indexi][indexi] -= dV_dp / timestep;
					f[indexi] -= problem.variables.pressure[indexi] *dV_dp / timestep;
				}

				// error reduction routine: volumetric error is damped and inserted to right hand side
				// if damping is not done, the solution method gets unstable!
				double maxErr = fabs(problem.variables.volErr.infinity_norm());
				double erri = fabs(problem.variables.volErr[subindexi]);
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
						f[indexi] += fac* (1-x_mi*(lofac/fac-1)/(x_lo-x_mi) + (lofac/fac-1)/(x_lo-x_mi)*erri/maxErr)* problem.variables.volErr[subindexi]*volume;
					else
						f[indexi] += fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxErr) * problem.variables.volErr[subindexi]*volume;
				}
			}
	  } // end grid traversal

//	    printmatrix(std::cout,A,"stiffnesmatrix","row");
//	    printvector(std::cout,f,"right hand side","row");

	    return;
	} // end function assemble


  template<class G, class RT>
  void Multiphysics2p2c<G, RT>::solve()
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
	void Multiphysics2p2c<G, RT>::totalVelocity(const RT t=0)
  {
  	//TODO remove DEBUG --->
  	int highc, highf;
  	double highval = 0.;
  	// <--- DEBUG

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
			int indexi = indexset.index(*it);
			int subindexi = mapHost2Sub(indexi);

			// get pressure  in element
			double pressi = this->problem.variables.pressure[indexi];

			// get absolute permeability
			FieldMatrix<ct,dim,dim> Ki(problem.soil.K(global,*it,local));

			double sati, lambdaI, fractionalWI;

			if (subindexi != -1) // cell is contained in subgrid
			{
				// total mobility and fractional flow factors
				sati = problem.variables.saturation[subindexi];
				std::vector<double> kr(problem.materialLaw.kr(sati, global, *it, local, T));
				double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], 1. - problem.variables.wet_X1[subindexi]);
				double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], problem.variables.nonwet_X1[subindexi]);
				lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
				fractionalWI = kr[0] / viscosityL / lambdaI;
			}
			else
			{
				// total mobility and fractional flow factors
				sati = 1.;
				std::vector<double> kr(problem.materialLaw.kr(sati, global, *it, local, T));
				double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], 0.);
				double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], 0.);
				lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
				fractionalWI = kr[0] / viscosityL / lambdaI;
			}

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
					int subindexj = mapHost2Sub(indexj);

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

					double satj, lambdaJ, fractionalWJ;

					if (subindexj != -1)
					{
						// total mobility and fractional flow factors
						double satj = problem.variables.saturation[subindexj];
						std::vector<double> kr = problem.materialLaw.kr(satj, nbglobal, *outside, nblocal, T);
						double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexj], 1. - problem.variables.wet_X1[subindexj]);
						double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexj], problem.variables.nonwet_X1[subindexj]);
						lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
						fractionalWJ = kr[0] / viscosityL / lambdaJ;
					}
					else
					{
						// total mobility and fractional flow factors
						double satj = 1.;
						std::vector<double> kr = problem.materialLaw.kr(satj, nbglobal, *outside, nblocal, T);
						double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexj], 0.);
						double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexj], 0.);
						lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
						fractionalWJ = kr[0] / viscosityL / lambdaJ;
					}

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

					//TODO rmove DEBUG --->
					if (vTotal.two_norm() > highval)
					{
						highc = indexi;
						highf = numberInSelf;
						highval = vTotal.two_norm();
					}
					// <--- DEBUG
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

						//TODO remove DEBUG --->
						if (vTotal.two_norm() > highval)
						{
							highc = indexi;
							highf = numberInSelf;
							highval = vTotal.two_norm();
						}
						// <--- DEBUG

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

		//TODO remove DEBUG --->
		std::cout << highval << " at cell "<<highc<<" and face " <<highf <<std::endl;
		// <---DEBUG

	return;
  }

	/*####################################################*
	 * 																										*
	 *     FUNCTION DEFINITIONS 2: TRANSPORT EQUATION			*
	 * 																										*
	 *####################################################*/

  template<class G, class RT>
  void Multiphysics2p2c<G,RT>::initialguess()
  {
		// iterate through leaf grid an evaluate c0 at cell center
		SubIterator endit = subgrid.template lend<0>(level_);
		for (SubIterator it = subgrid.template lbegin<0>(level_); it != endit; ++it)
		{
			EntityPointer epit = subgrid.template getHostEntity<0>(*it);
			int subindex = subindexset.index(*it);

			// get geometry information of cell
			GeometryType gt = it->geometry().type();
			const FieldVector<ct,dim>&
				local = ReferenceElements<ct,dim>::general(gt).position(0,0);
			FieldVector<ct,dimworld> global = it->geometry().global(local);

			// initial conditions
			double sat_0;
			double Z1_0;
			BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *epit, local); // get type of initial condition

			if (ictype == BoundaryConditions2p2c::saturation)// saturation initial condition
				sat_0 = problem.S0(global, *epit, local);
			else if(ictype == BoundaryConditions2p2c::concentration)			// saturation initial condition
			{
				Z1_0 = problem.Z1_0(global, *epit, local);
				double rho_l = problem.liquidPhase.density(T, 1e5, 0.);
				sat_0 = Z1_0 / rho_l;
				sat_0 /= Z1_0 / rho_l + (1 - Z1_0) * problem.materialLaw.nonwettingPhase.density(T, 1e5, 0.);
			}
			else
			{
				DUNE_THROW(Dune::NotImplemented, "Boundary condition " << ictype);
			}

			// initialize cell saturation
			this->problem.variables.saturation[subindex][0] = sat_0;
		}
	return;
  }//end function initialguess



  template<class G, class RT>
  void Multiphysics2p2c<G,RT>::transportInitial()
	{
    // iterate through grid an evaluate c0 at cell center
		Iterator endit = grid.template lend<0>(level_);
		for (Iterator it = grid.template lbegin<0>(level_); it != endit; ++it)
		{
			int index = indexset.index(*it);
			int subindex = mapHost2Sub(index);

			// get geometry information of cell
			GeometryType gt = it->geometry().type();
			const FieldVector<ct,dim>&
				local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
			FieldVector<ct,dimworld> global = it->geometry().global(local);

			// initial conditions
			double sat_0, C1_0, C2_0;
			Dune::BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *it, local);			// get type of initial condition

			if (subindex == -1) // this cell is not contained in subgrid
			{
				// initial conditions
				double sat_0, C1_0, C2_0;
				Dune::BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *it, local);			// get type of initial condition

				if (ictype == Dune::BoundaryConditions2p2c::saturation)  // saturation initial condition
				{
					std::cout<<"do not give saturtion initial conditions for 1p domain! index no: "<< index <<" location: "<< global[0] << ", " << global[1] << std::endl;
					sat_0 = problem.S0(global, *it, local);
					double dummy1, dummy2;
					satFlash(sat_0, problem.variables.pressure[index], T, problem.soil.porosity(global, *it, local), C1_0, C2_0, dummy1, dummy2);
				}
				else if (ictype == Dune::BoundaryConditions2p2c::concentration) // concentration initial condition
				{
					double Z1_0 = problem.Z1_0(global, *it, local);
					double dummy1,dummy2, dummy3;
					flashCalculation(Z1_0, problem.variables.pressure[index], T, problem.soil.porosity(global, *it, local), dummy1, C1_0, C2_0, dummy2, dummy3);
					if (dummy1 != 1.)
						std::cout<<"Z1_0 to low, gas phase appears! index no: "<< index <<" location: "<< global[0] << ", " << global[1] << std::endl;
				}
				problem.variables.totalConcentration[index + subindexset.size(0)] = C2_0;
			}
			else // this cell is contained in subgrid
			{
				if (ictype == Dune::BoundaryConditions2p2c::saturation)  // saturation initial condition
				{
					sat_0 = problem.S0(global, *it, local);
					satFlash(sat_0, problem.variables.pressure[index], T, problem.soil.porosity(global, *it, local), C1_0, C2_0, problem.variables.wet_X1[subindex][0], problem.variables.nonwet_X1[subindex][0]);
				}
				else if (ictype == Dune::BoundaryConditions2p2c::concentration) // concentration initial condition
				{
					double Z1_0 = problem.Z1_0(global, *it, local);
					flashCalculation(Z1_0, problem.variables.pressure[index], T, problem.soil.porosity(global, *it, local), sat_0, C1_0, C2_0, problem.variables.wet_X1[subindex][0], problem.variables.nonwet_X1[subindex][0]);
				}
				// initialize cell concentration
				problem.variables.totalConcentration[subindex] = C1_0;
				problem.variables.totalConcentration[index + subindexset.size(0)] = C2_0;
				this->problem.variables.saturation[subindex][0] = sat_0;
			}
		}
		return;
	} //end function transportInitial

  template<class G, class RT>
  int Multiphysics2p2c<G,RT>::concentrationUpdate(const RT t, RT& dt, RepresentationType& updateVec)
  {
		// initialize timestep dt very large
		dt = 1E100;

		// set update vector to zero
		updateVec = 0;

		int which;

		// iterate over grid
		Iterator endit = grid.template lend<0>(level_);
		for (Iterator it = grid.template lbegin<0>(level_); it != endit; ++it)
		{
			// get cell geometry informations
			GeometryType gt = it->geometry().type(); //geometry type
			const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
			FieldVector<ct,dimworld> global = it->geometry().global(local); // cell center in global coordinates
			double volume = it->geometry().integrationElement(local) * ReferenceElements<ct,dim>::general(gt).volume(); // cell volume, assume linear map here

			// cell index
			int indexi = indexset.index(*it);
			int subindexi = mapHost2Sub(indexi);

			// get source term
			updateVec[indexi + subindexset.size(0)] += problem.q(global, *it, local)[1] * volume;

			// porosity
			double poroI = problem.soil.porosity(global, *it, local);

			double rho_w_I, rho_n_I, fwI, fnI;
			double satI, Xw1_I, Xn1_I;
			if (subindexi == -1)
			{
				satI = 1.;
				fwI = 1.;
				fnI = 0.;
				rho_w_I = problem.liquidPhase.density(T, problem.variables.pressure[indexi], 0.);
				rho_n_I = problem.gasPhase.density(T, problem.variables.pressure[indexi], 0.);
			}
			else
			{
				satI = problem.variables.saturation[subindexi];
				std::vector<double> kr = problem.materialLaw.kr(satI, global, *it, local, T);
				Xw1_I = problem.variables.wet_X1[subindexi];
				Xn1_I = problem.variables.nonwet_X1[subindexi];
				double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexi], 1. - Xw1_I);
				double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexi], Xn1_I);
				double lambda = kr[0] / viscosityL + kr[1] / viscosityG;
				fwI = kr[0] / viscosityL / lambda;
				fnI = kr[1] / viscosityG / lambda;
				rho_w_I = problem.liquidPhase.density(T, problem.variables.pressure[indexi], 1. - Xw1_I);
				rho_n_I = problem.gasPhase.density(T, problem.variables.pressure[indexi], Xn1_I);
			}

			// some variables for time step calculation
			double sumfactor = 0;
			double sumfactor2 = 0;
			double sumDiff = 0;
			double sumDiff2 = 0;

			// run through all intersections with neighbors and boundary
			IntersectionIterator endis = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it);	is!=endis; ++is)
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

				// standardized velocity
				double velocityIJ = std::max(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / (volume), 0.0);

				// variables for timestep calculation
				double factor, factorC1, factorC2;


				if (is.neighbor()) // handle interior face
				{
					// access neighbor
					EntityPointer outside = is.outside();

					int indexj = elementmapper.map(*outside);
					int subindexj = mapHost2Sub(indexj);

					// neighbor geometry informations
					GeometryType nbgt = outside->geometry().type();
					const FieldVector<ct,dim>& nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
					FieldVector<ct,dimworld> nbglobal = outside->geometry().global(nblocal); // neighbor cell center in global coordinates

					// standardized velocity
					double velocityJI = std::max(-(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / volume), 0.0);

					// porosity
					double poroJ = problem.soil.porosity(nbglobal, *outside, nblocal);

					double rho_w_J, rho_n_J, fwJ, fnJ;
					double satJ, Xw1_J, Xn1_J;
					if (subindexj == -1) //neighbor cell lies outside of subdomain.
					{
						satJ = 1.;
						fwJ = 1.;
						fnJ = 0.;
						rho_w_J = problem.liquidPhase.density(T, problem.variables.pressure[indexj], 0.);
						rho_n_J = problem.gasPhase.density(T, problem.variables.pressure[indexj], 0.);

						if (subindexi == -1) //current cell lise outside of subdomain
							factorC2 =
									velocityJI * problem.variables.totalConcentration[indexj + subindexset.size(0)] / poroJ * numFlux(satJ, satI, fwJ, fwI)
								- velocityIJ * problem.variables.totalConcentration[indexi + subindexset.size(0)] / poroI * numFlux(satI, satJ, fwI, fwJ);
						else // current cell lies inside subdomian
						{
							factorC2 =
									velocityJI * problem.variables.totalConcentration[indexj + subindexset.size(0)] / poroJ * numFlux(satJ, satI, fwJ, fwI)
								- velocityIJ * rho_w_I * (1. - Xw1_I) * numFlux(satI, satJ, fwI, fwJ)
								- velocityIJ * rho_n_I * (1. - Xn1_I) * numFlux(1.-satI, 1.-satJ, fnI, fnJ);
									// actually, the last term must be zero!!! Otherwise the assumption of having only 1p transport in global domain is violated!
							factorC1 =
									velocityJI * (rho_w_J - problem.variables.totalConcentration[indexj + subindexset.size(0)] / poroJ) * numFlux(satJ, satI, fwJ, fwI)
								-	velocityIJ * rho_w_I * Xw1_I * numFlux(satI, satJ, fwI, fwJ)
								-	velocityIJ * rho_w_J * Xw1_J * numFlux(1.-satI, 1.-satJ, fnI, fnJ);
									// actually, the last term must be zero!!! Otherwise the assumption of having only 1p transport in global domain is violated!
						}
					}
					else // neighbor cell lies inside subdomain
					{
						satJ = problem.variables.saturation[subindexj];
						std::vector<double> kr = problem.materialLaw.kr(satJ, nbglobal, *outside, nblocal, T);
						Xw1_J = problem.variables.wet_X1[subindexj];
						Xn1_J = problem.variables.nonwet_X1[subindexj];
						double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[indexj], 1. - Xw1_J);
						double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[indexj], Xn1_J);
						double lambda = kr[0] / viscosityL + kr[1] / viscosityG;
						fwJ = kr[0] / viscosityL / lambda;
						fnJ = kr[1] / viscosityG / lambda;
						rho_w_J = problem.liquidPhase.density(T, problem.variables.pressure[indexj], 1. - Xw1_J);
						rho_n_J = problem.gasPhase.density(T, problem.variables.pressure[indexj], Xn1_J);

						if (subindexi == -1) // current cell lies outside subdomain
						{
							factorC2 =
								velocityJI * rho_w_J * (1. - Xw1_J) * numFlux(satJ, satI, fwJ, fwI)
							- velocityIJ * (problem.variables.totalConcentration[indexi + subindexset.size(0)] / poroI) * numFlux(satI, satJ, fwI, fwJ);
							+ velocityJI * rho_n_J * (1. - Xn1_J) * numFlux(1.-satJ, 1.-satI, fnJ, fnI);
								// actually, the last term must be zero!!! Otherwise the assumption of having only 1p transport in global domain is violated!
						}
						else
						{
							factorC2 =
										velocityJI * rho_w_J * (1. - Xw1_J) * numFlux(satJ, satI, fwJ, fwI)
									- velocityIJ * rho_w_I * (1. - Xw1_I) * numFlux(satI, satJ, fwI, fwJ)
									+ velocityJI * rho_n_J * (1. - Xn1_J) * numFlux(1.-satJ, 1.-satI, fnJ, fnI)
									- velocityIJ * rho_n_I * (1. - Xn1_I) * numFlux(1.-satI, 1.-satJ, fnI, fnJ);

							factorC1 =
										velocityJI * rho_w_J * Xw1_J * numFlux(satJ, satI, fwJ, fwI)
									- velocityIJ * rho_w_I * Xw1_I * numFlux(satI, satJ, fwI, fwJ)
									+ velocityJI * rho_n_J * Xn1_J * numFlux(1.-satJ, 1.-satI, fnJ, fnI)
									- velocityIJ * rho_n_I * Xn1_I * numFlux(1.-satI, 1.-satJ, fnI, fnJ);
						}
					}

					// for timestep control
					{
						double coeffW = fwI / satI;
						if(isinf(coeffW) || isnan(coeffW)) coeffW = 0;
						double coeffN = fnI / (1-satI);
						if(isinf(coeffN) || isnan(coeffN)) coeffN = 0;
						factor = (velocityJI - velocityIJ) * std::max(std::max(coeffW, coeffN),1.);
					}
				}

				else // handle boundary face
				{
					// cell center in global coordinates
					FieldVector<ct,dimworld> global = it->geometry().global(local);

					// face center in globel coordinates
					FieldVector<ct,dim> faceglobal = is.intersectionGlobal().global(facelocal);

					// distance vector between cell and face center
					FieldVector<ct,dimworld> distVec = global - faceglobal;
					double dist = distVec.two_norm();

					// standardized velocity
					double velocityJI = std::max(-(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / volume), 0.0);

					//get boundary conditions
					BoundaryConditions::Flags pressBCtype = problem.pbctype(faceglobal, *it, facelocalDim);
					if (pressBCtype == BoundaryConditions::dirichlet)
					{
						// get pressure on boundary
						double pressBound = problem.gPress(faceglobal, *it, facelocalDim);

						double satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound;
						BoundaryConditions2p2c::Flags bctype = problem.cbctype(faceglobal, *it, facelocalDim);
						if (bctype == BoundaryConditions2p2c::saturation)
						{
							if (subindexi == -1) std::cout << "boundary conditions for 1p domain is not allowed! index: "<< indexi << " position: " << faceglobal[0] << " ," << faceglobal[1] << std::endl;
							satBound = problem.gS(faceglobal, *it, facelocalDim);
							satFlash(satBound, pressBound, T, poroI, C1Bound, C2Bound, Xw1Bound, Xn1Bound);
						}
						if (bctype == BoundaryConditions2p2c::concentration)
						{
							double Z1Bound = problem.gZ(faceglobal, *it, facelocalDim);
							flashCalculation(Z1Bound, pressBound, T, poroI, satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound);
							if (subindexi == -1 && satBound != 1.) std::cout << "gZ in 1p domain too low, gas appears! index: "<< indexi << " position: " << faceglobal[0] << " ," << faceglobal[1] << std::endl;
						}

						// phase densities on boundary
						double rho_w_Bound = problem.liquidPhase.density(T, pressBound, (1. - Xw1Bound));
						double rho_n_Bound = problem.gasPhase.density(T, pressBound, Xn1Bound);

						// fractional flow factors on boundary
						double viscosityG = problem.gasPhase.viscosity(T, pressBound, Xn1Bound);
						double viscosityL = problem.liquidPhase.viscosity(T, pressBound, 1. - Xw1Bound);

						if (subindexi == -1) // cell does not lie in subdomain
						{
							// get saturation value at cell center
							double satI = 1.;

							// fractional flow factors on boundary
							std::vector<double> kr = problem.materialLaw.kr(satBound, global, *it, local, T);
							double lambda = kr[0] / viscosityL + kr[1] / viscosityG;
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
							factorC2 =
								 velocityJI * (1. - Xw1Bound) * rho_w_Bound * numFlux(satBound, satI, fwBound, fwI)
								- velocityIJ * problem.variables.totalConcentration[indexi + subindexset.size(0)] / poroI * numFlux(satI, satBound, fwI, fwBound);
						}
						else // cell lies in subdomain
						{
							// total mobility and fractional flow factors neighbor cell
							std::vector<double> kr = problem.materialLaw.kr(satBound, global, *it, local, T);
							double lambda = kr[0] / viscosityL + kr[1] / viscosityG;
							double fwBound = kr[0] / viscosityL / lambda;
							double fnBound = kr[1] / viscosityG / lambda;

							factorC1 =
								+ velocityJI * Xw1Bound * rho_w_Bound * numFlux(satBound, satI, fwBound, fwI)
								- velocityIJ * Xw1_I * rho_w_I * numFlux(satI, satBound, fwI, fwBound)
								+ velocityJI * Xn1Bound * rho_n_Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
								- velocityIJ * Xn1_I * rho_n_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound);
							factorC2 =
								 velocityJI * (1. - Xw1Bound) * rho_w_Bound * numFlux(satBound, satI, fwBound, fwI)
								- velocityIJ * (1. - Xw1_I) * rho_w_I * numFlux(satI, satBound, fwI, fwBound)
								+ velocityJI * (1. - Xn1Bound) * rho_n_Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
								- velocityIJ * (1. - Xn1_I) * rho_n_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound);
						}

						// for timestep control
						double coeffW = fwI / satI;
						if(isinf(coeffW) || isnan(coeffW)) coeffW = 0;
						double coeffN = fnI / (1-satI);
						if(isinf(coeffN) || isnan(coeffN)) coeffN = 0;
						factor = (velocityJI - velocityIJ) * std::max(std::max(coeffW, coeffN),1.);

					}
					else if (pressBCtype == BoundaryConditions::neumann)
					{
						FieldVector<RT,2> J = problem.J(faceglobal, *it, facelocalDim);
						double faceVol = integrationOuterNormal.two_norm();
						factorC1 = J[0] * faceVol / volume;
						factorC2 = J[1] * faceVol / volume;
					}

					else DUNE_THROW(NotImplemented, "there is no process boundary condition implemented");
				}

				// add to update vector
				updateVec[indexi + subindexset.size(0)] += factorC2;
				if (subindexi != -1)
					updateVec[subindexi] += factorC1;

				// TODO remove DEBUG --->
				factorcheck[indexi][numberInSelf] = factorC2;
				//<--- DEBUG

				// for time step calculation
				if (factor>=0)
					sumfactor += factor;
				else
					sumfactor2 += (-factor);
			} // end all intersections

			// compute timestep restriction
			sumfactor = std::max(sumfactor,sumfactor2) / poroI;
			sumDiff = std::max(sumDiff,sumDiff2);
			sumfactor = std::max(sumfactor,100*sumDiff);
			if ( 1./sumfactor < dt)
			{
				dt = 1./sumfactor;
				which= indexi;
			}

		} // end grid traversal
	//		printvector(std::cout,updateVec,"update","row");

		//TODO remove DEBUG --->
		for (Iterator hit = grid.template lbegin<0>(level_); hit != endit; ++hit)
		{
			int indexi = indexset.index(*hit);
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*hit);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*hit); is!=endit; ++is)
			{
				if (is.neighbor())
				{
					int indexj = indexset.index(*is.outside());
					double v1 = fabs(factorcheck[indexi][is.numberInSelf()][0]);
					double v2 = fabs(factorcheck[indexj][is.numberInNeighbor()][0]);
					double dev;
					if (v1 + v2 != 0.) dev = fabs( 2* (v1-v2)/(v1+v2));
					else dev = 0.;
					if (dev > 1e-10) std::cout<< "mass leak in cell " << indexi << " face " << is.numberInSelf() <<" " << dev << std::endl;
				}
			}
		}
		// <---DEBUG
		std::cout << " timestep restricting cell: " << which << std::endl;
		return which;
  } // end function "update"

  template<class G, class RT>
  void Multiphysics2p2c<G,RT>::flashCalculation(double Z1, double p, double temp, double poro, double& sat, double& C1, double& C2, double& Xw1, double& Xn1)
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

    double nu2 = 0;

    if (Z1 > Xn1 && Z1 < Xw1)
    	nu2 = -((K1-1)*Z1 + (K2-1)*(1-Z1)) / (K1-1) / (K2-1);

    else if (Z1 < Xn1)
    {
    	nu2 = 1;
    	Xn1 = Z1;
    }
    else if (Z1 > Xw1)
    {
    	nu2 = 0;
    	Xw1 = Z1;
		}

    double rho_w = problem.liquidPhase.density(temp, p, 1. - Xw1);
    double rho_n = problem.gasPhase.density(temp, p, Xn1);

    sat = (1-nu2) / rho_w;
    sat /= ((1-nu2)/rho_w + nu2/rho_n);

    C1 = poro * (Xw1 * sat * rho_w + Xn1 * (1-sat) * rho_n);
    C2 = poro * ((1. - Xw1) * sat * rho_w + (1. - Xn1) * (1-sat) * rho_n);

  } // end function flashCalculation

  template<class G, class RT>
  void Multiphysics2p2c<G,RT>::satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Xw1, double& Xn1)
  {
  	if (sat <= 0 || sat >= 1)
	  DUNE_THROW(RangeError,
	  	"Multiphysics2p2c :: saturation initial and boundary conditions may not equal zero or one!");
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

    C1  = poro* (sat * Xw1 * rho_w + (1-sat) * Xn1 * rho_n);
    C2  = poro* (sat * (1. - Xw1) * rho_w + (1-sat) * (1. - Xn1) * rho_n);
  }

  template<class G, class RT>
  void Multiphysics2p2c<G,RT>::postupdate(double t, double dt)
  {
  	int size = elementmapper.size();
  	// iterate through leaf grid an evaluate c0 at cell center
    SubIterator endit = subgrid.template lend<0>(level_);
    for (SubIterator it = subgrid.template lbegin<0>(level_); it != endit; ++it)
		{
    	EntityPointer epit = subgrid.template getHostEntity<0>(*it);
    	int indexi = indexset.index(*epit);
    	int subindex = subindexset.index(*it);

    	// get cell geometry informations
			GeometryType gt = it->geometry().type(); //geometry type
			const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
			FieldVector<ct,dimworld> global = it->geometry().global(local); // cell center in global coordinates

			double poro = problem.porosity(global, *epit, local);
			double rho_l = problem.liquidPhase.density(T, problem.variables.pressure[indexi][0],0.);
			double rho_g = problem.gasPhase.density(T, problem.variables.pressure[indexi][0],0.);

    	double Z1 = problem.variables.totalConcentration[subindex] / (problem.variables.totalConcentration[subindex] + problem.variables.totalConcentration[subindexset.size(0)+indexi]);
    	double C1 = problem.variables.totalConcentration[subindex][0];
    	double C2 = problem.variables.totalConcentration[subindexset.size(0)+indexi][0];
			flashCalculation(Z1, problem.variables.pressure[indexi], T, poro, problem.variables.saturation[subindex][0], C1, C2, problem.variables.wet_X1[subindex][0], problem.variables.nonwet_X1[subindex][0]);

			double nuw = problem.variables.saturation[subindex] * rho_l / (problem.variables.saturation[subindex] * rho_l + (1-problem.variables.saturation[subindex]) * rho_g);
			double massw = (problem.variables.totalConcentration[subindex][0] + problem.variables.totalConcentration[subindexset.size(0)+indexi][0]) * nuw;
			double massn = (problem.variables.totalConcentration[subindex][0] + problem.variables.totalConcentration[subindexset.size(0)+indexi][0]) * (1-nuw);
			double vol = massw / rho_l + massn / rho_g;
			problem.variables.volErr[subindex] = (vol - poro) / dt;
		}
		timestep = dt;
  }

}//end namespace Dune

#endif /*MULTIPHYSICS2P2C_HH_*/
