// $Id$

#ifndef DUNE_FVDIFFUSION_DEPRECATED_HH
#define DUNE_FVDIFFUSION_DEPRECATED_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/diffusion/diffusion_deprecated.hh"
#include "dumux/pardiso/pardiso.hh"
//#include "dumux/diffusion/problems/uniformproblem.hh"
//#include "dumux/transport/problems/simpleproblem.hh"

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Bernd Flemisch, Jochen Fritz; last changed by Markus Wolff
 */

namespace Dune {
//! \ingroup diffusion
//! Finite Volume Diffusion Model
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
template<class G, class RT, class VC> class FVDiffusion :
	public Diffusion< G, RT, VC> {
	template<int dim> struct ElementLayout {
		bool contains(GeometryType gt) {
			return gt.dim() == dim;
		}
	};

	enum {dim = G::dimension};
	enum {dimworld = G::dimensionworld};

	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::LevelGridView GV;
    typedef typename GV::IndexSet IS;
	typedef typename GV::template Codim<0>::Iterator Iterator;
	typedef typename G::template Codim<0>::HierarchicIterator
			HierarchicIterator;
	typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator
			IntersectionIterator;
	typedef typename G::ctype ct;
	typedef FieldMatrix<double,1,1> MB;
	typedef BCRSMatrix<MB> MatrixType;
	typedef FieldVector<double, 1> VB;
	typedef BlockVector<VB> Vector;
	typedef FieldVector<double,dim> R2;
	typedef BlockVector<R2> R3;
	typedef BlockVector<R3> LocVelType;

public:
	typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> >
			VelType;
	typedef BlockVector< FieldVector<RT,1> > RepresentationType;

	void assemble(const RT t);

	void solve();

	void pressure(const RT t=0) {
		assemble(t);
		solve();
		return;
	}

	void initializeMatrix();

	FVDiffusion(G& g, DiffusionProblem<G, RT, VC>& prob, int lev = -1, bool calcPressure = false) :
		Diffusion<G, RT, VC>(g, prob, lev == -1 ? g.maxLevel() : lev),
				elementmapper(g, g.levelIndexSet(this->level())), gridview(g.levelView(this->level())),
				indexset(gridview.indexSet()), A(g.size(
						this->level(), 0), g.size(this->level(), 0), (2*dim+1)
						*g.size(this->level(), 0), BCRSMatrix<MB>::random),
				f(g.size(this->level(), 0)), solverName_("BiCGSTAB"),
				preconditionerName_("SeqILU0") {
		initializeMatrix();
	}

	FVDiffusion(G& g, DiffusionProblem<G, RT, VC>& prob, std::string solverName,
			std::string preconditionerName, int lev = -1, bool calcPressure = false) :
		Diffusion<G, RT, VC>(g, prob, lev == -1 ? g.maxLevel() : lev),
				elementmapper(g, g.levelIndexSet(this->level())), gridview(g.levelView(this->level())),
				indexset(gridview.indexSet()), A(g.size(
						this->level(), 0), g.size(this->level(), 0), (2*dim+1)
						*g.size(this->level(), 0), BCRSMatrix<MB>::random),
				f(g.size(this->level(), 0)), solverName_(solverName),
				preconditionerName_(preconditionerName) {
		initializeMatrix();
	}

	EM elementmapper;
	const GV& gridview;
	const IS& indexset;

private:

	MatrixType A;
	RepresentationType f;
	std::string solverName_;
	std::string preconditionerName_;
};

template<class G, class RT, class VC> void FVDiffusion<G, RT, VC>::initializeMatrix() {
	// determine matrix row sizes
	Iterator eendit = this->diffproblem.variables.grid.template lend<0>(this->level());
	for (Iterator it = this->diffproblem.variables.grid.template lbegin<0>(this->level()); it != eendit; ++it) {
		// cell index
		int indexi = elementmapper.map(*it);

		// initialize row size
		int rowSize = 1;

		// run through all intersections with neighbors
		IntersectionIterator
				endit = IntersectionIteratorGetter<G, LevelTag>::end(*it);
		for (IntersectionIterator
				is = IntersectionIteratorGetter<G, LevelTag>::begin(*it); is
				!=endit; ++is)
			if (is->neighbor())
				rowSize++;
		A.setrowsize(indexi, rowSize);
	}
	A.endrowsizes();

	// determine position of matrix entries
	for (Iterator it = this->diffproblem.variables.grid.template lbegin<0>(this->level()); it != eendit; ++it) {
		// cell index
		int indexi = elementmapper.map(*it);

		// add diagonal index
		A.addindex(indexi, indexi);

		// run through all intersections with neighbors
		IntersectionIterator
				endit = IntersectionIteratorGetter<G, LevelTag>::end(*it);
		for (IntersectionIterator
				is = IntersectionIteratorGetter<G, LevelTag>::begin(*it); is
				!=endit; ++is)
			if (is->neighbor()) {
				// access neighbor
				EntityPointer outside = is->outside();
				int indexj = elementmapper.map(*outside);

				// add off diagonal index
				A.addindex(indexi, indexj);
			}
	}
	A.endindices();

	return;
}

template<class G, class RT, class VC> void FVDiffusion<G, RT, VC>::assemble(const RT t=0) {
	// initialization: set matrix A to zero
	A = 0;

	// find out whether gravity effects are relevant
	bool hasGravity = false;
	const FieldVector<ct,dim>& gravity = this->diffproblem.gravity();
	for (int k = 0; k < dim; k++)
		if (gravity[k] != 0)
			hasGravity = true;

	Iterator eendit = this->diffproblem.variables.grid.template lend<0>(this->level());
	for (Iterator it = this->diffproblem.variables.grid.template lbegin<0>(this->level()); it != eendit; ++it) {
		// cell geometry type
		GeometryType gt = it->geometry().type();

		// cell center in reference element
		const FieldVector<ct,dim>&local = ReferenceElements<ct,dim>::general(gt).position(0, 0);

		// get global coordinate of cell center
		FieldVector<ct,dim> global = it->geometry().global(local);

		// cell index
		int indexi = elementmapper.map(*it);

		// cell volume, assume linear map here
		double volume = it->geometry().integrationElement(local)*ReferenceElements<ct,dim>::general(gt).volume();

		// set right side to zero
		f[indexi] = volume*this->diffproblem.source(global, *it, local);

		// get absolute permeability
		FieldMatrix<ct,dim,dim> Ki(this->diffproblem.K(global, *it, local));

		//compute total mobility
		double lambdaI = 0, fractionalWI = 0;
		double sati = this->diffproblem.variables.saturation[indexi];
		lambdaI = this->diffproblem.materialLaw.mobTotal(sati);
		if (hasGravity)
			fractionalWI = this->diffproblem.materialLaw.fractionalW(sati);

		IntersectionIterator
				endit = IntersectionIteratorGetter<G, LevelTag>::end(*it);
		for (IntersectionIterator
				is = IntersectionIteratorGetter<G, LevelTag>::begin(*it); is
				!=endit; ++is) {

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
			integrationOuterNormal
			*= ReferenceElements<ct,dim-1>::general(gtf).volume();

			// get face volume
			double faceVol = 1;
			switch (G::dimension) {
				case 1: break;
				default: faceVol = is->intersectionGlobal().volume();
				break;
			}

			// compute directed permeability vector Ki.n
			FieldVector<ct,dim> Kni(0);
			Ki.umv(unitOuterNormal, Kni);

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

				// get absolute permeability
				FieldMatrix<ct,dim,dim> Kj(this->diffproblem.K(nbglobal, *outside, nblocal));

				// compute vectorized permeabilities
				FieldVector<ct,dim> Knj(0);
				Kj.umv(unitOuterNormal, Knj);
				double K_n_i = Kni * unitOuterNormal;
				double K_n_j = Knj * unitOuterNormal;
				double Kn = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
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
				double satj = this->diffproblem.variables.saturation[indexj];

				lambdaJ = this->diffproblem.materialLaw.mobTotal(satj);
				if (hasGravity)
				fractionalWJ = this->diffproblem.materialLaw.fractionalW(satj);

				// compute averaged total mobility
				// CAREFUL: Harmonic weightig can generate zero matrix entries,
				// use arithmetic weighting instead:
				double fractionalW = 0;
				double lambda = 0.5*(lambdaI + lambdaJ);
				if (hasGravity)
				fractionalW = 0.5*(fractionalWI + fractionalWJ);

				// update diagonal entry
				double entry = fabs(lambda*faceVol*(K*distVec)/(dist*dist));
				A[indexi][indexi] += entry;

				// set off-diagonal entry
				A[indexi][indexj] = -entry;

				if (hasGravity) {
					double factor = fractionalW*(this->diffproblem.materialLaw.wettingPhase.density())
					+ (1 - fractionalW)*(this->diffproblem.materialLaw.nonwettingPhase.density());
					f[indexi] -= factor*lambda*faceVol*(K*gravity);
				}

				if (this->diffproblem.capillarity) {
					// calculate saturation gradient
					FieldVector<ct,dim> satGradient = distVec;
					satGradient *= (satj - sati)/(dist*dist);

					// arithmetic average of the permeability
					K = ((Kni + Knj) *= 0.5);

					// capillary pressure w.r.t. saturation
					double pCI = this->diffproblem.materialLaw.pC(sati);
					double pCJ = this->diffproblem.materialLaw.pC(satj);

					// mobility of the nonwetting phase
					double lambdaN = 0.5*(this->diffproblem.materialLaw.mobN(1 - sati)
							+ this->diffproblem.materialLaw.mobN(1 - satj));

					// calculate capillary pressure gradient
					FieldVector<ct,dim> pCGradient = distVec;
					pCGradient *= -(pCJ - pCI)/(dist*dist);

					f[indexi] += lambdaN*faceVol*(K*pCGradient);
				}
			}
			// boundary face
			else
			{
				// center of face in global coordinates
				FieldVector<ct,dimworld>
				faceglobal = is->intersectionGlobal().global(facelocal);

				// compute total mobility
				double fractionalW = 1.;
				double lambda = lambdaI;
				if (hasGravity) fractionalW = fractionalWI;

				//get boundary condition for boundary face center
				BoundaryConditions::Flags bctype = this->diffproblem.bctype(faceglobal, *it, facelocalDim);
				if (bctype == BoundaryConditions::dirichlet)
				{
					FieldVector<ct,dimworld> distVec(global - faceglobal);
					double dist = distVec.two_norm();
					A[indexi][indexi] -= lambda*faceVol*(Kni*distVec)/(dist*dist);
					double g = this->diffproblem.dirichletPress(faceglobal, *it, facelocalDim);
					f[indexi] -= lambda*faceVol*g*(Kni*distVec)/(dist*dist);

					if (hasGravity) {
						double factor = fractionalW*(this->diffproblem.materialLaw.wettingPhase.density())
						+ (1 - fractionalW)*(this->diffproblem.materialLaw.nonwettingPhase.density());
						f[indexi] -= factor*lambda*faceVol*(Kni*gravity);
					}
					if (this->diffproblem.capillarity) {
						double satj = this->diffproblem.dirichletSat(faceglobal, *it, facelocalDim);

						// distance vector between barycenters
						FieldVector<ct,dimworld>
						distVec = global - faceglobal;

						// compute distance between cell centers
						double dist = distVec.two_norm();

						// calculate saturation gradient
						FieldVector<ct,dim> satGradient = distVec;
						satGradient *= (satj - sati)/(dist*dist);

						// capillary pressure w.r.t. saturation
						double pCI = this->diffproblem.materialLaw.pC(sati);
						double pCJ = this->diffproblem.materialLaw.pC(satj);

						// mobility of the nonwetting phase
						double lambdaN = 0.5*(this->diffproblem.materialLaw.mobN(1 - sati)
								+ this->diffproblem.materialLaw.mobN(1 - satj));

						// calculate capillary pressure gradient
						FieldVector<ct,dim> pCGradient = distVec;
						pCGradient *= -(pCJ - pCI)/(dist*dist);

						f[indexi] -= lambdaN*faceVol*(Kni*pCGradient);
					}
				}
				else
				{
					double J = this->diffproblem.neumannPress(faceglobal, *it, facelocalDim);
					f[indexi] -= faceVol*J;
				}
			}
		}
		// end all intersections
	} // end grid traversal
	return;
}

template<class G, class RT, class VC> void FVDiffusion<G, RT, VC>::solve() {
	MatrixAdapter<MatrixType,Vector,Vector> op(A);
	InverseOperatorResult r;

	if (preconditionerName_ == "SeqILU0") {
		SeqILU0<MatrixType,Vector,Vector> preconditioner(A, 1.0);
		if (solverName_ == "CG") {
			CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
			solver.apply((this->diffproblem.variables.pressure), f, r);
		} else if (solverName_ == "BiCGSTAB") {
			BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
			solver.apply((this->diffproblem.variables.pressure), f, r);
		} else
			DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
					<< preconditionerName_<< " and "<< solverName_ << ".");
	} else if (preconditionerName_ == "SeqPardiso") {
		SeqPardiso<MatrixType,Vector,Vector> preconditioner(A);
		if (solverName_ == "Loop") {
			LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
			solver.apply(this->diffproblem.variables.pressure, f, r);
		} else
			DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
					<< preconditionerName_<< " and "<< solverName_ << ".");
	} else
		DUNE_THROW(NotImplemented, "FVDiffusion :: solve : preconditioner "
				<< preconditionerName_ << ".");
//	printmatrix(std::cout, A, "global stiffness matrix", "row", 11, 3);
//	printvector(std::cout, f, "right hand side", "row", 200, 1, 3);
//	printvector(std::cout, (this->diffproblem.variables.pressure), "pressure", "row", 200, 1, 3);
	return;
}

}
#endif
