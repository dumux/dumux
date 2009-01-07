// $Id: capillarydiffusion.hh 899 2008-12-08 19:02:34Z lauser $

#ifndef DUNE_GRAVITYPART_HH
#define DUNE_GRAVITYPART_HH

#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/diffusion/diffusionproblem.hh"

//! \ingroup transport
//! \defgroup diffPart Diffusive transport
/**
* @file
* @brief  Base class for defining the diffusive part of an advection-diffusion equation
* @author Bernd Flemisch, last changed by Markus Wolff
*/
namespace Dune
{
/*!\ingroup diffPart
* @brief  Base class for defining the diffusive part of an advection-diffusion equation
*/
template<class G, class RT, class VC>
class GravityPart : public DiffusivePart<G,RT>
{
	enum{dim = G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
	typedef Dune::FieldVector<RT, dim> FieldVector;
	typedef BlockVector< Dune::FieldVector<RT,1> > SatType;

public:
	virtual FieldVector operator() (const Entity& entity, const int numberInSelf,
			const RT satIntersection, const FieldVector& satGradient, const RT time,
			const RT satI, const RT satJ) const
{
		// cell geometry type
		GeometryType gt = entity.geometry().type();

		// cell center in reference element
		const FieldVector& local = ReferenceElements<RT,dim>::general(gt).position(0,0);

		// get global coordinate of cell center
		FieldVector global = entity.geometry().global(local);

		// get absolute permeability of cell
		FieldMatrix<RT,dim,dim> K(problem.K(global,entity,local));

		IntersectionIterator endis = entity.ilevelend();
		IntersectionIterator is = entity.ilevelbegin();
		for (; is != endis; ++is)
		{
			if(is->numberInSelf() == numberInSelf)
				break;
		}

		// get geometry type of face
		GeometryType gtf = is->intersectionSelfLocal().type();

		// center in face's reference element
		const Dune::FieldVector<RT,dim-1>& facelocal = ReferenceElements<RT,dim-1>::general(gtf).position(0,0);

		FieldVector unitOuterNormal = is->unitOuterNormal(facelocal);
		//std::cout<<"unitOuterNormaldiff"<<unitOuterNormal<<std::endl;

		FieldVector nbglobal(-1e100);
		if (is->neighbor()) {
			// access neighbor
			EntityPointer outside = is->outside();

			// compute factor in neighbor
			GeometryType nbgt = outside->geometry().type();
			const FieldVector& nblocal = ReferenceElements<RT,dim>::general(nbgt).position(0,0);

			// neighbor cell center in global coordinates
			nbglobal = outside->geometry().global(nblocal);

			// take arithmetic average of absolute permeability
			K += problem.K(nbglobal, *outside, nblocal);
			K *= 0.5;
		}

		// set result to (rho_w - rho_n)*g
		FieldVector helpresult = gravity;

		// set result to K*(rho_w - rho_n)*g
		FieldVector result(0);
		K.umv(helpresult, result);

		//get lambda_bar = lambda_n*f_w
		double mobbarI=constRel.mobN(1-satI)*constRel.fractionalW(satI);
		double mobbarJ=constRel.mobN(1-satJ)*constRel.fractionalW(satJ);

		// set result to f_w*lambda_n*K*(rho_w - rho_n)*g
		//result *= (mobbarI+mobbarJ)*0.5;
		if (global[dim-1] > nbglobal[dim-1])
			result *= mobbarI;
		else
			result *= mobbarJ;

		return result;
}

	GravityPart (DiffusionProblem<G, RT, VC>& prob)
	: problem(prob), constRel(problem.materialLaw), wettingPhase(constRel.wettingPhase),
	nonwettingPhase(constRel.nonwettingPhase)
	{
		double rhoDiff = wettingPhase.density() - nonwettingPhase.density();
		gravity = problem.gravity();
		gravity *= -rhoDiff;
	}

private:
	DiffusionProblem<G, RT, VC>& problem;
	TwoPhaseRelations& constRel;
	const Medium& wettingPhase;
	const Medium& nonwettingPhase;
	FieldVector gravity;
};
}

#endif
