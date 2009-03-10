// $Id:$

#ifndef DUNE_FVSHALLOWWATER_HH
#define DUNE_FVSHALLOWWATER_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/bvector.hh>
#include "dumux/shallowwater/shallowwater.hh"
#include "dumux/shallowwater/shallownumericalflux.hh"
#include "dumux/shallowwater/shallowproblemplain.hh"
#include "dumux/shallowwater/shallowvariableclass.hh"

namespace Dune
{
//! \ingroup transport                                                                                                                                                                                                                                                                                                                                                                                            
//! The finite volume model for the solution of the transport equation
template<class G, class DT, class VC> class FVShallowWater :
	public ShallowWater< G, DT, VC>
{
	template<int dim> struct ElementLayout
	{
		bool contains(Dune::GeometryType gt)
		{
			return gt.dim() == dim;
		}
	};

	enum
	{	dim = G::dimension};
	enum
	{	dimworld = G::dimensionworld};

	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::template Codim<0>::EntityPointer EntityPointer;

	typedef typename G::LeafGridView GV;
	typedef typename GV::IndexSet IS;
	typedef typename GV::template Codim<0>::Iterator Iterator;
	typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator
			IntersectionIterator;

	typedef Dune::FieldVector<DT,dim> SlopeType;
	typedef Dune::FieldVector<DT,dim+1> SystemType;
	typedef Dune::BlockVector<Dune::FieldVector<DT,dim+1> > SolutionType;
	typedef Dune::FieldVector<DT,dim> LocalPosition;
	typedef Dune::FieldVector<DT,dimworld> GlobalPosition;

public:
	typedef Dune::BlockVector<Dune::FieldVector<DT,dim+1> > RepresentationType;

	int update(const DT t, DT& dt, SolutionType& updateVec, DT& cFLFac);

	void initialize();

	FVShallowWater(G& g, ShallowProblemBase<G, DT, VC>& problem,
			NumericalFlux<G,DT>& numFl = *(new FirstOrderUpwind<G,DT>)) :
		ShallowWater<G, DT, VC>(g, problem),
				elementmapper(g, g.leafIndexSet()), gridview(g.leafView()),
				indexset(gridview.indexSet()), numFlux(numFl)
	{
	}

private:
	EM elementmapper;
	const GV& gridview;
	const IS& indexset;
	NumericalFlux<G,DT>& numFlux;

};

template<class G, class DT, class VC> int FVShallowWater<G, DT, VC>::update(
		const DT t, DT& dt, SolutionType& updateVec, DT& cFLFac = 1)
{
	// initialize dt very large, why?
	dt = 1e100;
	updateVec = 0;
	// initialize and declare variables

	// compute update vector
	Iterator eendit = this->grid.template leafend<0>();
	for (Iterator it = this->grid.template leafbegin<0>(); it != eendit; ++it)
	{
		DT dist;
		DT bottomElevationI;
		DT bottomElevationJ;
		DT wDepthI;
		DT wDepthJ;
		DT wDepthFace;
		DT wDepthFaceI;
		DT wDepthFaceJ;
		SlopeType velI;
		SlopeType velJ;
		SlopeType bottomSlopeTerm(0); //slope term resulting from divergence form
		SlopeType bottomSlope(0); //real bottomslope
		SystemType bottomSlopeVector(0);
		SystemType flux(0);
		SlopeType divergenceTerm(0);
		SystemType summedFluxes(0);

		// cell geometry type
		Dune::GeometryType gt = it->geometry().type();

		// cell center in reference element
		const LocalPosition &localPos =
				Dune::ReferenceElements<DT,dim>::general(gt).position(0, 0);

		GlobalPosition globalPos = it->geometry().global(localPos);

		// cell volume, assume linear map here
		double volume = it->geometry().integrationElement(localPos)
				*Dune::ReferenceElements<DT,dim>::general(gt).volume();

		// cell index
		int globalIdxI = elementmapper.map(*it);

		//get bottomslopevector for entity
		bottomSlope =this->problem.surface.calcBottomSlopes(globalPos, *it,
				localPos);

		// get waterdepth at cell center
		wDepthI = this->problem.variables.wDepth[globalIdxI];

		// get velocity at cell center
		velI = this->problem.variables.velocity[globalIdxI];

		// run through all intersections with neighbors and boundary
		IntersectionIterator endit =
				IntersectionIteratorGetter<G, LeafTag>::end(*it);
		for (IntersectionIterator is =
				IntersectionIteratorGetter<G, LeafTag>::begin(*it); is !=endit; ++is)
		{

			// local number of facet
			int numberInSelf = is->numberInSelf();

			// get geometry type of face
			Dune::GeometryType gtf = is->intersectionSelfLocal().type();

			// center in face's reference element
			const Dune::FieldVector<DT,dim-1>& faceLocalPos =
					Dune::ReferenceElements<DT,dim-1>::general(gtf).position(0, 0);

			// center of face inside volume reference element
			const LocalPosition& faceLocalDim =
					Dune::ReferenceElements<DT,dim>::general(gtf).position(is->numberInSelf(), 1);

			//get normal vector of face
			Dune::FieldVector<DT,dimworld> nVec = is->outerNormal(faceLocalPos);

			// get normal vector scaled with volume
			Dune::FieldVector<DT,dimworld> nVecScaled =
					is->integrationOuterNormal(faceLocalPos);
			nVecScaled*=Dune::ReferenceElements<DT,dim-1>::general(gtf).volume();
			

			// handle interior face
			if (is->neighbor())
			{
				// access neighbor
				EntityPointer outside = is->outside();
				int globalIdxJ = elementmapper.map(*outside);

				Dune::GeometryType nbgt = outside->geometry().type();
				const LocalPosition& nbLocalPos =
						Dune::ReferenceElements<DT,dim>::general(nbgt).position(0, 0);

				// neighbor cell center in global coordinates
				Dune::FieldVector<DT,dimworld> nbGlobalPos = outside->geometry().global(nbLocalPos);

				// distance vector between barycenters
				Dune::FieldVector<DT,dimworld> distVec = globalPos - nbGlobalPos;

				//compute distance between cell centers
				//DT dist = distVec.two_norm();

				// get waterdepth at neighbor cell center
				wDepthJ = this->problem.variables.wDepth[globalIdxJ];

				// get velocity at neighbor cell center
				velJ = this->problem.variables.velocity[globalIdxJ];

				//fluxVector for a given scheme (Upwindor different will be catched from numericalflux)

				flux = numFlux(velI, velJ, wDepthI, wDepthJ,
						nVecScaled, nVec);
			
				std::cout<<"cell "<<globalIdxI<<"flux_interior of face "
					<<numberInSelf<<"=" <<flux<<std::endl;

			}

			// handle boundary face
			if (is->boundary())
			{
				// center of face in global coordinates
				GlobalPosition faceGlobalPos = is->intersectionGlobal().global(faceLocalPos);

				// distance vector between barycenters
				Dune::FieldVector<DT,dimworld> distVec = globalPos
						- faceGlobalPos;

				//compute distance between cell centers
				//DT dist = distVec.two_norm();

				//get boundary type
				BoundaryConditions::Flags bctype = this->problem.bctype(
						faceGlobalPos, *it, faceLocalDim);


				if (bctype == BoundaryConditions::dirichlet)
				{
					// get waterdepth at boundary
					wDepthFace = this->problem.dirichlet(faceGlobalPos, *it,
							faceLocalDim);
					velFace(0);

					DT conti = velFace * wDepthFace;
					DT xMomentum = velFace * wDepthFace*velFace;
					xMomentum += gravity*wDepthFace*0.5;
							
					flux[0] = conti;
					flux[1] = xMomentum;
					flux[2] = yMomentum;

					
					//std::cout<<"flux_dirichlet of face "<<numberInSelf<<"="
					//<<flux<<std::endl;


				}
				if (bctype == BoundaryConditions::neumann) //no flow condition at side walls

				{

					flux = this->problem.neumann(faceGlobalPos, *it,
					faceLocalDim);

					//std::cout<<"flux_neumann of face "<<numberInSelf<<"="<<flux
					//<<std::endl;
				}
			}
			summedFluxes += flux;

		}
			
		bottomSlopeVector[0]=0;
			bottomSlopeVector[1]=9.81;
			bottomSlopeVector[1]*=bottomSlope[0];
			bottomSlopeVector[1]*=(wDepthI);
			bottomSlopeVector[2]=9.81;
			bottomSlopeVector[2]*=bottomSlope[1];
			bottomSlopeVector[2]*=(wDepthI);

		//std::cout<<"BottomSlope Vector of cell="<<globalIdxI<<"= "
		//<<bottomSlopeVector<<std::endl;

		//Quelle einbinden und in Systemvektor konvertieren
		DT sourceTerm = this->problem.setSource(globalPos, *it, localPos);

		SystemType sourceTermVector(0);
		sourceTermVector[0]=sourceTerm;
		sourceTermVector[1]=0;
		sourceTermVector[2]=0;
		sourceTermVector *= volume;

		//std::cout<<"Source Term Vector of cell="<<globalIdxI<<"= "
		//<<sourceTermVector<<std::endl;

		updateVec[globalIdxI]=summedFluxes;		
		updateVec[globalIdxI]-= sourceTermVector;
		updateVec[globalIdxI]-=bottomSlopeVector;
		updateVec[globalIdxI]/= (-volume);
		
		//std::cout<<"update for cell"<<globalIdxI<<" = "<<updateVec[globalIdxI]
		//<<std::endl;

		//add cfl-criterium
		//calculate timestep for every cell and take the minimum

	}

	// end grid traversal

	return 0;
}

template<class G, class DT, class VC> void FVShallowWater<G, DT, VC>::initialize()
{

	Iterator eendit = this->grid.template leafend<0>();
	for (Iterator it = this->grid.template leafbegin<0>(); it != eendit; ++it)
	{
		// get geometry type
		Dune::GeometryType gt = it->geometry().type();

		// get cell center in reference element
		const LocalPosition &localPos =
				Dune::ReferenceElements<DT,dim>::general(gt).position(0, 0);

		// get global coordinate of cell center
		GlobalPosition globalPos = it->geometry().global(localPos);

		int globalIdx = elementmapper.map(*it);

		DT initialWaterDepth=this->problem.setInitWDepth(globalPos, *it,
				localPos);
		SlopeType initialVelocity=this->problem.setInitVel(globalPos, *it,
				localPos);

		// initialize cell values
		this->problem.variables.wDepth[globalIdx] = initialWaterDepth;
		this->problem.variables.velocity[globalIdx] = initialVelocity;
		this->problem.variables.globalSolution[globalIdx][0]
				= initialWaterDepth;
		this->problem.variables.globalSolution[globalIdx][1]
				= initialWaterDepth*initialVelocity[0];
		this->problem.variables.globalSolution[globalIdx][1]
				= initialWaterDepth*initialVelocity[1];
	}
	return;
}

}
#endif

// end namespace


