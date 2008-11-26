// $Id$
#ifndef DUNE_FVDIFFUSIONVELOCITYUPSCALED_HH
#define DUNE_FVDIFFUSIONVELOCITYUPSCALED_HH

#include "dumux/diffusion/fv/fvdiffusion.hh"
#include <../../../../dune-subgrid/subgrid/subgrid.hh>

namespace Dune
{

template<class G, class RT, class VC>
class FVDiffusionVelocityUpscaled: public FVDiffusion<G, RT, VC>
{

	enum
	{
		dim = G::dimension
	};
	enum
	{
		dimworld = G::dimensionworld
	};

typedef	typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::ctype ct;
	typedef typename G::LevelGridView GV;
	typedef typename GV::IndexSet IS;
	typedef typename GV::template Codim<0>::Iterator Iterator;
	typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator
	IntersectionIterator;
	typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
	typedef typename G::template Codim<0>::EntityPointer EntityPointer;

public:
	FVDiffusionVelocityUpscaled(G& g, FractionalFlowProblem<G, RT, VC>& prob, int lev = -1)
	: FVDiffusion<G,RT,VC>(g, prob, lev)
	{}

	void calcTotalVelocity(const RT t, int lev = -1) const
	{
		if (lev<0)
		lev=this->diffproblem.variables.translevel;
		//		std::cout<<"lev = "<<lev<<std::endl;

		if (this->level()<=lev)
		DUNE_THROW(NotImplemented,"this->level()=<lev");
		const int nFine = (int)pow(2, ((this->level()-lev)*(dim-1)) ); //fine element edges per coarse element edge
		const int blocksize = 2 * dim;
		int hitcount[blocksize];
		const IS& coarseIset = this->grid.levelIndexSet(lev);

		typedef FieldVector<ct,dim> FineScaleVelType; // velocity for fine scale elemnt edge
		typedef BlockVector<FineScaleVelType> FineScaleonCoarseEdgeType; // vector holding velocities of all fine edges situated on one coarse edge
		BlockVector<FineScaleonCoarseEdgeType> fineOnCoarse(blocksize); // vector for all faces of the coarse scale element

		FineScaleonCoarseEdgeType fineVelocity(blocksize);

		for (int i=0; i<blocksize; i++)
		{
			fineOnCoarse[i].resize(nFine);
		}
		fineVelocity=0;
		fineOnCoarse=0;

		Iterator eendit = this->grid.template lend<0>(lev);
		for (Iterator cit = this->grid.template lbegin<0>(lev); cit != eendit; ++cit)
		{
			for (int i=0; i<blocksize; i++) hitcount[i] = 0;
			int coarseIndex = coarseIset.index(*cit);

			HierarchicIterator endit = cit-> hend(this->level());
			for (HierarchicIterator it = cit->hbegin(this->level()); it != endit; ++it)
			{
				//only iterat through difflevel!!!
				if (it->level() != this->level()) continue;

				//get some cell properties
				GeometryType gt = it->geometry().type();
				const FieldVector<ct,dim>&
				local = ReferenceElements<ct,dim>::general(gt).position(0,0);
				FieldVector<ct,dimworld> global = it->geometry().global(local); //global coordinates of cell center
				int indexi = this->elementmapper.map(*it); // index of fine-scale cell


				IntersectionIterator isend = it->ilevelend();
				for (IntersectionIterator is = it->ilevelbegin(); is!=isend; ++is)
				{
					GeometryType gtf = is->intersectionSelfLocal().type();
					const FieldVector<ct,dim-1>& facelocal
					= ReferenceElements<ct,dim-1>::general(gtf).position(0,0);
					FieldVector<ct,dim> faceglobal
					= is->intersectionGlobal().global(facelocal); // global coordinate of face center

					if (is->neighbor())
					{
						// neighbor's properties
						EntityPointer outside = is->outside();
						int indexj = this->elementmapper.map(*outside);// neigbor's fine-scale cell index

						GeometryType nbgt = outside->geometry().type();
						const FieldVector<ct,dim>&
						nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
						FieldVector<ct,dimworld>
						nbglobal = outside->geometry().global(nblocal); // global coordinate of neighbor's cell center


						// check if fine scale face *is lies on any face *cis of the coarse scale element
						IntersectionIterator cisend = cit->ilevelend();
						for (IntersectionIterator cis = cit->ilevelbegin(); cis != cisend; ++cis)
						{
							// checking if fine scale element face *is lies on coarse-scale element face *cis
							GeometryType cgtf = cis->intersectionSelfLocal().type();
							FieldVector<ct,dim-1> cfacelocal = ReferenceElements<ct,dim-1>::general(cgtf).position(0,0);
							bool isInCis;
							if ( cis->neighbor() )
							{
								/* if face *is lies in face *cis, *is must be defined both in the inside and outside
								 element of *cis! for the inside element this is true anyway, because of the hierarchic
								 Iteration. */
								isInCis = cis->outside()->geometry().checkInside( cis->outside()->geometry().local(faceglobal));
							}
							else isInCis = false;
							if ( isInCis )
							{
								int coarseNumberInSelf = cis->numberInSelf();
								int hits = hitcount[coarseNumberInSelf];

								EntityPointer coutside = cis->outside();
								int coarseIndexj = this->diffproblem.variables.transmapper.map(*coutside);

								// get pressure and permeability and total mobility in fine-scale element
								double pressi = this->diffproblem.variables.pressure[indexi];

								FieldMatrix<ct,dim,dim> Ki = this->diffproblem.soil.K(global,*it,local);

								double sati = this->diffproblem.variables.saturation[coarseIndex];
								double lambdaI = this->diffproblem.materialLaw.mobTotal(sati,global,*it,local);

								double faceVol[2*dim];
								// run through all intersections with neighbors and boundary

								// get some face properties
								int numberInSelf = is->numberInSelf();
								switch(G::dimension)
								{
									case 1: faceVol[numberInSelf] = 1;
									default: faceVol[numberInSelf] = is->intersectionGlobal().volume(); // volume of face
								}

								FieldVector<ct,dim> unitOuterNormal
								= is->unitOuterNormal(facelocal); // normal vector of unit length

								// handle interior face

								// get neighbor pressure and permeability
								double pressj = this->diffproblem.variables.pressure[indexj];

								FieldMatrix<ct,dim,dim> Kj = this->diffproblem.soil.K(nbglobal,*outside,nblocal);

								// compute vectorized permeabilities
								FieldVector<ct,dim> Kni(0);
								FieldVector<ct,dim> Knj(0);
								Ki.umv(unitOuterNormal, Kni);
								Kj.umv(unitOuterNormal, Knj);
								// compute permeability normal to intersection and take harmonic mean
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

								// distance between cell centers
								FieldVector<ct,dimworld>
								distVec = global - nbglobal;
								double dist = distVec.two_norm();
								// get averaged total mobility
								double satj = this->diffproblem.variables.saturation[coarseIndexj];
								double lambdaJ = this->diffproblem.materialLaw.mobTotal(satj,nbglobal,*outside,nblocal);
								double mob = 0.5 * (lambdaI + lambdaJ);
								// compute total velocity with Darcy's Law
								fineVelocity[numberInSelf] = K;
								fineVelocity[numberInSelf] *= mob * (pressi-pressj) / dist;
								fineOnCoarse[coarseNumberInSelf][hits] = fineVelocity[numberInSelf];
								hitcount[coarseNumberInSelf] += 1;
							}
						}
					}
					// boundary face

					else
					{
						// check if fine scale face *is lies on any face *cis of the coarse scale element
						IntersectionIterator cisend = cit->ilevelend();
						for (IntersectionIterator cis = cit->ilevelbegin(); cis != cisend; ++cis)
						{
							// checking if fine scale element face *is lies on coarse-scale element face *cis
							GeometryType cgtf = cis->intersectionSelfLocal().type();
							FieldVector<ct,dim-1> cfacelocal = ReferenceElements<ct,dim-1>::general(cgtf).position(0,0);
							bool isInCis;
							if (cis->boundary() && is->boundary())
							{
								/* if face *is lies in face *cis and
								 *cis lies on a boundary, *is must lie on a boundry too!
								 to be sure, that it is the same boundary, check if the normals are the same. */
								isInCis = ( cis->unitOuterNormal(cfacelocal) * is->unitOuterNormal(facelocal) )> 0.999;
							}
							else isInCis = false;
							if ( isInCis )
							{
								int coarseNumberInSelf = cis->numberInSelf();
								int hits = hitcount[coarseNumberInSelf];

								// get pressure and permeability and total mobility in fine-scale element
								double pressi = this->diffproblem.variables.pressure[indexi];

								FieldMatrix<ct,dim,dim> Ki = this->diffproblem.soil.K(global,*it,local);

								double sati = this->diffproblem.variables.saturation[coarseIndex];
								double lambdaI = this->diffproblem.materialLaw.mobTotal(sati,global,*it,local);

								double faceVol[2*dim];
								// run through all intersections with neighbors and boundary

								// get some face properties
								int numberInSelf = is->numberInSelf();
								switch(G::dimension)
								{
									case 1: faceVol[numberInSelf] = 1;
									default: faceVol[numberInSelf] = is->intersectionGlobal().volume(); // volume of face
								}

								FieldVector<ct,dim> unitOuterNormal
								= is->unitOuterNormal(facelocal); // normal vector of unit length

								const FieldVector<ct,dim>&
								facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
								//get boundary condition for boundary face center
								BoundaryConditions::Flags bctype = this->diffproblem.bctypePress(faceglobal, *it, facelocalDim);
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
									double lambda = lambdaI;

									double g = this->diffproblem.gPress(faceglobal, *it, facelocalDim);

									FieldVector<ct,dim> vTotal(Kni);
									vTotal *= lambda*(g-pressi)/dist;
									fineVelocity[numberInSelf] = vTotal;
								}
								else
								{
									double J = this->diffproblem.JPress(faceglobal, *it, facelocalDim);
									FieldVector<ct,dimworld> unitOuterNormal
									= is->unitOuterNormal(facelocal);
									fineVelocity[numberInSelf] = unitOuterNormal;
									fineVelocity[numberInSelf] *= J;
								}
								fineOnCoarse[coarseNumberInSelf][hits] = fineVelocity[numberInSelf];
								hitcount[coarseNumberInSelf] += 1;
							}
						}
					}

				}
				// end intersection traversal

				//				// check for conservativity
				//				if (dim == 2)
				//				{
				//					double diff = fabs(fineVelocity[0][0]*faceVol[0]
				//							- fineVelocity[1][0]*faceVol[1]
				//							+ fineVelocity[2][1]*faceVol[2]
				//							- fineVelocity[3][1]*faceVol[3])
				//					/(fabs(fineVelocity[0][0]*faceVol[0])
				//							+ fabs(fineVelocity[1][0]*faceVol[1])
				//							+ fabs(fineVelocity[2][1]*faceVol[2])
				//							+ fabs(fineVelocity[3][1]*faceVol[3]));
				//					if (diff> 1e-6)
				//					{
				//						std::cout << "NOT conservative!!! diff = " << diff << ", indexi = " << indexi << std::endl;
				//						std::cout << fineVelocity[0][0]*faceVol[0] << ", "
				//						<< fineVelocity[1][0]*faceVol[1] << ", "
				//						<< fineVelocity[2][1]*faceVol[2] << ", "
				//						<< fineVelocity[3][1]*faceVol[3] << std::endl;
				//					}
				//				}
			}// end hierarchic iteration
			// evaluate mean velocities and standard deviations at coarse element edges
			// and write them to the velocity struct.
			for (int i = 0; i<blocksize; i++)
			{
				FieldVector<ct,dim> mean(0);
				// the mean
				for (int j = 0; j<nFine; j++)
				{
					mean += fineOnCoarse[i][j];
				}
				mean/= nFine;

				// write to velocity struct
				this->diffproblem.variables.velocity[coarseIndex][i] = mean;
			}

		} // end grid traversal

		return;
	}// end method totalVelocity

};
}
#endif
