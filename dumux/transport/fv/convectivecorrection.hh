// $Id:$
#ifndef DUNE_CONVECTIVECORRECTION_HH
#define DUNE_CONVECTIVECORRECTION_HH

#include "dumux/transport/fv/convectivepart.hh"
#include "dumux/fractionalflow/fractionalflowproblem.hh"
//#include "dumux/material/property_baseclasses.hh"

namespace Dune
{

template<class G, class RT, class VC>
class ConvectiveCorrection: public ConvectivePart<G, RT>
{
	enum
	{
		dim = G::dimension
	};
typedef	typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename Entity::Geometry Geometry;
	typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
	typedef FieldVector<RT, dim> VectorType;
	typedef typename G::ctype ct;

public:
	virtual double operator() (const Entity& entity, const RT satI, VectorType faceGlobal) const
	{
		double result(0);

		int size = problem.soil.getDispersionSat().M();

		int index = problem.variables.transmapper.map(entity);

		// element geometry
		const Geometry& geometry = entity.geometry();

		GeometryType gt = entity.type();
		// cell center in reference element
		const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0, 0);

		// get global coordinate of cell center
		const FieldVector<ct,dim> global = geometry.global(local);

		FieldVector<double,dim*2> fluxVector(0);

		VectorType unitOuterNormal(0);

		IntersectionIterator endis = entity.ilevelend();
		IntersectionIterator is = entity.ilevelbegin();
		for (; is != endis; ++is)
		{
			// get geometry type of face
			GeometryType gtf = is->intersectionSelfLocal().type();

			// center in face's reference element
			const FieldVector<ct,dim-1>& faceLocal = ReferenceElements<RT,dim-1>::general(gtf).position(0,0);

			// center of face in global coordinates
			Dune::FieldVector<ct,dim> faceGlobalCheck = is->intersectionGlobal().global(faceLocal);

			unitOuterNormal = is->unitOuterNormal(faceLocal);

			if (faceGlobal == faceGlobalCheck )
			{
				int faceNumber = is->numberInSelf();
				for (int i=0;i<size;i++)
				{
					if (satI == problem.soil.Sr_w(global, entity, local))
					{
						break;
					}
					if (problem.soil.getMSat()[index][i][faceNumber]>= satI)
					{
						//					std::cout<<"sat = "<<sat<<"satD = "<<problem.soil.getDispersionSatInterface()[indexI][i][numberInSelf]<<std::endl;
						double satdiff1 = problem.soil.getMSat()[index][i][faceNumber] - satI;
						double satdiff2 = 1e100;
						if (i)
						{
							satdiff2 = satI - problem.soil.getMSat()[index][i-1][faceNumber];
						}
						if (satdiff1 < satdiff2)
						{
							result += problem.soil.getM()[index][i][faceNumber];
//													std::cout<<"m = "<<problem.soil.getM()[index][i][faceNumber]<<std::endl;
//													std::cout<<"sat = "<<satI<<"satM = "<<problem.soil.getMSat()[index][i][faceNumber]<<std::endl;
						}
						else
						{
							result += problem.soil.getM()[index][i-1][faceNumber];
							//						std::cout<<"D = "<<problem.soil.getDispersionInterface()[indexI][i-1][numberInSelf]<<std::endl;
							//						std::cout<<"sat = "<<sat<<"satD = "<<problem.soil.getDispersionSatInterface()[indexI][i][numberInSelf]<<std::endl;
						}
						break;
					}
					if (i==(size-1))
					{
						result += problem.soil.getM()[index][i][faceNumber];
						break;
					}
				}
			}
		}
		double direction = 0;
		for (int i=0;i<dim;i++)
		{
			direction += unitOuterNormal[i];
		}

		return result;//*direction;

	}

	ConvectiveCorrection (FractionalFlowProblem<G, RT, VC>& prob)
	: problem(prob)

	{}

private:
	FractionalFlowProblem<G, RT, VC>& problem;
};
}

#endif
