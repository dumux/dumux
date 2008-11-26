// $Id$
#ifndef DUNE_DISPERSIVECORRECTION_HH
#define DUNE_DISPERSIVECORRECTION_HH

#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/fractionalflow/fractionalflowproblem.hh"
//#include "dumux/material/property_baseclasses.hh"

namespace Dune
{

template<class G, class RT, class VC>
class DispersiveCorrection: public DiffusivePart<G, RT>
{
	enum
	{
		dim = G::dimension
	};
typedef	typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
	typedef FieldVector<RT, dim> VectorType;
	typedef typename G::ctype ct;

public:
	virtual VectorType operator() (const Entity& entity, const int numberInSelf,
			const RT satIntersection, const VectorType& satGradient, const RT time,
			const RT satI, const RT satJ) const
	{
		VectorType result(0);
		VectorType helpresult(0);

		int size = problem.soil.getDispersionSat().M();

		int indexI = problem.variables.transmapper.map(entity);

		GeometryType gt = entity.type();
		// cell center in reference element
		const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0, 0);

		// get global coordinate of cell center
		const FieldVector<ct,dim> global = entity.geometry().global(local);

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
		const FieldVector<ct,dim-1>& facelocal = ReferenceElements<RT,dim-1>::general(gtf).position(0,0);

		VectorType unitOuterNormal = is->unitOuterNormal(facelocal);
		//		std::cout<<"unitOuterNormaldiff"<<unitOuterNormal<<std::endl;

		if (is->neighbor())
		{
			// access neighbor
			EntityPointer outside = is->outside();

			int indexJ = problem.variables.transmapper.map(*outside);

			for (int i=0;i<size;i++)
			{
				if (problem.soil.getDispersionSat()[indexI][i]> satI)
				{
					double satdiff1 = problem.soil.getDispersionSat()[indexI][i] - satI;
					double satdiff2 = 1e100;
					if (i)
					{
						satdiff2 = satI - problem.soil.getDispersionSat()[indexI][i-1];
					}
					if (satdiff1 < satdiff2)
					{
						problem.soil.getDispersion()[indexI][i].mv(satGradient,result);
						//						std::cout<<"result1 = "<<problem.soil.getDispersion()[indexI][i]<<std::endl;
					}
					else
					{
						problem.soil.getDispersion()[indexI][i-1].mv(satGradient,result);
						//						std::cout<<"result2 = "<<problem.soil.getDispersion()[indexI][i-1]<<std::endl;
						//						std::cout<<"satGrad = "<<satGradient<<std::endl;
					}
					break;
				}
				if (i==(size-1))
				{
					problem.soil.getDispersion()[indexI][i].mv(satGradient,helpresult);
					break;
				}
			}
			for (int i=0;i<size;i++)
			{
				if (problem.soil.getDispersionSat()[indexJ][i]> satJ)
				{
					double satdiff1 = problem.soil.getDispersionSat()[indexJ][i] - satJ;
					double satdiff2 = 1e100;
					if (i)
					{
						satdiff2 = satJ - problem.soil.getDispersionSat()[indexJ][i-1];
					}
					if (satdiff1 < satdiff2)
					{
						problem.soil.getDispersion()[indexJ][i].mv(satGradient,helpresult);
						//						std::cout<<"helpresult1 = "<<problem.soil.getDispersion()[indexJ][i]<<std::endl;
					}
					else
					{
						//						std::cout<<"D = "<<problem.soil.getDispersion()[indexJ][i-1]<<std::endl;
						//						std::cout<<"satG = "<<satGradient<<std::endl;
						problem.soil.getDispersion()[indexJ][i-1].mv(satGradient,helpresult);
						//						std::cout<<"helpresult2 = "<<problem.soil.getDispersion()[indexJ][i-1]<<std::endl;
					}
					break;
				}
				if (i==(size-1))
				{
					problem.soil.getDispersion()[indexJ][i].mv(satGradient,helpresult);
					break;
				}
			}

			for (int i = 0; i< dim;i++)
			{
				if (result[i]*helpresult[i]==0)
					result[i]=0;
				else
					result[i] = 2*result[i]*helpresult[i]/(result[i]+helpresult[i]);
			}
			//			for (int i=0;i<dim;i++)
			//			{
			//				result[i]*=unitOuterNormal[i];
			//			}

		}
		else
		{
			for (int i=0;i<size;i++)
			{
				if (problem.soil.getDispersionSat()[indexI][i]> satJ)
				{
					double satdiff1 = problem.soil.getDispersionSat()[indexI][i] - satJ;
					double satdiff2 = 1e100;
					if (i)
					{
						satdiff2 = satJ - problem.soil.getDispersionSat()[indexI][i-1];
					}
					if (satdiff1 < satdiff2)
					{
						problem.soil.getDispersion()[indexI][i].mv(satGradient,result);
						//						std::cout<<"result1Bound = "<<problem.soil.getDispersion()[indexI][i]<<std::endl;
					}
					else
					{
						problem.soil.getDispersion()[indexI][i-1].mv(satGradient,result);
						//						std::cout<<"result2Bound = "<<problem.soil.getDispersion()[indexI][i]<<std::endl;
					}
					break;
				}
				if (i==(size-1))
				{
					problem.soil.getDispersion()[indexI][i].mv(satGradient,result);
					break;
				}
			}
			result *= 1;

			//			for (int i=0;i<dim;i++)
			//			{
			//				result[i]*=unitOuterNormal[i];
			//			}
		}

		return result;

	}

	DispersiveCorrection (FractionalFlowProblem<G, RT, VC>& prob)
	: problem(prob)

	{}

private:
	FractionalFlowProblem<G, RT, VC>& problem;
};
}

#endif
