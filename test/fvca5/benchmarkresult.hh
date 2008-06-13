#ifndef DUNE_BENCHMARKRESULT_HH
#define DUNE_BENCHMARKRESULT_HH

namespace Dune
{

template<int dim>
struct ElementLayout
{
	bool contains (GeometryType gt)
	{
		return gt.dim() == dim;
	}
}; 

template<int dim>
struct FaceLayout
{
	bool contains (GeometryType gt)
	{
		return gt.dim() == dim-1;
	}
}; 

struct BenchmarkResult 
{
	double relativeL2Error;
	double ergrad;
	double ervell2;
	double uMin; 
	double uMax;
	double flux0;
	double flux1;
	double fluy0;
	double fluy1;
	double sumf;
	double sumflux;
	double exactflux0;
	double exactflux1;
	double exactfluy0;
	double exactfluy1;
	double errflx0;
	double errflx1;
	double errfly0;
	double errfly1;
	double erflm;
	double ener1;
	double ener2;
	double eren;
	double uMean;
	
	template<class GridType, class ProblemType, class SolutionType>
	void evaluate(const GridType& grid, ProblemType& problem, 
							SolutionType& solution, bool pureNeumann = false)
	{
	    typedef typename GridType::Traits::template Codim<0>::Entity Entity;
	    typedef typename Entity::Geometry Geometry;
	    typedef typename GridType::Traits::LevelIndexSet IS;
	    typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
	    typedef typename IntersectionIteratorGetter<GridType,LeafTag>::IntersectionIterator IntersectionIterator;
	    typedef MultipleCodimMultipleGeomTypeMapper<GridType,IS,ElementLayout> EM;
		typedef MultipleCodimMultipleGeomTypeMapper<GridType,IS,FaceLayout> FM;
	    typedef typename GridType::ctype ct; 
	    
	    enum{dim = GridType::dimension};

	    const IS& indexset(grid.levelIndexSet(grid.maxLevel()));
	    EM elementmapper(grid, grid.levelIndexSet(grid.maxLevel()));
	    FM facemapper(grid, grid.levelIndexSet(grid.maxLevel()));
	    
	    uMean = 0;
	    double domainVolume = 0;
	    Iterator eendit = indexset.template end<0,All_Partition>();
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	      {
	    	// get entity 
	    	const Entity& element = *it; 
	    	
	    	// get volume 
	    	double volume = element.geometry().volume();
	    	
			// cell index
			int indexi = elementmapper.map(element);
			
			// get approximate solution value 
			uMean += volume*(*solution)[indexi];
			
			// add to domainVolume 
			domainVolume += volume;
	      }
	    uMean /= domainVolume; 
	    
	    if (pureNeumann) {
	    	for (int i = 0; i < (*solution).size(); i++)
	    	(*solution)[i] -= uMean;
	    }
	    
	    
	    uMin = 1e100;
	    uMax = -1e100;
	    flux0 = 0;
	    flux1 = 0;
	    fluy0 = 0;
	    fluy1 = 0;
	    sumf = 0;
	    sumflux = 0;
	    exactflux0 = 0;
	    exactflux1 = 0;
	    exactfluy0 = 0;
	    exactfluy1 = 0;
	    erflm = 0;
	    ener1 = 0;
	    ener2 = 0;
	    double numerator = 0; 
	    double denominator = 0;
	    double numeratorGrad = 0; 
	    double denominatorGrad = 0;
	    double numeratorFlux = 0; 
	    double denominatorFlux = 0;
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	      {
	    	// get entity 
	    	const Entity& element = *it; 
	    	
	    	// element geometry
	    	const Geometry& geometry = element.geometry();
	    	
	    	// cell geometry type
			GeometryType gt = geometry.type();
			
			const typename CRShapeFunctionSetContainer<double,double,dim>::value_type& 
				sfs=CRShapeFunctions<double,double,dim>::general(gt,1);
				
			// cell center in reference element
			const FieldVector<ct,dim>& 
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0);
			
			// get global coordinate of cell center
			FieldVector<ct,dim> global = geometry.global(local);
			
			// get exact solution value 
			double exactValue = problem.exact(global);
			
			// cell index
			int indexi = elementmapper.map(element);
			
			// get approximate solution value 
			double approximateValue = (*solution)[indexi];
			
			// update uMin and uMax
			uMin = std::min(uMin, approximateValue);
			uMax = std::max(uMax, approximateValue);
			
			// cell volume, assume linear map here
			double volume = geometry.volume();

			// update sumf 
			sumf += volume*problem.q(global, element, local);
				
			// get the absolute permeability
			FieldMatrix<double,dim,dim> K = problem.K(global, element, local);

			numerator += volume*(exactValue - approximateValue)*(exactValue - approximateValue);
			denominator += volume*exactValue*exactValue;

			FieldVector<ct,2*dim> fluxVector;
			FieldVector<ct,dim> exactGradient; 
			IntersectionIterator endis = IntersectionIteratorGetter<GridType,LeafTag>::end(element);
			for (IntersectionIterator is = IntersectionIteratorGetter<GridType,LeafTag>::begin(element); is!=endis; ++is)
			{
				// get geometry type of face
				GeometryType gtf = is->intersectionSelfLocal().type();
			  
				// local number of facet 
				int i = is->numberInSelf();

				// global number of face 
			    int faceIndex = facemapper.template map<1>(element, i);

			    const FieldVector<double,dim>& faceLocal = sfs[i].position();
				FieldVector<double,dim> faceGlobal = geometry.global(faceLocal);
				double faceVol = is->intersectionGlobal().volume();
			  
				// center in face's reference element
				const FieldVector<double,dim-1>& faceLocalNm1 = ReferenceElements<double,dim-1>::general(gtf).position(0,0);
				
				// get normal vector 
				FieldVector<double,dim> unitOuterNormal = is->unitOuterNormal(faceLocalNm1);
				
				// get the approximate solution on the face
				double approximateFace = (*solution.pressTrace)[faceIndex];
				
				// get the exact gradient
				exactGradient = problem.exactGrad(faceGlobal);
				
				// get the negative exact velocity 
				FieldVector<double,dim> KGrad(0);
				K.umv(exactGradient, KGrad);
				
				// calculate the exact normal velocity
				double exactFlux = KGrad*unitOuterNormal;
				
				// get the approximate normalvelocity
				double approximateFlux = (*solution.normalVelocity)[indexi][i];
				
				// calculate the difference in the normal velocity 
				double fluxDiff = exactFlux + approximateFlux;
				
				// update mean value error 
				erflm = std::max(erflm, fabs(fluxDiff));
				
				numeratorFlux += volume*fluxDiff*fluxDiff;
				denominatorFlux += volume*exactFlux*exactFlux;
				
				// calculate the fluxes through the element faces 
				exactFlux *= faceVol;
				approximateFlux *= faceVol;
				fluxVector[i] = approximateFlux;
				
				//if (is.boundary()) {
				if (!is->neighbor()) {
					if (fabs(faceGlobal[1]) < 1e-6) {
						fluy0 += approximateFlux;
						exactfluy0 += exactFlux;
						ener2 += -approximateFlux*approximateFace;
					}
					else if (fabs(faceGlobal[1] - 1.0) < 1e-6) {
						fluy1 += approximateFlux;
						exactfluy1 += exactFlux;
						ener2 += -approximateFlux*approximateFace;
					}
					else if (faceGlobal[0] < 1e-6) {
						flux0 += approximateFlux;
						exactflux0 += exactFlux;
						ener2 += -approximateFlux*approximateFace;
					}
					else if (fabs(faceGlobal[0] - 1.0) < 1e-6) {
						flux1 += approximateFlux;
						exactflux1 += exactFlux;
						ener2 += -approximateFlux*approximateFace;
					}
				}
			}
			
			// calculate velocity on reference element 
			FieldVector<ct,dim> refVelocity; 
			if (geometry.corners() == 3) {
				refVelocity[0] = 1.0/3.0*(fluxVector[0] + fluxVector[2] - 2.0*fluxVector[1]);
				refVelocity[1] = 1.0/3.0*(fluxVector[0] + fluxVector[1] - 2.0*fluxVector[2]);
			}
			else {
				refVelocity[0] = 0.5*(fluxVector[1] - fluxVector[0]);
				refVelocity[1] = 0.5*(fluxVector[3] - fluxVector[2]);
			}
			
			// get the transposed Jacobian of the element mapping 
			const FieldMatrix<ct,dim,dim>& jacobianInv = geometry.jacobianInverseTransposed(local);
			FieldMatrix<ct,dim,dim> jacobianT(jacobianInv);
			jacobianT.invert();
			
			// calculate the element velocity by the Piola transformation
			FieldVector<ct,dim> elementVelocity(0);
			jacobianT.umtv(refVelocity, elementVelocity);
			elementVelocity /= geometry.integrationElement(local);
			
			// get the approximate gradient 
			FieldVector<ct,dim> approximateGradient;
			K.solve(approximateGradient, elementVelocity);
			
			// get the exact gradient
			exactGradient = problem.exactGrad(global);
			
			// the difference between exact and approximate gradient 
			FieldVector<ct,dim> gradDiff(exactGradient);
			gradDiff += approximateGradient; 

			// add to energy 
			ener1 += volume*(approximateGradient*elementVelocity);
			
			numeratorGrad += volume*(gradDiff*gradDiff);
			denominatorGrad += volume*(exactGradient*exactGradient);
	      }
	    
	    relativeL2Error = sqrt(numerator/denominator);
	    ergrad = sqrt(numeratorGrad/denominatorGrad);
	    ervell2 = sqrt(numeratorFlux/denominatorFlux);
	    sumflux = flux0 + flux1 + fluy0 + fluy1 - sumf;
	    errflx0 = fabs((flux0 + exactflux0)/exactflux0);
	    errflx1 = fabs((flux1 + exactflux1)/exactflux1);
	    errfly0 = fabs((fluy0 + exactfluy0)/exactfluy0);
	    errfly1 = fabs((fluy1 + exactfluy1)/exactfluy1);
	    eren = fabs(ener1 - ener2)/std::max(ener1, ener2);
	    
	    return;
	}
};

}

#endif
