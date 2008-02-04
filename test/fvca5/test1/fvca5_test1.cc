#include "config.h"
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fe/fediffusion.hh"
#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "fvca5test1problem.hh"
 
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
	  
template<class GridType, class ProblemType, class SolutionType>
double relativeL2Error(const GridType& grid, const ProblemType& problem, 
						const SolutionType& solution)
{
    typedef typename GridType::Traits::template Codim<0>::Entity Entity;
    typedef typename GridType::Traits::LevelIndexSet IS;
    typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
    typedef MultipleCodimMultipleGeomTypeMapper<GridType,IS,ElementLayout> EM;
    typedef typename GridType::ctype ct; 
    
    enum{dim = GridType::dimension};

    const IS& indexset(grid.levelIndexSet(grid.maxLevel()));
    EM elementmapper(grid, grid.levelIndexSet(grid.maxLevel()));
    
    double numerator = 0; 
    double denominator = 0;
    Iterator eendit = indexset.template end<0,All_Partition>();
    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
      {
    	// get entity 
    	const Entity& element = *it; 
    	
    	// cell geometry type
		GeometryType gt = element.geometry().type();
		
		// cell center in reference element
		const FieldVector<ct,dim>& 
		  local = ReferenceElements<ct,dim>::general(gt).position(0,0);
		
		// get global coordinate of cell center
		FieldVector<ct,dim> global = element.geometry().global(local);
		
		// get exact solution value 
		double exactValue = problem.exact(global);
		
		// cell index
		int indexi = elementmapper.map(element);
		
		// get approximate solution value 
		double approximateValue = solution[indexi];
		
		// cell volume, assume linear map here
		double volume = element.geometry().integrationElement(local)
			*ReferenceElements<ct,dim>::general(gt).volume();

		numerator += volume*(exactValue - approximateValue)*(exactValue - approximateValue);
		denominator += volume*exactValue*exactValue;
      }
    
    return sqrt(numerator/denominator);
}
}

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef double NumberType; 
    typedef Dune::UGGrid<dim> GridType; 

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( argv[1] );

    // grid reference 
    GridType& grid = *gridPtr;

    Dune::FVCA5Test1Problem<GridType, NumberType> problem;

    Dune::Timer timer;
    timer.reset();
    //Dune::FEDiffusion<GridType, NumberType> diffusion(grid, problem);
    //Dune::FVDiffusion<GridType, NumberType> diffusion(grid, problem);
    Dune::MimeticDiffusion<GridType, NumberType> diffusion(grid, problem);
    
    diffusion.pressure();
    std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
    //printvector(std::cout, *diffusion, "pressure", "row", 200, 1, 3);
    
    std::cout << "relative discrete L2 error erl2 = " << relativeL2Error(grid, problem, *diffusion) << std::endl;
    
    diffusion.vtkout("fvca5_test1", 0);


    
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
