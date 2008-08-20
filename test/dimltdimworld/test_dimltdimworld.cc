#include "config.h"
#include <iostream>
#include <iomanip>
#include "dumux/onedinndgrid/onedinndgrid.hh"
#include <dune/grid/common/gridinfo.hh> 
 
int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=1;
    const int dimworld=3;

    // create a grid object
    typedef double NumberType; 
    typedef Dune::OneDInNDGrid<dimworld> GridType; 
    typedef GridType::ctype CType;

    Dune::FieldVector<CType, dimworld> ll(1.0);
    ll[1] = 1.0; ll[dimworld-1] = 2.0;
    Dune::FieldVector<CType, dimworld> ur(2.0);
    ur[1] = 3.0; ur[dimworld-1] = 5.0;
    GridType grid(2, ll, ur);
 
    grid.globalRefine(1); 
    
    typedef GridType::Codim<0>::LeafIterator Iterator;
    Iterator ebeginit = grid.leafbegin<0>();
    Iterator eendit = grid.leafend<0>();
    Iterator lastit(eendit);
    for (Iterator it = ebeginit; it != eendit; ++it) 
      lastit = it;
    

    Dune::GeometryType gtype = ebeginit->geometry().type();
    const Dune::FieldVector<CType,dim>& beginlocal = Dune::ReferenceElements<CType,dim>::general(gtype).position(0,1);
    const Dune::FieldVector<CType,dim>& endlocal = Dune::ReferenceElements<CType,dim>::general(gtype).position(1,1);
    Dune::FieldVector<CType,dimworld> lowerLeft = ebeginit->geometry().global(beginlocal);
    Dune::FieldVector<CType,dimworld> upperRight = lastit->geometry().global(endlocal); 

    Dune::FieldVector<CType,dimworld> unitDirection = upperRight - lowerLeft;
    unitDirection /= unitDirection.two_norm();
    std::cout << "lowerLeft = " << lowerLeft << ", upperRight = " << upperRight 
	      << ", unitDirection = " << unitDirection << std::endl;

    for (Iterator it = ebeginit; it != eendit; ++it) {
	// cell geometry type
	Dune::GeometryType gt = it->geometry().type();
	
	// cell center in reference element
	const Dune::FieldVector<CType,dim>& 
	  local = Dune::ReferenceElements<CType,dim>::general(gt).position(0,0);
	
	// get global coordinate of cell center
	Dune::FieldVector<CType,dimworld> global = it->geometry().global(local);

	std::cout << "element center: local = " << local << ", global = " << global << std::endl;

	// intersection iterator type
	typedef GridType::Codim<0>::LeafIntersectionIterator IntersectionIterator;

	IntersectionIterator isend = it->ileafend(); 
	for (IntersectionIterator is = it->ileafbegin(); is!=isend; ++is)
	  {
	    const Dune::FieldVector<CType,dim>& 
	      localInRef = Dune::ReferenceElements<CType,dim>::general(gt).position(is->numberInSelf(),1);

	    // get global coordinate of face
	    Dune::FieldVector<CType,dimworld> faceglobal = it->geometry().global(localInRef);

	    const Dune::FieldVector<CType,dim-1> facelocal(0);

	    Dune::FieldVector<CType,dimworld> unitOuterNormal = is->unitOuterNormal(facelocal);

	    std::cout << "local face id: " << is->numberInSelf() << ", faceglobal = " << faceglobal 
		      << ", normal " << unitOuterNormal << std::endl;
	  }
    }	
    
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
