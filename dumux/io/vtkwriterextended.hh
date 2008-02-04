#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune { 


	template<class GridImp, class IS=typename GridImp::template Codim<0>::LeafIndexSet>
	class VTKWriterExtended : VTKWriter<GridImp, IS> {
	public: 
		typedef VTKWriter<GridImp, IS> VTKWriter;
		typedef typename VTKWriter::VTKFunction VTKFunction;
		
	    void addFaceData (VTKFunction* p)
	      {
	        //celldata.push_back(p);
	      }

	    template<class V>
	    void addFaceData (const V& v, std::string name)
	      {
	        //VTKFunction* p = new P0VectorWrapper<V>(grid,is,v,name);
	        //celldata.push_back(p);
	      }
	    
	    void addCellData (VTKFunction* p)
	      {
	        this->addCellData(p);
	      }

	    template<class V>
	    void addCellData (const V& v, std::string name)
	      {
	    	this->addCellData(v, name);
	      }

	    void write (const char* name, VTKOptions::OutputType ot = VTKOptions::ascii)
	      {
	    	this->write(name, ot);
	      }
	    
	    VTKWriterExtended (const GridImp& g, VTKOptions::DataMode dm = VTKOptions::conforming) :
	      VTKWriter(g, dm)
	      { }

	};
	
}


