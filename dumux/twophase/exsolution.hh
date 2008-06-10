#ifndef DUNE_EXSOLUTION_HH
#define DUNE_EXSOLUTION_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/utility/intersectiongetter.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/linearlaw.hh>

/**
 * @file
 * @brief  Basic class for embedding of an exact solution
 * @author Markus Wolff
 */

namespace Dune {

template<class G, class RT> class ExSolution {
	typedef typename G::ctype DT;
	enum {n=G::dimension,m=2};
	typedef BlockVector<FieldVector<RT, 1> > BV;
	typedef BlockVector<FieldVector<RT, m> > BVu;	
	typedef BlockVector<FieldVector<RT, (m+1)> > BVuEx;
	typedef typename G::Traits::LeafIndexSet IS;
	typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;	
	typedef typename G::Traits::template Codim<0>::Entity Entity;

	//Layout-helper-class
	template<int dim> struct P1Layout {
		bool contains(Dune::GeometryType gt) {
			if (gt.dim() == 0)
				return true;
			return false;
		}
	};

	typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P1Layout>
			VertexMapper;
	
	
	//private accsess functions
	void uExInit() {
		uEx.resize(mapper.size());
		uEx=0;
		return;
	}
protected:
	void uExIn(BV &uExIn) {
		uEx=uExIn;
		return;
	}

	void uExInVertex(RT &value, int &ElementIndex, int VariableIndex) {
		uEx[ElementIndex][VariableIndex]=value;
		return;
	}

public:

	ExSolution(const G &g) :
		grid(g), mapper(grid),uEx(0),error(0),elementvolume(0) {
			uExInit();
	}

	virtual ~ExSolution() {
	}

	//public access function
	BV& uExOut() const{
		return uEx;
	}
	
	RT uExOutVertex(int &ElementIndex, int VariableIndex) const{
		return uEx[ElementIndex][VariableIndex];
	}
	
	void calcSatError (BVu &Approx)
	    {
		int size=mapper.size();
		error.resize(size);
		elementvolume.resize(size);
		
		error=0;
		elementvolume=0;
		
		const IS& indexset(grid.leafIndexSet());
		typedef MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout> VMapper;
		VMapper vertexmapper(grid,grid.leafIndexSet());
		
		 Iterator eendit = indexset.template end<0, All_Partition>();
		      for (Iterator it = indexset.template begin<0, All_Partition>(); it != eendit; ++it)
			{
			  // get entity 
			  const Entity& entity = *it;

			  // get geometry type
			  Dune::GeometryType gt = it->geometry().type();
			  
			  double cellvolume = entity.geometry().volume();

			   // cell center in reference element
			  const Dune::FieldVector<DT,n> 
			    cellLocal = Dune::ReferenceElements<DT,n>::general(gt).position(0,0);
				
			  // get global coordinate of cell center
			  const Dune::FieldVector<DT,n> cellGlobal = it->geometry().global(cellLocal);

			  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
			    sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt, n);
			  int size = sfs.size();

			  // begin loop over vertices
			  for(int i=0; i < size; i++)
			    {
			      // local coordinate of vertex 
			      const FieldVector<DT,n> vertexLocal = sfs[i].position();
			      
			      // get global coordinate of vertex
			      const FieldVector<DT,n> vertexGlobal = it->geometry().global(vertexLocal);
			      
			      int globalId = vertexmapper.template map<n>(entity, sfs[i].entity());
			      
			      elementvolume[globalId]+=cellvolume/size;
//			      std::cout<<"elementvolume = "<<elementvolume[globalId]<<std::endl;
			    }
			}

		double globalvolume = elementvolume.one_norm();
//		std::cout<<"globalvolume = "<<globalvolume<<std::endl;
			
		for (int i=0;i<size;i++){
			error[i]=uEx[i][1]-Approx[i][1];
		}
			
	    double diffNorm = error.two_norm();
	    std::cout<<"diffNorm = "<<diffNorm<<std::endl;
	    
	    for (int i=0;i<size;i++)
	    	uEx[i][2] = diffNorm * pow((elementvolume[i]/globalvolume),0.5);
	    
	      return;
	    }

protected:
	const G &grid;
	VertexMapper mapper;
	BVuEx uEx;
	BV error;
	BV elementvolume;
};
}
#endif
