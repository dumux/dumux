#ifndef DUNE_FVELEMENTGEOMETRY_HH
#define DUNE_FVELEMENTGEOMETRY_HH

namespace Dune
{

template<class G>
class FVElementGeometry  
{
public:
	enum{dim = G::dimension};
	enum{maxNC = (dim < 3 ? 4 : 8)};
	enum{maxNF = (dim < 3 ? 4 : 12)};
	typedef typename G::ctype DT;
   typedef typename G::Traits::template Codim<0>::Entity Entity;
   typedef typename Entity::Geometry Geometry;
    
	struct SDValues 
	{
		FieldVector<DT, maxNC> shape;                    /* values of shape functions at ip              */
		FieldVector<FieldVector<DT, dim>, maxNC> grad;              /* derivatives of shape functions at ip */
		FieldMatrix<DT, dim, dim> J;                         /* jacobian at ip                               */
		FieldMatrix<DT, dim, dim> Jinv;                  /* inverse of jacobian at ip                    */
		DT detJ;                                    /* det of jacobian at ip                                */
	};

	struct SubControlVolume /* FV intersected with element                  */
	{
		int co;                                                 /* # of corner                                                  */
		FieldVector<DT, dim> local;                   /* node position                                                */
		FieldVector<DT, dim> global;                   /* node position                                                */
	    DT volume;                                  /* volume (area) of scv                                 */
	};                                     

	struct SubControlVolumeFace
	{
		int i,j;                                                /* scvf seperates corner i and j of elem*/
		FieldVector<DT, dim> ip_local;                 /* integration point in local coords    */
		FieldVector<DT, dim> ip_global;                    /* integration point in global coords       */
		FieldVector<DT, dim> normal;                   /* normal on face at ip pointing to CV j*/
		SDValues sdv;                                  /* shape fcts, deriv. etc. at scv-faces */
	};

	struct BoundaryFace {
		int co;                                                 /* corresponding corner                                 */
		int side;                                               /* boundary side of element                             */
		FieldVector<DT, dim> ip_local;                 /* integration point in local coords    */
		FieldVector<DT, dim-1> param;            /* local side coordinates                       */
		FieldVector<DT, dim> normal;                   /* normal on face at ip pointing to CV j*/
		DT area;                                    /* area of boundary face                                */
		SDValues sdv;                                  /* shape fcts, deriv. etc. at b-faces   */
	}; 

   FieldVector<DT, dim> cellLocal;
   FieldVector<DT, dim> cellGlobal;
   DT cellVolume;
   SubControlVolume subContVol[maxNC];
   SubControlVolumeFace subContVolFace[maxNF];
    
	FVElementGeometry<G>(const Entity& e)
	{
      const Geometry& geometry = e.geometry();
      GeometryType gt = geometry.type();
       
       const typename LagrangeShapeFunctionSetContainer<DT,DT,dim>::value_type& 
       	sfs=LagrangeShapeFunctions<DT,DT,dim>::general(gt, 1);
       	
       // get cell volume
       cellVolume = geometry.volume();
       
       // cell center in reference element
 		 cellLocal = ReferenceElements<DT,dim>::general(gt).position(0,0);
 		
 	  
 		 // get global coordinate of cell center
 		 cellGlobal = geometry.global(cellLocal);
   	 
 		 int nodes = sfs.size();
 		 for (int i = 0; i < nodes; i++) {
 			 // assuming rectangular grids:
 			 subContVol[i].volume = 0.25*cellVolume;
 			 
 			 subContVol[i].local = sfs[i].position();
 			 
 			 subContVol[i].global = geometry.global(subContVol[i].local);
 		 }

 		 int faces = e.template count<dim-1>();
 		 for (int k = 0; k < faces; k++) {
 			 int idx0 = ReferenceElements<DT,dim>::general(gt).subEntity(k, dim-1, 0, dim);
 			 int idx1 = ReferenceElements<DT,dim>::general(gt).subEntity(k, dim-1, 1, dim);
 			 subContVolFace[k].i = std::min(idx0, idx1);
 			 subContVolFace[k].j = std::max(idx0, idx1); 
 		 }
 		 
	}

};

}


#endif



