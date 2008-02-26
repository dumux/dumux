#ifndef DUNE_FVELEMENTGEOMETRY_HH
#define DUNE_FVELEMENTGEOMETRY_HH

namespace Dune
{

template<class G>
class FVElementGeometry  
{
	enum{dim = G::dimension};
	enum{maxNC = (dim < 3 ? 4 : 8)};
	enum{maxNE = (dim < 3 ? 4 : 12)};
	enum{maxNF = (dim < 3 ? 0 : 6)};
	enum{maxCOS = (dim < 3 ? 2 : 4)};
	enum{maxBF = (dim < 3 ? 8 : 24)};
	typedef typename G::ctype DT;
   typedef typename G::Traits::template Codim<0>::Entity Entity;
   typedef typename Entity::Geometry Geometry;
   typedef FieldVector<DT,dim> FV;
   typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;
   
   DT quadrilateralArea(const FV& p0, const FV& p1, const FV& p2, const FV& p3)
   {
	   return 0.5*fabs((p3[0] - p1[0])*(p2[1] - p0[1]) - (p3[1] - p1[1])*(p2[0] - p0[0]));
   }

   FV crossProduct(const FV& a, const FV& b)
   {
	   FV c;
	   c[0] = a[1]*b[2] - a[2]*b[1];
	   c[1] = a[2]*b[0] - a[0]*b[2];
	   c[2] = a[0]*b[1] - a[1]*b[0];
	   
	   return c;
   }
   
   DT pyramidVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3, const FV& p4)
   {
	   FV a = p2 - p0;
	   FV b = p3 - p1;
	   FV h = p4 - p0;
	   FV n = crossProduct(a, b);
   	
	   return(1.0/6.0*(n*h));
   }
   
   DT prismVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3, const FV& p4, const FV& p5)
   {
	   FV a = p4 - p0;
	   FV b = p1 - p3;
	   FV c = p1 - p0;
	   FV d = p2 - p0;
	   FV e = p5 - p0;
	   FV m = crossProduct(a, b);
	   FV n = m + crossProduct(c, d);
   	
	   return(fabs(1.0/6.0*(n*e)));
   }

   DT hexahedronVolume (const FV& p0, const FV& p1, const FV& p2, const FV& p3, 
		   				const FV& p4, const FV& p5, const FV& p6, const FV& p7)
   {
	   return(prismVolume(p0,p1,p2,p4,p5,p6)
			   + prismVolume(p0,p2,p3,p4,p6,p7));
   }
	  
   FV normalOfQuadrilateral3D(const FV& p0, const FV& p1, const FV& p2, const FV& p3)   
   {
	   FV a = p2 - p0;
	   FV b = p3 - p1;
	   FV normal = crossProduct(a, b);
	   normal *= 0.5; 
	   
	   return normal;
   }

   DT quadrilateralArea3D(const FV& p0, const FV& p1, const FV& p2, const FV& p3)   
   {
	   return (normalOfQuadrilateral3D(p0, p1, p2, p3).two_norm());
   }
   
   void getFaceIndices(int nNodes, int k, int& leftFace, int& rightFace)
   {
	   int edgeToFaceTet[2][6] = {
			   {2, 0, 3, 1, 2, 0},
	   		   {3, 3, 1, 2, 0, 1}
	   };
	   int edgeToFacePyramid[2][8] = {
			   {1, 2, 3, 4, 4, 1, 2, 3},
	   		   {0, 0, 0, 0, 1, 2, 3, 4}
	   };
	   int edgeToFacePrism[2][9] = {
			   {1, 2, 3, 3, 1, 2, 4, 4, 4},
	   		   {0, 0, 0, 1, 2, 3, 1, 2, 3}
	   };
	   int edgeToFaceHex[2][12] = {
			   {0, 2, 3, 1, 4, 1, 0, 5, 2, 4, 5, 3},
	   		   {2, 1, 0, 3, 0, 4, 5, 1, 4, 3, 2, 5}
	   };

	  switch (nNodes) {
	  case 4:
		  leftFace = edgeToFaceTet[0][k];
		  rightFace = edgeToFaceTet[1][k];
		  break;
	  case 5:
		  leftFace = edgeToFacePyramid[0][k];
		  rightFace = edgeToFacePyramid[1][k];
		  break;
	  case 6:
		  leftFace = edgeToFacePrism[0][k];
		  rightFace = edgeToFacePrism[1][k];
		  break;
	  case 8:
		  leftFace = edgeToFaceHex[0][k];
		  rightFace = edgeToFaceHex[1][k];
		  break;
	  default:
		  DUNE_THROW(NotImplemented, "FVElementGeometry :: getFaceIndices for nNodes = " << nNodes); 
		  break;
	  }
	  return;
   }
   
   void getEdgeIndices(int nNodes, int face, int node, int& leftEdge, int& rightEdge)
   {
	   int faceAndNodeToLeftEdgeTet[4][4] = {
			   {-1,  1,  1,  4},
			   { 2, -1,  2,  3},
			   { 0,  0, -1,  3},
			   { 0,  0,  1, -1}
	   };
	   int faceAndNodeToRightEdgeTet[4][4] = {
			   {-1,  4,  5,  5},
			   { 3, -1,  5,  5},
			   { 3,  4, -1,  4},
			   { 2,  1,  2, -1}
	   };
	   int faceAndNodeToLeftEdgePyramid[5][5] = {
			   { 3,  0,  1,  2, -1},
			   { 0,  0, -1, -1,  4},
			   {-1,  1,  1, -1,  5},
			   {-1, -1,  2,  2,  6},
			   { 3, -1, -1,  3,  4}
	   };
	   int faceAndNodeToRightEdgePyramid[5][5] = {
			   { 0,  1,  2,  3, -1},
			   { 4,  5, -1, -1,  5},
			   {-1,  5,  6, -1,  6},
			   {-1, -1,  6,  7,  7},
			   { 4, -1, -1,  7,  7}
	   };
	   int faceAndNodeToLeftEdgePrism[5][6] = {
			   { 0,  0,  1, -1, -1, -1},
			   { 0,  0, -1,  3,  4, -1},
			   {-1,  1,  1, -1,  4,  5},
			   { 2, -1,  2,  3, -1,  5},
			   {-1, -1, -1,  6,  6,  7}
	   };
	   int faceAndNodeToRightEdgePrism[5][6] = {
			   { 2,  1,  2, -1, -1, -1},
			   { 3,  4, -1,  6,  6, -1},
			   {-1,  4,  5, -1,  7,  7},
			   { 3, -1,  5,  8, -1,  8},
			   {-1, -1, -1,  8,  7,  8}
	   };
	   int faceAndNodeToLeftEdgeHex[6][8] = {
			   { 0, -1,  4, -1,  6, -1,  2, -1},
			   {-1,  5, -1,  3, -1,  1, -1,  7},
			   { 8,  1, -1, -1,  0, 10, -1, -1},
			   {-1, -1,  2,  9, -1, -1, 11,  3},
			   { 4,  8,  9,  5, -1, -1, -1, -1},
			   {-1, -1, -1, -1, 10,  7,  6, 11}
	   };
	   int faceAndNodeToRightEdgeHex[6][8] = {
			   { 4, -1,  2, -1,  0, -1,  6, -1},
			   {-1,  1, -1,  5, -1,  7, -1,  3},
			   { 0,  8, -1, -1, 10,  1, -1, -1},
			   {-1, -1,  9,  3, -1, -1,  2, 11},
			   { 8,  5,  4,  9, -1, -1, -1, -1},
			   {-1, -1, -1, -1,  6, 10, 11,  7}
	   };

	  switch (nNodes) {
	  case 4:
		  leftEdge = faceAndNodeToLeftEdgeTet[face][node];
		  rightEdge = faceAndNodeToRightEdgeTet[face][node];
		  break;
	  case 5:
		  leftEdge = faceAndNodeToLeftEdgePyramid[face][node];
		  rightEdge = faceAndNodeToRightEdgePyramid[face][node];
		  break;
	  case 6:
		  leftEdge = faceAndNodeToLeftEdgePrism[face][node];
		  rightEdge = faceAndNodeToRightEdgePrism[face][node];
		  break;
	  case 8:
		  leftEdge = faceAndNodeToLeftEdgeHex[face][node];
		  rightEdge = faceAndNodeToRightEdgeHex[face][node];
		  break;
	  default:
		  DUNE_THROW(NotImplemented, "FVElementGeometry :: getFaceIndices for nNodes = " << nNodes); 
		  break;
	  }
	  return;
   }
public:
	int boundaryFaceIndex(int face, int nodeInFace) const
	{
		return (face*maxCOS + nodeInFace);
	}

	struct SubControlVolume //!< FV intersected with element                  
	{
		FV local; //!< local node position                                                
		FV global; //!< global node position                                                
	    DT volume; //!< volume of scv                                 
	};                                     

	struct SubControlVolumeFace
	{
		int i,j; //!< scvf seperates corner i and j of elem
		FV ipLocal; //!< integration point in local coords    
		FV ipGlobal; //!< integration point in global coords       
		FV normal; //!< normal on face at ip pointing to CV j with length equal to |scvf| 
		FieldVector<FV, maxNC> grad; //!< derivatives of shape functions at ip 
	};

	struct BoundaryFace {
		FV ipLocal; //!< integration point in local coords    
		FV ipGlobal; //!< integration point in global coords    
		DT area; //!< area of boundary face                                
	}; 

   FV cellLocal; //!< local coordinate of cell center 
   FV cellGlobal; //!< global coordinate of cell center 
   DT cellVolume; //!< cell volume 
   SubControlVolume subContVol[maxNC]; //!< data of the sub control volumes 
   SubControlVolumeFace subContVolFace[maxNE]; //!< data of the sub control volume faces 
   BoundaryFace boundaryFace[maxBF]; //!< data of the boundary faces 
   FV edgeCoord[maxNE]; //!< global coordinates of the edge centers 
   FV faceCoord[maxNF]; //!< global coordinates of the face centers
   int nNodes; //!< number of nodes 
   int nEdges; //!< number of edges 
   int nFaces; //!< number of faces (0 in < 3D)
    
   FVElementGeometry<G>()
   {}
   
   void update(const Entity& e)
	{
		const Geometry& geometry = e.geometry();
		GeometryType gt = geometry.type();
       
		const typename ReferenceElementContainer<DT,dim>::value_type& 
      	referenceElement = ReferenceElements<DT,dim>::general(gt);

      	const typename LagrangeShapeFunctionSetContainer<DT,DT,dim>::value_type& 
       	sfs=LagrangeShapeFunctions<DT,DT,dim>::general(gt, 1);
       	
       	cellVolume = geometry.volume();
       	cellLocal = referenceElement.position(0,0);
       	cellGlobal = geometry.global(cellLocal);
   	 
       	nNodes = referenceElement.size(dim);
       	nEdges = referenceElement.size(dim-1);
       	nFaces = dim < 3 ? 0 : referenceElement.size(1);
 		 
       	// corners:
       	for (int node = 0; node < nNodes; node++) { 
 			 subContVol[node].local = referenceElement.position(node, dim);
 			 subContVol[node].global = geometry.global(subContVol[node].local);
       	} 
 		 
       	// edges:
       	for (int edge = 0; edge < nEdges; edge++) { 
 			 edgeCoord[edge] = geometry.global(referenceElement.position(edge, dim-1));
       	} 
 		 
       	// faces:
       	for (int face = 0; face < nFaces; face++) { 
 			 faceCoord[face] = geometry.global(referenceElement.position(face, 1));
       	} 
 		 
       	// fill sub control volume data:
//       	switch (nEdges) {
//       	case 1: // 1D
//       		subContVol[0].volume = 0.5*cellVolume;
//       		break;
//       	case 3: // 2D, triangle 
//       		subContVol[0].volume = quadrilateralArea(subContVol[0].global, edgeCoord[2], cellGlobal, edgeCoord[1]);
//       		subContVol[1].volume = quadrilateralArea(subContVol[1].global, edgeCoord[0], cellGlobal, edgeCoord[2]);
//       		subContVol[2].volume = quadrilateralArea(subContVol[2].global, edgeCoord[1], cellGlobal, edgeCoord[0]);
//       		break;
//       	case 4: // 2D, quadrilateral
       		subContVol[0].volume = quadrilateralArea(subContVol[0].global, edgeCoord[2], cellGlobal, edgeCoord[0]);
       		subContVol[1].volume = quadrilateralArea(subContVol[1].global, edgeCoord[1], cellGlobal, edgeCoord[2]);
       		subContVol[2].volume = quadrilateralArea(subContVol[2].global, edgeCoord[0], cellGlobal, edgeCoord[3]);
       		subContVol[3].volume = quadrilateralArea(subContVol[3].global, edgeCoord[3], cellGlobal, edgeCoord[1]);
//       		break;
//       	case 6: // 3D, tetrahedron
//       		for (int k = 0; k < nNodes; k++)
//       			subContVol[k].volume = 0.25*cellVolume;
//       		break;
//       	case 8: // 3D, pyramid
//       		subContVol[0].volume = hexahedronVolume(subContVol[0].global, edgeCoord[0], faceCoord[0], edgeCoord[3], 
//       												edgeCoord[4], faceCoord[1], cellGlobal, faceCoord[4]);
//       		subContVol[1].volume = hexahedronVolume(subContVol[1].global, edgeCoord[1], faceCoord[0], edgeCoord[0], 
//       												edgeCoord[5], faceCoord[2], cellGlobal, faceCoord[1]);
//      		subContVol[2].volume = hexahedronVolume(subContVol[2].global, edgeCoord[2], faceCoord[0], edgeCoord[1], 
//       												edgeCoord[6], faceCoord[3], cellGlobal, faceCoord[2]);
//      		subContVol[3].volume = hexahedronVolume(subContVol[3].global, edgeCoord[3], faceCoord[0], edgeCoord[2], 
//       												edgeCoord[7], faceCoord[4], cellGlobal, faceCoord[3]);
//      		subContVol[4].volume = cellVolume - subContVol[0].volume - subContVol[1].volume 
//      								- subContVol[2].volume - subContVol[3].volume;
//       		break;
//      	case 9: // 3D, prism
//       		subContVol[0].volume = hexahedronVolume(subContVol[0].global, edgeCoord[0], faceCoord[0], edgeCoord[2], 
//       												edgeCoord[3], faceCoord[1], cellGlobal, faceCoord[3]);
//       		subContVol[1].volume = hexahedronVolume(subContVol[1].global, edgeCoord[1], faceCoord[0], edgeCoord[0], 
//       												edgeCoord[4], faceCoord[2], cellGlobal, faceCoord[1]);
//       		subContVol[2].volume = hexahedronVolume(subContVol[2].global, edgeCoord[2], faceCoord[0], edgeCoord[1], 
//       												edgeCoord[5], faceCoord[3], cellGlobal, faceCoord[2]);
//       		subContVol[3].volume = hexahedronVolume(edgeCoord[3], faceCoord[1], cellGlobal, faceCoord[3], subContVol[3].global, 
//       												edgeCoord[6], faceCoord[4], edgeCoord[8]);
//       		subContVol[4].volume = hexahedronVolume(edgeCoord[4], faceCoord[2], cellGlobal, faceCoord[1], subContVol[4].global, 
//       												edgeCoord[7], faceCoord[4], edgeCoord[6]);
//       		subContVol[5].volume = hexahedronVolume(edgeCoord[5], faceCoord[3], cellGlobal, faceCoord[2], subContVol[5].global, 
//       												edgeCoord[8], faceCoord[4], edgeCoord[7]);
//       		break;
//       	case 12: // 3D, hexahedron
//       		subContVol[0].volume = hexahedronVolume(subContVol[0].global, edgeCoord[8], faceCoord[4], edgeCoord[4], 
//       												edgeCoord[0], faceCoord[2], cellGlobal, faceCoord[0]);
//       		subContVol[1].volume = hexahedronVolume(subContVol[1].global, edgeCoord[5], faceCoord[4], edgeCoord[8], 
//       												edgeCoord[1], faceCoord[1], cellGlobal, faceCoord[2]);
//       		subContVol[2].volume = hexahedronVolume(subContVol[2].global, edgeCoord[4], faceCoord[4], edgeCoord[9], 
//       												edgeCoord[2], faceCoord[0], cellGlobal, faceCoord[3]);
//       		subContVol[3].volume = hexahedronVolume(subContVol[3].global, edgeCoord[9], faceCoord[4], edgeCoord[5], 
//       												edgeCoord[3], faceCoord[3], cellGlobal, faceCoord[1]);
//       		subContVol[4].volume = hexahedronVolume(edgeCoord[0], faceCoord[2], cellGlobal, faceCoord[0], 
//       												subContVol[4].global, edgeCoord[10], faceCoord[5], edgeCoord[6]);
//       		subContVol[5].volume = hexahedronVolume(edgeCoord[1], faceCoord[1], cellGlobal, faceCoord[2], 
//       												subContVol[5].global, edgeCoord[7], faceCoord[5], edgeCoord[10]);
//       		subContVol[6].volume = hexahedronVolume(edgeCoord[2], faceCoord[0], cellGlobal, faceCoord[3], 
//       												subContVol[6].global, edgeCoord[6], faceCoord[5], edgeCoord[11]);
//       		subContVol[7].volume = hexahedronVolume(edgeCoord[3], faceCoord[3], cellGlobal, faceCoord[1], 
//       												subContVol[7].global, edgeCoord[11], faceCoord[5], edgeCoord[7]);
//       		break;
//       	default:
//       		DUNE_THROW(NotImplemented, "FVElementGeometry for nEdges = " << nEdges);
//       	}
       	
//       	for (int k = 0; k < nNodes; k++)
//       		std::cout << "node " << k << ", volume = " << subContVol[k].volume 
//       			<< ", local = " << subContVol[k].local << ", global = " << subContVol[k].global << std::endl;

       	// fill sub control volume face data:
       	for (int k = 0; k < nEdges; k++) { // begin loop over edges / sub control volume faces
 			 int i = referenceElement.subEntity(k, dim-1, 0, dim);
 			 int j = referenceElement.subEntity(k, dim-1, 1, dim);
 			 if (nEdges == 4 && (i == 2 || j == 2)) 
 				 std::swap(i, j);
 			 subContVolFace[k].i = i;
 			 subContVolFace[k].j = j; 
 			 
 			 // calculate the local integration point and the face normal
 			 FV ipLocal;
 			 FV diffVec;
// 			 switch (dim) {
//			  case 1:
//				  subContVolFace[k].ipLocal = 0.5;
//				  subContVolFace[k].normal = 1.0;
//				  break;
//			  case 2:
				  ipLocal = referenceElement.position(k, dim-1) + cellLocal;
				  ipLocal *= 0.5;
				  subContVolFace[k].ipLocal = ipLocal;
				  diffVec = cellGlobal - edgeCoord[k];
				  subContVolFace[k].normal[0] = diffVec[1];
				  subContVolFace[k].normal[1] = -diffVec[0];
//				  break;
//			  case 3: 
//				  int leftFace; 
//				  int rightFace; 
//				  getFaceIndices(nNodes, k, leftFace, rightFace);
//				  ipLocal = referenceElement.position(k, dim-1) + cellLocal 
//				  				+ referenceElement.position(leftFace, 1)
//				  				+ referenceElement.position(rightFace, 1);
//				  ipLocal *= 0.25;
//				  subContVolFace[k].ipLocal = ipLocal;
//				  subContVolFace[k].normal = normalOfQuadrilateral3D(edgeCoord[k], faceCoord[rightFace], 
//						  											cellGlobal, faceCoord[leftFace]);
//				  break;
//			  }

 			 // get the global integration point and the Jacobian inverse
 			 subContVolFace[k].ipGlobal = geometry.global(ipLocal);
 			 FieldMatrix<DT,dim,dim> jacInvT = geometry.jacobianInverseTransposed(ipLocal);

			  
//			  std::cout << "SCV Face " << k << ", i = " << i << ", j = " << j 
//			  			<< ", ipLocal = " << ipLocal << ", ipGlobal = " << subContVolFace[k].ipGlobal << ", normal = " << subContVolFace[k].normal 
//			  			<< std::endl;
			  
 			 // calculate the shape function gradients 
 			 for (int node = 0; node < nNodes; node++) {
	        	  FV grad(0),temp;
	        	  for (int l = 0; l < dim; l++) 
	        		  temp[l] = sfs[node].evaluateDerivative(0, l, subContVolFace[k].ipLocal);
	        	  jacInvT.umv(temp, grad);
	        	  subContVolFace[k].grad[node] = grad;
	          }
       	} // end loop over edges / sub control volume faces

       	// fill boundary face data: 
       	IntersectionIterator endit = e.ileafend();
       	for (IntersectionIterator it = e.ileafbegin(); it != endit; ++it)
       		if (it.boundary())
       		{
       			int face = it.numberInSelf();
       			int nNodesOfFace = referenceElement.size(face, 1, dim);
       			for (int nodeInFace = 0; nodeInFace < nNodesOfFace; nodeInFace++)
       			{
       				int nodeInElement = referenceElement.subEntity(face, 1, nodeInFace, dim);
        			int bfIndex = boundaryFaceIndex(face, nodeInFace);
//           			switch (dim) {
//           			case 1: 
//           				boundaryFace[bfIndex].ipLocal = referenceElement.position(nodeInElement, dim);
//           				boundaryFace[bfIndex].area = 1.0;
//           				break;
//           			case 2: 
           				boundaryFace[bfIndex].ipLocal = referenceElement.position(nodeInElement, dim) 
           							+ referenceElement.position(face, 1);
           				boundaryFace[bfIndex].ipLocal *= 0.5;
           				//! \todo should be 0.5 instead of 0.25
           				boundaryFace[bfIndex].area = 0.25*it.intersectionGlobal().volume(); 
//           				break;
//           			case 3:
//           				int leftEdge;
//           				int rightEdge;
//           				getEdgeIndices(nNodes, face, nodeInElement, leftEdge, rightEdge);
//           				boundaryFace[bfIndex].ipLocal = referenceElement.position(nodeInElement, dim) 
//           												+ referenceElement.position(face, 1) 
//           												+ referenceElement.position(leftEdge, dim-1)
//           												+ referenceElement.position(rightEdge, dim-1);
//           				boundaryFace[bfIndex].ipLocal *= 0.25;
//           				//! \todo should be 1.0 instead of 0.5
//           				boundaryFace[bfIndex].area = 0.5*quadrilateralArea3D(subContVol[nodeInElement].global, 
//           												edgeCoord[rightEdge], faceCoord[face], edgeCoord[leftEdge]);
//           				break;
//           	       	default:
//           	       		DUNE_THROW(NotImplemented, "FVElementGeometry for dim = " << dim);
//           			}       				
       				boundaryFace[bfIndex].ipGlobal = geometry.global(boundaryFace[bfIndex].ipLocal);
       				
//      			  std::cout << "boundary face " << face << ", node = " << nodeInElement << ", ipLocal = " 
//      			  	<< boundaryFace[bfIndex].ipLocal << ", ipGlobal = " << boundaryFace[bfIndex].ipGlobal 
//      			  	<< ", area = " << boundaryFace[bfIndex].area << std::endl;

       			}
       		}
	}

};

}


#endif



