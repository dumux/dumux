#include<dune/grid/common/mcmgmapper.hh>

class GlobalNodeIdCompare {
public:
	bool operator() (const unsigned& globalId1, const unsigned& globalId2) const
	{
		return (globalId1 < globalId2);
	}
};

class VertexCompare {
public:
	template <class HostVertexPointer>
	bool operator() (const HostVertexPointer& p1, const HostVertexPointer& p2) const
	{
		char stringP1[100];
		sprintf(stringP1, "%e%e%e", p1->geometry()[0][0], p1->geometry()[0][1], p1->geometry()[0][2]);
		char stringP2[100];
		sprintf(stringP2, "%e%e%e", p2->geometry()[0][0], p2->geometry()[0][1], p2->geometry()[0][2]);
		return strcmp(stringP1, stringP2) < 0;
	}
};

template <class GridType>
class VertexOnLine
{
	public:
	typedef typename GridType::template Codim<GridType::dimension>::EntityPointer VertexPointer;
	typedef typename GridType::template Codim<GridType::dimension-1>::EntityPointer LinePointer;
    VertexPointer nodePoint;
    unsigned globalId;
    std::vector<double> parameter;
    std::vector<unsigned> indexVertexVectorOnLine;
    std::vector<LinePointer> lineVectorOnLine;
    std::vector<unsigned> indexVertexVectorOutLine;
    std::vector<LinePointer> lineVectorOutLine;

    typedef std::vector<VertexOnLine<GridType> > VectorVertexOnLineType;
    bool boundary();
    double length(VectorVertexOnLineType& vectorVertexOnLine);
    Dune::FieldVector<double,GridType::dimension> unitPD(VectorVertexOnLineType& vectorNeighbor, unsigned index); // PD means positive direction
    template<class VectorNeighbor> Dune::FieldVector<double,GridType::dimension> normal(VectorNeighbor& vectorNeighbor, unsigned index);
    Dune::FieldVector<double,GridType::dimension> unitPDBF(VectorVertexOnLineType& vectorNeighbor); // BF means boundary face
    Dune::FieldVector<double,GridType::dimension> normalBF(VectorVertexOnLineType& vectorNeighbor);

    VertexOnLine(VertexPointer node0, std::vector<double> param)
	: nodePoint(node0)
	{
    	parameter = param;
	}

    VertexOnLine()
	{
	}

};

template <class GridType>
bool VertexOnLine<GridType>::boundary()
{
	if(indexVertexVectorOnLine.size()==1)
		return true;
	return false;
};

template <class GridType>
double VertexOnLine<GridType>::length(VectorVertexOnLineType& vectorVertexOnLine)
{
	double l=0.0;
	for(unsigned i=0; i < indexVertexVectorOnLine.size(); i++)
	{
		Dune::FieldVector<double,GridType::dimension> distVec;
		distVec = vectorVertexOnLine[indexVertexVectorOnLine[i]].nodePoint->geometry()[0] - nodePoint->geometry()[0];
		double dist = distVec.two_norm();
		l += dist/2.0;
	}
	return l;
};

template <class GridType>
template <class VectorNeighbor>
Dune::FieldVector<double,GridType::dimension> VertexOnLine<GridType>::normal(VectorNeighbor& vectorNeighbor, unsigned index)
{
	Dune::FieldVector<double,GridType::dimension> normalVec;
	normalVec = vectorNeighbor[index].nodePoint->geometry()[0] - nodePoint->geometry()[0];
	double dist = normalVec.two_norm();
	normalVec *=1/dist;
	return(normalVec);
};

// this unitPD only makes sense as long as the pipe is along the coordinate axis
template <class GridType>
Dune::FieldVector<double,GridType::dimension> VertexOnLine<GridType>::unitPD(VectorVertexOnLineType& vectorNeighbor, unsigned index)
{
	Dune::FieldVector<double,GridType::dimension> unitPozitiveDirectionVector;
	unitPozitiveDirectionVector = vectorNeighbor[index].nodePoint->geometry()[0] - nodePoint->geometry()[0];
	for (unsigned k = 0; k < GridType::dimension; k++)
	    {
			if(unitPozitiveDirectionVector[k]<0) unitPozitiveDirectionVector[k] *= -1;
	    }
	double dist = unitPozitiveDirectionVector.two_norm();
	unitPozitiveDirectionVector *=1/dist;
	return(unitPozitiveDirectionVector);
};

template <class GridType>
Dune::FieldVector<double,GridType::dimension> VertexOnLine<GridType>::normalBF(VectorVertexOnLineType& vectorNeighbor)
{
	if(boundary())
	{
		Dune::FieldVector<double,GridType::dimension> normalVec;
		normalVec = nodePoint->geometry()[0] - vectorNeighbor[indexVertexVectorOnLine[0]].nodePoint->geometry()[0];
		double dist = normalVec.two_norm();
		normalVec *=1/dist;
		return(normalVec);
	}
	else
	{
    	Dune::Exception exception;
    	exception.message("this node is not a boundary !!!!");
    	throw exception;
	}
};

// this unitPD only makes sense as long as the pipe is along the coordinate axis
template <class GridType>
Dune::FieldVector<double,GridType::dimension> VertexOnLine<GridType>::unitPDBF(VectorVertexOnLineType& vectorNeighbor)
{
	if(boundary())
	{
		Dune::FieldVector<double,GridType::dimension> unitPozitiveDirectionVector;
		unitPozitiveDirectionVector = nodePoint->geometry()[0] - vectorNeighbor[indexVertexVectorOnLine[0]].nodePoint->geometry()[0];
		for (unsigned k = 0; k < GridType::dimension; k++)
		    {
				if(unitPozitiveDirectionVector[k]<0) unitPozitiveDirectionVector[k] *= -1;
		    }
		double dist = unitPozitiveDirectionVector.two_norm();
		unitPozitiveDirectionVector *=1/dist;
		return(unitPozitiveDirectionVector);
	}
	else
	{
    	Dune::Exception exception;
    	exception.message("this node is not a boundary !!!!");
    	throw exception;
	}
};

template <class GridType>
class VertexOutLine
{
	public:
	typedef typename GridType::template Codim<GridType::dimension>::EntityPointer VertexPointer;
	typedef typename GridType::template Codim<GridType::dimension-1>::EntityPointer LinePointer;
	unsigned globalId;
    VertexPointer nodePoint;
    std::vector<unsigned> indexVertexVectorOnLine;
    std::vector<LinePointer> lineVectorOutLine;
    VertexOutLine(VertexPointer node0)
	: nodePoint(node0)
	{
	}

    VertexOutLine()
	{
	}

};

template<class G, class GP, class VM, class MapperGlobalNodeIDtoPipeNodeOnOutlineIndexType, class VerticeVectOnline, class VerticeVectOutline>                                      /*@\label{tc:tra0}@*/
void isolate (G& grid, GP& gridP, VM& vertexmapper, MapperGlobalNodeIDtoPipeNodeOnOutlineIndexType& mapGlobalNodeIDtoPipeNodeOnlineIndex, MapperGlobalNodeIDtoPipeNodeOnOutlineIndexType& mapGlobalNodeIDtoPipeNodeOutlineIndex,  VerticeVectOnline& vertexVectorOnLine, VerticeVectOutline& vertexVectorOutLine)
{
	// first we extract the dimensions of the grid
	const int dim = G::dimension;                        /*@\label{tc:dim}@*/
	typedef typename G::ctype ct;                        /*@\label{tc:ct}@*/

	typedef typename G::template Codim<0>::LeafIterator HostElementIterator;
	typedef typename G::template Codim<dim>::EntityPointer HostVertexPointer;
	typedef typename G::template Codim<dim-1>::EntityPointer HostLinePointer;
	typedef typename G::template Codim<0>::LeafIntersectionIterator HostIntersectionIterator;
	typedef typename G::Traits::template Codim<0>::Entity Entity;

    std::map<HostVertexPointer, unsigned, VertexCompare> mapVertexToIndexOnLine;
    std::map<HostVertexPointer, unsigned, VertexCompare> mapVertexToIndexOutLine;

	typedef typename std::map<HostVertexPointer, unsigned>::const_iterator mapVertexIterator;

	unsigned indexMapVerticesOnLine = 0;
	unsigned indexMapVerticesOutLine = 0;
	for (HostElementIterator it = grid.template leafbegin<0>(); it!=grid.template leafend<0>(); ++it)
	{
		// cell geometry type
		Dune::GeometryType gt = it->geometry().type();
		const Entity& element = *it;
		for (HostIntersectionIterator is = it->ileafbegin(); is!= it->ileafend(); ++is)
		{
			// get geometry type of face
			Dune::GeometryType gtf = is.intersectionSelfLocal().type();
			// get local id of line on element
			int localIdFaceonElement = is.numberInSelf();
		    int numLeaf = Dune::ReferenceElements<ct,dim>::general(gt).size(localIdFaceonElement,dim-2,dim-1);

		    // start line loop
			for (int i=0; i<numLeaf; i++ )
			{
				int localIdLineonElement= Dune::ReferenceElements<ct,dim>::general(gt).subEntity(localIdFaceonElement, dim-2, i, dim-1);

//				HostLinePointer linePointer = element.template entity<dim-1>(localIdLineonElement);

				// get local id of point on line
				int	locaIdPoint0onElement = Dune::ReferenceElements<ct,dim>::general(gt).subEntity(localIdLineonElement, dim-1, 0, dim);
				int	locaIdPoint1onElement = Dune::ReferenceElements<ct,dim>::general(gt).subEntity(localIdLineonElement, dim-1, 1, dim);
				// get vertex pointers
//				HostVertexPointer vertexPoint0 = element.template entity<dim>(locaIdPoint0onElement);
//				HostVertexPointer vertexPoint1 = element.template entity<dim>(locaIdPoint1onElement);
				// get parameters for nodes
				std::vector<double>& param0 = gridP.parameters(*element.template entity<dim>(locaIdPoint0onElement));
				std::vector<double>& param1 = gridP.parameters(*element.template entity<dim>(locaIdPoint1onElement));

				double NodeParameter[1];
				NodeParameter[0]= param0[0];
				NodeParameter[1]= param1[0];

				bool addpoint0onLine = false;
				bool addpoint1onLine = false;
				bool addpoint0outLine = false;
				bool addpoint1outLine = false;
				if ( NodeParameter[0]!=0  || NodeParameter[1]!=0)
				{
					if(NodeParameter[0]!=0)
					{
						mapVertexIterator current  = mapVertexToIndexOnLine.find(element.template entity<dim>(locaIdPoint0onElement));
						if( current == mapVertexToIndexOnLine.end() ) //check wheter the point is already added or not
						{
							VertexOnLine<G> vertex(element.template entity<dim>(locaIdPoint0onElement), param0); //create the vertex class
							vertex.globalId = vertexmapper.template map<dim>(*it, locaIdPoint0onElement); // find the global id of the point on the big grid
							mapGlobalNodeIDtoPipeNodeOnlineIndex[vertex.globalId] = indexMapVerticesOnLine; // map global id of point to index online
							mapVertexToIndexOnLine[element.template entity<dim>(locaIdPoint0onElement)] = (indexMapVerticesOnLine++); // map vertex pointer to index online
							vertexVectorOnLine.push_back(vertex); // add vertex object to vertexVectorOnline
							addpoint0onLine = true;
						}
					}
					else
					{
						mapVertexIterator current  = mapVertexToIndexOutLine.find(element.template entity<dim>(locaIdPoint0onElement));
						if( current == mapVertexToIndexOutLine.end() )
						{
							VertexOutLine<G> vertex(element.template entity<dim>(locaIdPoint0onElement));
							vertex.globalId = vertexmapper.template map<dim>(*it, locaIdPoint0onElement);
							mapGlobalNodeIDtoPipeNodeOutlineIndex[vertex.globalId] = indexMapVerticesOutLine;
							mapVertexToIndexOutLine[element.template entity<dim>(locaIdPoint0onElement)] = (indexMapVerticesOutLine++);
							vertexVectorOutLine.push_back(vertex);
							addpoint0outLine = true;
						}
					}

					if(NodeParameter[1]!=0)
					{
						mapVertexIterator current  = mapVertexToIndexOnLine.find(element.template entity<dim>(locaIdPoint1onElement));
						if( current == mapVertexToIndexOnLine.end() )
						{
							VertexOnLine<G> vertex(element.template entity<dim>(locaIdPoint1onElement), param1);
							vertex.globalId = vertexmapper.template map<dim>(*it, locaIdPoint1onElement);
							mapGlobalNodeIDtoPipeNodeOnlineIndex[vertex.globalId] = indexMapVerticesOnLine;
							mapVertexToIndexOnLine[element.template entity<dim>(locaIdPoint1onElement)] = (indexMapVerticesOnLine++);
							vertexVectorOnLine.push_back(vertex);
							addpoint1onLine = true;
						}
					}
					else
					{
						mapVertexIterator current  = mapVertexToIndexOutLine.find(element.template entity<dim>(locaIdPoint1onElement));
						if( current == mapVertexToIndexOutLine.end() )
						{
							VertexOutLine<G> vertex(element.template entity<dim>(locaIdPoint1onElement));
							vertex.globalId = vertexmapper.template map<dim>(*it, locaIdPoint1onElement);
							mapGlobalNodeIDtoPipeNodeOutlineIndex[vertex.globalId] = indexMapVerticesOutLine;
							mapVertexToIndexOutLine[element.template entity<dim>(locaIdPoint1onElement)] = (indexMapVerticesOutLine++);
							vertexVectorOutLine.push_back(vertex);
							addpoint1outLine = true;
						}
					}
					// fill neighborhood
					// fill indexVertexVectorOnLine Vector and lineVectorOnline Vector in vertexVectorOnLine Vector
					if ( NodeParameter[0]!=0  && NodeParameter[1]!=0 )
					{
						mapVertexIterator current0  = mapVertexToIndexOnLine.find(element.template entity<dim>(locaIdPoint0onElement));
						mapVertexIterator current1  = mapVertexToIndexOnLine.find(element.template entity<dim>(locaIdPoint1onElement));
						unsigned indexVectorPoint0 = current0->second;
						unsigned indexVectorPoint1 = current1->second;
						if (addpoint0onLine || addpoint1onLine) // if at least one of the point is new
						{
							vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOnLine.push_back(indexVectorPoint1);
							vertexVectorOnLine[indexVectorPoint1].indexVertexVectorOnLine.push_back(indexVectorPoint0);
							vertexVectorOnLine[indexVectorPoint0].lineVectorOnLine.push_back(element.template entity<dim-1>(localIdLineonElement));
							vertexVectorOnLine[indexVectorPoint1].lineVectorOnLine.push_back(element.template entity<dim-1>(localIdLineonElement));
						}
						else if(!addpoint0onLine && !addpoint1onLine) // if none of the points are new (fill the missing part where two points are there but no neighborhood is stored)
						{
							bool neighbourhoodAdded = false;
							for (unsigned i=0 ; i < vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOnLine.size(); i++)
							{
								if(vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOnLine[i] == indexVectorPoint1)
								{
									neighbourhoodAdded = true;
									break;
								}
							}
							if (!neighbourhoodAdded)
							{
								vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOnLine.push_back(indexVectorPoint1);
								vertexVectorOnLine[indexVectorPoint1].indexVertexVectorOnLine.push_back(indexVectorPoint0);
								vertexVectorOnLine[indexVectorPoint0].lineVectorOnLine.push_back(element.template entity<dim-1>(localIdLineonElement));
								vertexVectorOnLine[indexVectorPoint1].lineVectorOnLine.push_back(element.template entity<dim-1>(localIdLineonElement));
							}
						}
					}
					// end filling indexVertexVectorOnLine Vector and lineVectorOnline Vector in vertexVectorOnLine Vector

					// fill indexVertexVectorOutLine Vector and lineVectorOutline Vector in vertexVectorOutLine Vector and in vertexVectorOnLine Vector
					if ( !(NodeParameter[0]!=0  && NodeParameter[1]!=0) )
					{
						if (NodeParameter[0]!=0)
						{
							mapVertexIterator current0  = mapVertexToIndexOnLine.find(element.template entity<dim>(locaIdPoint0onElement));
							mapVertexIterator current1  = mapVertexToIndexOutLine.find(element.template entity<dim>(locaIdPoint1onElement));
							unsigned indexVectorPoint0 = current0->second;
							unsigned indexVectorPoint1 = current1->second;
							if (addpoint0onLine || addpoint1outLine)
							{
								vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOutLine.push_back(indexVectorPoint1);
								vertexVectorOutLine[indexVectorPoint1].indexVertexVectorOnLine.push_back(indexVectorPoint0);
								vertexVectorOnLine[indexVectorPoint0].lineVectorOutLine.push_back(element.template entity<dim-1>(localIdLineonElement));
								vertexVectorOutLine[indexVectorPoint1].lineVectorOutLine.push_back(element.template entity<dim-1>(localIdLineonElement));
							}
							else if(!addpoint0onLine && !addpoint1onLine)
							{
								bool neighbourhoodAdded = false;
								for (unsigned i=0 ; i < vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOutLine.size(); i++)
								{
									if(vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOutLine[i] == indexVectorPoint1)
									{
										neighbourhoodAdded = true;
										break;
									}
								}
								if (!neighbourhoodAdded)
								{
									vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOutLine.push_back(indexVectorPoint1);
									vertexVectorOutLine[indexVectorPoint1].indexVertexVectorOnLine.push_back(indexVectorPoint0);
									vertexVectorOnLine[indexVectorPoint0].lineVectorOutLine.push_back(element.template entity<dim-1>(localIdLineonElement));
									vertexVectorOutLine[indexVectorPoint1].lineVectorOutLine.push_back(element.template entity<dim-1>(localIdLineonElement));
								}
							}
						}
						else if (NodeParameter[1]!=0)
						{
							mapVertexIterator current0  = mapVertexToIndexOutLine.find(element.template entity<dim>(locaIdPoint0onElement));
							mapVertexIterator current1  = mapVertexToIndexOnLine.find(element.template entity<dim>(locaIdPoint1onElement));
							unsigned indexVectorPoint0 = current0->second;
							unsigned indexVectorPoint1 = current1->second;
							if (addpoint0outLine || addpoint1onLine)
							{
								vertexVectorOutLine[indexVectorPoint0].indexVertexVectorOnLine.push_back(indexVectorPoint1);
								vertexVectorOnLine[indexVectorPoint1].indexVertexVectorOutLine.push_back(indexVectorPoint0);
								vertexVectorOutLine[indexVectorPoint0].lineVectorOutLine.push_back(element.template entity<dim-1>(localIdLineonElement));
								vertexVectorOnLine[indexVectorPoint1].lineVectorOutLine.push_back(element.template entity<dim-1>(localIdLineonElement));
							}
							else if(!addpoint0onLine && !addpoint1onLine)
							{
								bool neighbourhoodAdded = false;
								for (unsigned i=0 ; i < vertexVectorOnLine[indexVectorPoint1].indexVertexVectorOutLine.size(); i++)
								{
									if(vertexVectorOnLine[indexVectorPoint1].indexVertexVectorOutLine[i] == indexVectorPoint0)
									{
										neighbourhoodAdded = true;
										break;
									}
								}
								if (!neighbourhoodAdded)
								{
									vertexVectorOutLine[indexVectorPoint0].indexVertexVectorOnLine.push_back(indexVectorPoint1);
									vertexVectorOnLine[indexVectorPoint1].indexVertexVectorOutLine.push_back(indexVectorPoint0);
									vertexVectorOutLine[indexVectorPoint0].lineVectorOutLine.push_back(element.template entity<dim-1>(localIdLineonElement));
									vertexVectorOnLine[indexVectorPoint1].lineVectorOutLine.push_back(element.template entity<dim-1>(localIdLineonElement));
								}
							}
						}
					}
					// end filling indexVertexVectorOutLine Vector and lineVectorOutline Vector in vertexVectorOutLine Vector and in vertexVectorOnLine Vector
				}
			} // end line loop
		} //end intersection loop
	}
}
