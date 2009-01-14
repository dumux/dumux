#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include"dune/common/mpihelper.hh" // An initializer of MPI
#include"dune/common/exceptions.hh" // We use exceptions
#include <dune/grid/uggrid.hh>
#include <dune/grid/albertagrid.hh>
//#include "gridcheck.cc"
//#include "checkgeometryinfather.cc"
//#include "checkintersectionit.cc"
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/subgrid/subgrid.hh>

using namespace Dune;

class VertexCompare {
public:
    template <class HostVertexPointer>
    bool operator() (const HostVertexPointer& p1, const HostVertexPointer& p2) const
    {
        char stringP1[100];
        sprintf(stringP1, "%e%e%e", p1->geometry().corner(0][0], p1->geometry()[0][1], p1->geometry()[0][2));
        char stringP2[100];
        sprintf(stringP2, "%e%e%e", p2->geometry().corner(0][0], p2->geometry()[0][1], p2->geometry()[0][2));
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
    std::vector<double> parameter;
    std::vector<unsigned> indexVertexVectorOnLine;
    std::vector<LinePointer> lineVectorOnLine;
    std::vector<unsigned> indexVertexVectorOutLine;
    std::vector<LinePointer> lineVectorOutLine;

    typedef std::vector<VertexOnLine<GridType> > VectorVertexOnLineType;
    bool boundary();
    double length(VectorVertexOnLineType& vectorVertexOnLine);
    template<class VectorNeighbor> Dune::FieldVector<double,GridType::dimension> normal(VectorNeighbor& vectorNeighbor, unsigned index);
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
        distVec = vectorVertexOnLine[indexVertexVectorOnLine[i]].nodePoint->geometry().corner(0] - nodePoint->geometry()[0);
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
    normalVec = vectorNeighbor[index].nodePoint->geometry().corner(0] - nodePoint->geometry()[0);
    double dist = normalVec.two_norm();
    normalVec *=1/dist;
    return(normalVec);
};

template <class GridType>
class VertexOutLine
{
    public:
    typedef typename GridType::template Codim<GridType::dimension>::EntityPointer VertexPointer;
    typedef typename GridType::template Codim<GridType::dimension-1>::EntityPointer LinePointer;
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


template<class G, class GP, class VerticeVectOnline, class VerticeVectOutline>                                      /*@\label{tc:tra0}@*/
void isolate (G& grid, GP& gridP, VerticeVectOnline& vertexVectorOnLine, VerticeVectOutline& vertexVectorOutLine)
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

                HostLinePointer linePointer = element.template entity<dim-1>(localIdLineonElement);

                // get local id of point on line
                int    locaIdPoint0onElement = Dune::ReferenceElements<ct,dim>::general(gt).subEntity(localIdLineonElement, dim-1, 0, dim);
                int    locaIdPoint1onElement = Dune::ReferenceElements<ct,dim>::general(gt).subEntity(localIdLineonElement, dim-1, 1, dim);
                // get vertex pointers
                HostVertexPointer vertexPoint0 = element.template entity<dim>(locaIdPoint0onElement);
                HostVertexPointer vertexPoint1 = element.template entity<dim>(locaIdPoint1onElement);
                // get parameters for nodes
                std::vector<double>& param0 = gridP.parameters(*vertexPoint0);
                std::vector<double>& param1 = gridP.parameters(*vertexPoint1);

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
                        mapVertexIterator current  = mapVertexToIndexOnLine.find(vertexPoint0);
                        if( current == mapVertexToIndexOnLine.end() )
                        {
                            VertexOnLine<G> vertex(vertexPoint0, param0);
                            mapVertexToIndexOnLine[vertexPoint0] = (indexMapVerticesOnLine++);
                            vertexVectorOnLine.push_back(vertex);
                            addpoint0onLine = true;
                        }
                    }
                    else
                    {
                        mapVertexIterator current  = mapVertexToIndexOutLine.find(vertexPoint0);
                        if( current == mapVertexToIndexOutLine.end() )
                        {
                            mapVertexToIndexOutLine[vertexPoint0] = (indexMapVerticesOutLine++);
                            vertexVectorOutLine.push_back(vertexPoint0);
                            addpoint0outLine = true;
                        }
                    }

                    if(NodeParameter[1]!=0)
                    {
                        mapVertexIterator current  = mapVertexToIndexOnLine.find(vertexPoint1);
                        if( current == mapVertexToIndexOnLine.end() )
                        {
                            VertexOnLine<G> vertex(vertexPoint1, param1);
                            mapVertexToIndexOnLine[vertexPoint1] = (indexMapVerticesOnLine++);
                            vertexVectorOnLine.push_back(vertex);
                            addpoint1onLine = true;
                        }
                    }
                    else
                    {
                        mapVertexIterator current  = mapVertexToIndexOutLine.find(vertexPoint1);
                        if( current == mapVertexToIndexOutLine.end() )
                        {
                            mapVertexToIndexOutLine[vertexPoint1] = (indexMapVerticesOutLine++);
                            vertexVectorOutLine.push_back(vertexPoint1);
                            addpoint1outLine = true;
                        }
                    }
                    // fill indexVertexVectorOnLine Vector and lineVectorOnline Vector in vertexVectorOnLine Vector
                    if ( NodeParameter[0]!=0  && NodeParameter[1]!=0 )
                    {
                        mapVertexIterator current0  = mapVertexToIndexOnLine.find(vertexPoint0);
                        mapVertexIterator current1  = mapVertexToIndexOnLine.find(vertexPoint1);
                        unsigned indexVectorPoint0 = current0->second;
                        unsigned indexVectorPoint1 = current1->second;
                        if (addpoint0onLine || addpoint1onLine)
                        {
                            vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOnLine.push_back(indexVectorPoint1);
                            vertexVectorOnLine[indexVectorPoint1].indexVertexVectorOnLine.push_back(indexVectorPoint0);
                            vertexVectorOnLine[indexVectorPoint0].lineVectorOnLine.push_back(linePointer);
                            vertexVectorOnLine[indexVectorPoint1].lineVectorOnLine.push_back(linePointer);
                        }
                        else if(!addpoint0onLine && !addpoint1onLine) //fill the missing part where two points are there but no neighborhoot is stored
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
                                vertexVectorOnLine[indexVectorPoint0].lineVectorOnLine.push_back(linePointer);
                                vertexVectorOnLine[indexVectorPoint1].lineVectorOnLine.push_back(linePointer);
                            }
                        }
                    }
                    // end filling indexVertexVectorOnLine Vector and lineVectorOnline Vector in vertexVectorOnLine Vector

                    // fill indexVertexVectorOutLine Vector and lineVectorOutline Vector in vertexVectorOutLine Vector and in vertexVectorOnLine Vector
                    if ( !(NodeParameter[0]!=0  && NodeParameter[1]!=0) )
                    {
                        if (NodeParameter[0]!=0)
                        {
                            mapVertexIterator current0  = mapVertexToIndexOnLine.find(vertexPoint0);
                            mapVertexIterator current1  = mapVertexToIndexOutLine.find(vertexPoint1);
                            unsigned indexVectorPoint0 = current0->second;
                            unsigned indexVectorPoint1 = current1->second;
                            if (addpoint0onLine || addpoint1outLine)
                            {
                                vertexVectorOnLine[indexVectorPoint0].indexVertexVectorOutLine.push_back(indexVectorPoint1);
                                vertexVectorOutLine[indexVectorPoint1].indexVertexVectorOnLine.push_back(indexVectorPoint0);
                                vertexVectorOnLine[indexVectorPoint0].lineVectorOutLine.push_back(linePointer);
                                vertexVectorOutLine[indexVectorPoint1].lineVectorOutLine.push_back(linePointer);
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
                                    vertexVectorOnLine[indexVectorPoint0].lineVectorOutLine.push_back(linePointer);
                                    vertexVectorOutLine[indexVectorPoint1].lineVectorOutLine.push_back(linePointer);
                                }
                            }
                        }
                        else if (NodeParameter[1]!=0)
                        {
                            mapVertexIterator current0  = mapVertexToIndexOutLine.find(vertexPoint0);
                            mapVertexIterator current1  = mapVertexToIndexOnLine.find(vertexPoint1);
                            unsigned indexVectorPoint0 = current0->second;
                            unsigned indexVectorPoint1 = current1->second;
                            if (addpoint0outLine || addpoint1onLine)
                            {
                                vertexVectorOutLine[indexVectorPoint0].indexVertexVectorOnLine.push_back(indexVectorPoint1);
                                vertexVectorOnLine[indexVectorPoint1].indexVertexVectorOutLine.push_back(indexVectorPoint0);
                                vertexVectorOutLine[indexVectorPoint0].lineVectorOutLine.push_back(linePointer);
                                vertexVectorOnLine[indexVectorPoint1].lineVectorOutLine.push_back(linePointer);
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
                                    vertexVectorOutLine[indexVectorPoint0].lineVectorOutLine.push_back(linePointer);
                                    vertexVectorOnLine[indexVectorPoint1].lineVectorOutLine.push_back(linePointer);
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

int main(int argc, char** argv)
{
  try{
    // set up the grid
    const int dim = 3;
    typedef double NumberType;

//    typedef Dune::ALUSimplexGrid<dim,dim> GridType;
//    typedef Dune::UGGrid<dim> GridType;
    typedef Dune::ALUCubeGrid<dim,dim> GridType;

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( "grids/rectangle_inner_manual_3d.dgf" );

    // grid reference
    GridType& grid = *gridPtr;

    typedef VertexOnLine<GridType> VertexType1;
    typedef VertexOutLine<GridType> VertexType2;

    std::vector<VertexType1> vertexVectorOnLine;
    std::vector<VertexType2> vertexVectorOutLine;

    isolate(grid, gridPtr, vertexVectorOnLine, vertexVectorOutLine);

    std::cout << "number of vertices on line = "<< vertexVectorOnLine.size() << std::endl;
    std::cout << "number of vertices out line = "<< vertexVectorOutLine.size() << std::endl;

    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
    {
        std::cout << "vertice on line coord: " <<vertexVectorOnLine[k].nodePoint->geometry().corner(0) << std::endl;
        std::cout << "       vertice on line: " << k << std::endl;
        for (unsigned n = 0; n < vertexVectorOnLine[k].parameter.size(); n++)
        {
            std::cout << "       parameters: " << vertexVectorOnLine[k].parameter[n] << std::endl;
        }
        std::cout << "           boundary: " << vertexVectorOnLine[k].boundary() << std::endl;
        std::cout << "             length: " << vertexVectorOnLine[k].length(vertexVectorOnLine) << std::endl;
        std::cout << "             on line size: " << vertexVectorOnLine[k].lineVectorOnLine.size() << std::endl;
        for (unsigned m = 0; m < vertexVectorOnLine[k].indexVertexVectorOnLine.size(); m++)
        {
            std::cout << "                neighbour vertices on line: " << vertexVectorOnLine[k].indexVertexVectorOnLine[m] << std::endl;
            std::cout << "                normal vector: " << vertexVectorOnLine[k].normal(vertexVectorOnLine, vertexVectorOnLine[k].indexVertexVectorOnLine[m]) << std::endl;
        }
        std::cout << "             out line size: " << vertexVectorOnLine[k].lineVectorOutLine.size() << std::endl;
        for (unsigned m = 0; m < vertexVectorOnLine[k].indexVertexVectorOutLine.size(); m++)
        {
            std::cout << "                neighbour vertices out line: " << vertexVectorOnLine[k].indexVertexVectorOutLine[m] << std::endl;
            std::cout << "                normal vector: " << vertexVectorOnLine[k].normal(vertexVectorOutLine, vertexVectorOnLine[k].indexVertexVectorOutLine[m]) << std::endl;
        }
    }

    for (unsigned k = 0; k < vertexVectorOutLine.size(); k++)
    {
        std::cout << "vertice out line coord: " <<vertexVectorOutLine[k].nodePoint->geometry().corner(0) << std::endl;
        std::cout << "       vertice out line: " << k << std::endl;
        std::cout << "          out line size: " << vertexVectorOutLine[k].lineVectorOutLine.size() << std::endl;
        for (unsigned m = 0; m < vertexVectorOutLine[k].indexVertexVectorOnLine.size(); m++)
        {
            std::cout << "             neighbour vertices on line: " << vertexVectorOutLine[k].indexVertexVectorOnLine[m] << std::endl;
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
