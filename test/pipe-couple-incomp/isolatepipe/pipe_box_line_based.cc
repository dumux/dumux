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
        sprintf(stringP1, "%e%e", p1->geometry().corner(0][0], p1->geometry()[0][1));
        char stringP2[100];
        sprintf(stringP2, "%e%e", p2->geometry().corner(0][0], p2->geometry()[0][1));
        return strcmp(stringP1, stringP2) < 0;
    }
};


template <class VertexPointer>
struct Vertice
{
    VertexPointer nodePointer;
    bool boundary;
    std::vector<int> neighbourEdges;

    Vertice(VertexPointer nodeP)
        : nodePointer(nodeP)
    {
        boundary = true;
    }

    Vertice()
    {

    }
};

template <class GridType>
class Edge
{
public:
    typedef typename GridType::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<GridType::dimension>::EntityPointer VertexPointer;
    IntersectionIterator isIt;
    std::vector<Vertice<VertexPointer> > vertices;
    Edge(IntersectionIterator isI, VertexPointer node0, VertexPointer node1)
        : isIt(isI), vertices(2, Vertice<VertexPointer>(node0))
    {
        vertices[1] = *(new Vertice<VertexPointer>(node1));
    }

    Edge()
    {

    }

};

template<class G, class GP, class EdgeV, class MapperVetex>                                      /*@\label{tc:tra0}@*/
void isolate (G& grid, GP& gridP, EdgeV& edgeVector, MapperVetex& vertexToIndex)
{
    // first we extract the dimensions of the grid
    const int dim = G::dimension;                        /*@\label{tc:dim}@*/
    typedef typename G::ctype ct;                        /*@\label{tc:ct}@*/

    typedef typename G::template Codim<0>::LeafIterator HostElementIterator;
    typedef typename G::template Codim<dim>::EntityPointer HostVertexPointer;
    typedef typename G::template Codim<0>::LeafIntersectionIterator HostIntersectionIterator;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    //    typedef Edge<HostIntersectionIterator, HostVertexPointer> Edge;
    typedef typename SubGrid<dim,G>::template Codim<0>::LevelIterator ElementIterator;

    SubGrid<dim,G> subGrid(grid);

    //    std::map<HostVertexPointer, int, VertexCompare> vertexToIndex;
    //    std::vector<Edge> edgeVector;

    typedef typename std::map<HostVertexPointer, int>::const_iterator mapIterator;

    int vertices = 0;
    // create subgrid
    {
        subGrid.createBegin();
        for (HostElementIterator it = grid.template leafbegin<0>(); it!=grid.template leafend<0>(); ++it)
        {
            // cell geometry type
            Dune::GeometryType gt = it->geometry().type();

            const Entity& element = *it;
            bool addElem;
            for (HostIntersectionIterator is = it->ileafbegin(); is!= it->ileafend(); ++is)
            {
                // get geometry type of face
                Dune::GeometryType gtf = is.intersectionSelfLocal().type();
                // get local id of line on element
                int localIdLineonElement = is.numberInSelf();
                // get local id of point on line
                int    locaIdPoint0onElement = Dune::ReferenceElements<ct,dim>::general(gt).subEntity(localIdLineonElement, dim-1, 0, dim);
                int    locaIdPoint1onElement = Dune::ReferenceElements<ct,dim>::general(gt).subEntity(localIdLineonElement, dim-1, 1, dim);
                // get vertex pointers
                HostVertexPointer leftVertex = element.template entity<dim>(locaIdPoint0onElement);
                HostVertexPointer rightVertex = element.template entity<dim>(locaIdPoint1onElement);
                // get parameters for nodes
                std::vector<double>& paramleft = gridP.parameters(*leftVertex);
                std::vector<double>& paramright = gridP.parameters(*rightVertex);
                double NodeParameter[1];
                NodeParameter[0]= paramleft[0];
                NodeParameter[1]= paramright[0];
                if ( NodeParameter[0]!=0  && NodeParameter[1]!=0 && !(NodeParameter[0]>0 && NodeParameter[1]>0 && NodeParameter[0]!= NodeParameter[1]) )
                {
                    std::cout << "Mapper size " << vertexToIndex.size() << std::endl;
                    std::cout << "first param: " << NodeParameter[0] << "second param: " << NodeParameter[1] << std::endl;
                    // add element to subgrid
                    addElem = true;
                    // add line to line index
                    int oldVertices = vertices;
                    mapIterator cur0  = vertexToIndex.find(leftVertex);
                    if( cur0 == vertexToIndex.end() )
                        vertexToIndex[leftVertex] = (vertices++);
                    mapIterator cur1  = vertexToIndex.find(rightVertex);
                    if( cur1 == vertexToIndex.end() )
                        vertexToIndex[rightVertex] = (vertices++);
                    // if new vertice(s) is found add edge to edge vector
                    if (oldVertices < vertices)
                    {
                        Edge<G> edge(is, leftVertex, rightVertex);
                        edgeVector.push_back(edge);
                    }
                    else // add missing edges (when two nodes are already in vertexmapper but the line is not in edgevector)
                    {
                        bool found = false;
                        for (unsigned k = 0; k < edgeVector.size(); k++)
                        {
                            int nodeIdEv[2];
                            int nodeId[2];
                            mapIterator cur0Ev  = vertexToIndex.find(edgeVector[k].vertices[0].nodePointer);
                            mapIterator cur1Ev  = vertexToIndex.find(edgeVector[k].vertices[1].nodePointer);
                            nodeIdEv[0] = cur0Ev->second;
                            nodeIdEv[1] = cur1Ev->second;
                            mapIterator curLeft  = vertexToIndex.find(leftVertex);
                            mapIterator curRight  = vertexToIndex.find(rightVertex);
                            nodeId[0] = curLeft->second;
                            nodeId[1] = curRight->second;

                            if ( (nodeIdEv[0]==nodeId[0] && nodeIdEv[1]==nodeId[1]) || (nodeIdEv[0]==nodeId[1] && nodeIdEv[1]==nodeId[0]) )
                            {
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                        {
                            Edge<G> edge(is, leftVertex, rightVertex);
                            edgeVector.push_back(edge);
                        }
                    }
                    std::cout << "Mapper size " << vertexToIndex.size() << std::endl;
                }
            }
            if (addElem) subGrid.add(it);
        }
        subGrid.createEnd();

        // fill edge structure
        for (unsigned k = 0; k < edgeVector.size(); k++)
        {
            int nodeId[2];
            mapIterator cur0  = vertexToIndex.find(edgeVector[k].vertices[0].nodePointer);
            mapIterator cur1  = vertexToIndex.find(edgeVector[k].vertices[1].nodePointer);
            nodeId[0] = cur0->second;
            nodeId[1] = cur1->second;
            for (unsigned j = k+1; j < edgeVector.size(); j++)
            {
                int nodeIdNb[2];
                mapIterator cur0Nb  = vertexToIndex.find(edgeVector[j].vertices[0].nodePointer);
                mapIterator cur1Nb  = vertexToIndex.find(edgeVector[j].vertices[1].nodePointer);
                nodeIdNb[0] = cur0Nb->second;
                nodeIdNb[1] = cur1Nb->second;

                for (int n = 0; n < 2; n++)
                    for (int m = 0; m < 2; m++)
                    {
                        if (nodeId[n] == nodeIdNb[m])
                        {
                            edgeVector[k].vertices[n].boundary = false;
                            edgeVector[j].vertices[m].boundary = false;
                            edgeVector[k].vertices[n].neighbourEdges.push_back(j);
                            edgeVector[j].vertices[m].neighbourEdges.push_back(k);
                        }
                    }
            }
        }


        std::cout << "number of Lines = "<< edgeVector.size() << std::endl;

        //        for (unsigned k = 0; k < edgeVector.size(); k++)
        //        {
        //            std::cout << "edge: " << k << std::endl;
        //            for (unsigned j = 0; j < 2; j++)
        //            {
        //                std::cout << "  numberInself: " << edgeVector[k].isIt->numberInSelf() << std::endl;
        //                std::cout << "    vertice: " <<edgeVector[k].vertices[j].nodePointer->geometry().corner(0) << std::endl;
        //                std::cout << "        boundary: " << edgeVector[k].vertices[j].boundary << std::endl;
        //                for (unsigned n = 0; n < edgeVector[k].vertices[j].neighbourEdges.size(); n++)
        //                    std::cout << "            neighbour edges: " << edgeVector[k].vertices[j].neighbourEdges[n] << std::endl;
        //            }
        //        }

        //        std::cout << "Mapper size " << vertexToIndex.size() << std::endl;
        //        for (mapIterator mi = vertexToIndex.begin(); mi != vertexToIndex.end(); mi++)
        //        {
        //            std::cout << "Vertex " << mi->first->geometry().corner(0) << ", index " << mi->second << std::endl;
        //        }

    }

    //    // check subgrid
    //    std::cout << "Check subgrid with dim=" << dim << std::endl;
    //
    //    ElementIterator eIt    = subGrid.template lbegin<0>(0);
    //    ElementIterator eEndIt = subGrid.template lend<0>(0);
    //
    //    for (; eIt!=eEndIt; ++eIt) {
    //
    //        typename SubGrid<dim,G>::template Codim<0>::Entity::LevelIntersectionIterator nIt    = eIt->ilevelbegin();
    //        typename SubGrid<dim,G>::template Codim<0>::Entity::LevelIntersectionIterator nEndIt = eIt->ilevelend();
    //
    //        for (; nIt!=nEndIt; ++nIt) {
    //
    //            Dune::FieldVector<double,dim> v0 = nIt.intersectionGlobal()[0];
    //            Dune::FieldVector<double,dim> v1 = nIt.intersectionGlobal()[1];
    //
    //            std::cout << "v0: " << v0 << std::endl;
    //            std::cout << "v1: " << v1 << std::endl;
    //
    //            if (nIt.boundary())
    //                std::cout << "Boundary!\n";
    //            else
    //                std::cout << "Not a boundary!\n";
    //
    //        }
    //
    //    }


}

int main(int argc, char** argv)
{
    try{
        // set up the grid
        const int dim = 2;
        typedef double NumberType;

        //    typedef Dune::ALUSimplexGrid<dim,dim> GridType;
        typedef Dune::UGGrid<dim> GridType;

        // create grid pointer, GridType is defined by gridtype.hh
        Dune::GridPtr<GridType> gridPtr( "grids/rectangle_inner_manual.dgf" );

        // grid reference
        GridType& grid = *gridPtr;

        typedef GridType::Codim<dim>::EntityPointer HostVertexPointer;
        typedef GridType::Codim<0>::LeafIntersectionIterator HostIntersectionIterator;

        typedef Edge<GridType> Edges;

        typedef std::map<HostVertexPointer, int, VertexCompare>  MapperVertex;
        MapperVertex vertexToIndex;

        std::vector<Edges> edgeVector;

        isolate(grid, gridPtr, edgeVector, vertexToIndex);

        std::cout << "number of Lines = "<< edgeVector.size() << std::endl;

        for (unsigned k = 0; k < edgeVector.size(); k++)
        {
            std::cout << "edge: " << k << std::endl;
            for (unsigned j = 0; j < 2; j++)
            {
                std::cout << "  numberInself: " << edgeVector[k].isIt.numberInSelf() << std::endl;
                std::cout << "    vertice: " <<edgeVector[k].vertices[j].nodePointer->geometry().corner(0) << std::endl;
                std::cout << "        boundary: " << edgeVector[k].vertices[j].boundary << std::endl;
                for (unsigned n = 0; n < edgeVector[k].vertices[j].neighbourEdges.size(); n++)
                    std::cout << "            neighbour edges: " << edgeVector[k].vertices[j].neighbourEdges[n] << std::endl;
            }
        }

        typedef std::map<HostVertexPointer, int>::const_iterator mapIterator;

        std::cout << "Mapper size " << vertexToIndex.size() << std::endl;
        for (mapIterator mi = vertexToIndex.begin(); mi != vertexToIndex.end(); mi++)
        {
            std::cout << "Vertex " << mi->first->geometry().corner(0) << ", index " << mi->second << std::endl;
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
