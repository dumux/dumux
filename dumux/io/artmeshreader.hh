/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Artmesh reader: Reads Artmesh files (ASCII) and constructs a UG grid
 *  for modeling lower dimensional discrete fracture-matrix problems.
 */

#ifndef DUMUX_IO_ARTMESH_READER_HH
#define DUMUX_IO_ARTMESH_READER_HH

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/uggrid.hh>

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

// plot the vertices, edges, elements
#define PLOT 0
// TODO: this is a verbosity level

namespace Dumux
{
#ifdef HAVE_UG
typedef Dune::UGGrid<2> GridType; // Currently implemented for 2D fractured systems
#else
typedef Dune::YaspGrid<2> GridType;
#endif

typedef Dune::FieldVector<double, 3> Coordinates;
typedef std::vector<Coordinates> VerticesVector;
typedef Dune::FieldVector<int, 3> EdgePoints;
typedef std::vector<EdgePoints> EdgesVector;
typedef Dune::FieldVector<int, 4> Faces;
typedef std::vector<Faces> FacesVector;

typedef double NumberType; // TODO: from parameter system?

/*
The number of
    boundary vertices,
    boundary edges (connecting two boundary vertices),
    boundary faces (connecting n boundary edges) and
    boundary elements (should be 0)
are read from the boundarty description file and stored here.
*/
int bVertices, bEdges, bFaces, bElements;

// Store the vertex coordinates, edge infos, face infos.
Dune::FieldVector<double, 2> vertexPosition;

template<int dim>
struct P1Layout
{
    bool contains (Dune::GeometryType gt)
    {
        return (gt.dim() == 0);
    }
};

/*!
 * \brief Reads in Artmesh files.
 */
template <class GridType>
class ArtReader {

public:
    /*!
     * \brief Reads the file formats obtained with the Artmesh generator
     * 
     * \param
     */
    void read_art_file(const char *artFileName)
    {
        std::cout << "Opening " << artFileName << std::endl;
        std::ifstream inFile(artFileName);

        //intialize
        Coordinates coordinate_;
        EdgePoints edge_points_;
        Faces faces_;
        std::string jump;

        while (inFile >> jump)
        {
            //do nothing with the first lines in the code
            //start reading information only after getting to %Points
            if (jump == "VertexNumber:")
            {
                inFile >> jump;
                double dummy = atof(jump.c_str());
                vertexNumber = dummy;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "EdgeNumber:")
            {
                inFile >> jump;
                double dummy = atof(jump.c_str());
                edgeNumber = dummy;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "FaceNumber:")
            {
                inFile >> jump;
                double dummy = atof(jump.c_str());
                faceNumber = dummy;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "ElementNumber:")
            {
                inFile >> jump;
                double dummy = atof(jump.c_str());
                elementNumber = dummy;
                break;
            }
        }
        //// finished reading the header information: total number of verts, etc..
        ///////////////////////////////////////////////////////////////////////////
        while (inFile >> jump)
        {
            //jump over the lines until the ones with the vertices "% Vertices: x y z"
            if (jump == "Vertices:")
            {
                break;
            }
        }
        while (inFile >> jump)
        {
            //jskip words until "z" from the line "% Vertices: x y z"
            if (jump == "z")
            {
                std::cout << "Start reading the vertices" << std::endl;
                break;
            }
        }

        while (inFile >> jump)
        {
            if (jump == "$")
            {
                std::cout << "Finished reading the vertices" << std::endl;
                break;
            }
            double dummy = atof(jump.c_str());
            coordinate_[0] = dummy;
            for (int k=1; k<3; k++)
            {
                inFile >> jump;
                dummy = atof(jump.c_str());
                coordinate_[k] = dummy;
            }
            vertices_.push_back(coordinate_);
        }
    /////Reading Edges
        while (inFile >> jump)
        {
            //jump over the characters until the ones with the edges
            if (jump == "Points):")
            {
                std::cout << std::endl << "Start reading the edges" << std::endl;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "$")
            {
                std::cout << "Finished reading the edges" << std::endl;
                break;
            }
            double dummy = atof(jump.c_str());
            edge_points_[0] = dummy;
            for (int k=1; k<3; k++)
            {
                inFile >> jump;
                dummy = atof(jump.c_str());
                edge_points_[k] = dummy;
            }
            edges_vector_.push_back(edge_points_);
        }

    /////Reading Faces
        while (inFile >> jump)
        {
            //jump over the characters until the ones with the edges
            if (jump == "Edges):")
            {
                std::cout << std::endl << "Start reading the elements" << std::endl;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "$")
            {
                std::cout << "Finished reading the elements" << std::endl;
                break;
            }
            double dummy = atof(jump.c_str());
            faces_[0] = dummy;
            for (int k=1; k<4; k++)
            {
                inFile >> jump;
                dummy = atof(jump.c_str());
                faces_[k] = dummy;
            }
            faces_vector_.push_back(faces_);
        }
    }
    void outputARTtoScreen()
    {
    ////////OUTPUT for verification
    //////////////////////////////////////////////////////////////////
            std::cout << std::endl << "printing VERTICES" << std::endl;
            for (int i = 0; i < vertices_.size(); i++)
            {
                for (int j=0; j < 3; j++)
                std::cout << vertices_[i][j]<<"\t";
                std::cout << std::endl;
            }

            std::cout << std::endl << "printing EDGES" << std::endl;
            for (int i = 0; i < edges_vector_.size(); i++)
            {
                for (int j=0; j < 3; j++)
                std::cout<<edges_vector_[i][j]<<"\t";
                std::cout<<std::endl;
            }
            ////Printing Faces
            std::cout << std::endl << "printing FACES" << std::endl;
            for (int i = 0; i < faces_vector_.size(); i++)
            {
                for (int j=0; j < 4; j++)
                {
                    std::cout<<faces_vector_[i][j]<<" ";
                }
                std::cout<<std::endl;
            }

            std::cout << std::endl << "Total number of vertices "<<vertexNumber<<std::endl;
            std::cout << "Total number of edges: "<<edgeNumber<<std::endl;
            std::cout << "Total number of faces: "<<faceNumber<<std::endl;
            std::cout << "Total number of elements: "<<elementNumber<<std::endl;

            if (vertices_.size() != vertexNumber )
            {
                std::cout << "The total given number of vertices: "<<vertexNumber
                        <<" is not the same with the read number of entries: "
                        <<vertices_.size()<<"" << std::endl;
            }
            if (edges_vector_.size() != edgeNumber )
            {
                std::cout << "The total given number of edges: "<<edgeNumber
                        <<" is not the same with the read number of entries: "
                        <<edges_vector_.size()<<"" << std::endl;
            }
            if (faces_vector_.size() != faceNumber )
            {
                std::cout << "The total given number of faces: "<<faceNumber
                        <<" is not the same with the read number of entries: "
                        <<faces_vector_.size()<<"" << std::endl;
            }
    }

    /*!
     * \brief Create the UG grid from the read in data
     */
    GridType* createGrid()
    {
        // set up the grid factory
        Dune::FieldVector<double,2> position;
        Dune::GridFactory<GridType> factory;
        bVertices = vertexNumber;
        bEdges = edgeNumber;
        bFaces = faceNumber;
        VerticesVector verts(bVertices);
        EdgesVector edges(bEdges);
        FacesVector elems(bFaces);

        // Plot the vertices
    #if PLOT
        std::cout << "*================*"<<std::endl;
        std::cout << "* Vertices"<<std::endl;
        std::cout << "*================*"<<std::endl;
    #endif
        for (int k=0; k<bVertices; k++)
        {
            verts[k][0] = vertices_[k][0];
            verts[k][1] = vertices_[k][1];
            verts[k][2] = vertices_[k][2];

            //Printing the vertices vector
    #if PLOT
            std::cout << verts[k][0]<<"\t\t";
            std::cout << verts[k][1]<<"\t\t";
            std::cout << verts[k][2]<<std::endl;
    #endif
            for (int i=0; i<2; i++)
            {
                position[i]=verts[k][i];
            }
            factory.insertVertex(position);
        }
    #if PLOT
        std::cout << "*================*"<<std::endl;
        std::cout << "* Edges"<<std::endl;
        std::cout << "*================*"<<std::endl;
    #endif
        for (int k=0; k<bEdges; k++)
        {
            edges[k][0]=edges_vector_[k][0];
            edges[k][1]=edges_vector_[k][1];
            edges[k][2]=edges_vector_[k][2];
    #if PLOT
            //Printing the Edge vector
            std::cout<<edges[k][0]<<"\t\t";
            std::cout<<edges[k][1]<<"\t";
            std::cout<<edges[k][2]<<std::endl;
    #endif
        }

        // Plot the Elements (Faces)
    #if PLOT
        std::cout << "*================*"<<std::endl;
        std::cout << "* Faces"<<std::endl;
        std::cout << "*================*"<<std::endl;
    #endif

        for (int k=0; k<faceNumber; k++)
        {
            elems[k][0]=faces_vector_[k][1];
            elems[k][1]=faces_vector_[k][2];
            elems[k][2]=faces_vector_[k][3];
    #if PLOT
            std::cout<<elems[k][0]<<"\t";
            std::cout<<elems[k][1]<<"\t";
            std::cout<<elems[k][2]<<std::endl;
    #endif
        }

        //***********************************************************************//
        //Create the Elements in Dune::GridFactory //
        //***********************************************************************//

        for (int i=0; i<bFaces; i++)
        {
            std::vector<unsigned int> nodeIdx(3);
            Dune::FieldVector<double,3> point(0);
    #if PLOT
            std::cout << "====================================="<<std::endl;
            std::cout << "globalElemIdx "<<i<<std::endl;
            std::cout << "====================================="<<std::endl;
    #endif
            int edgeIdx = 0;
            //first node of the element - from first edge Node 1
            nodeIdx[0] = edges[elems[i][edgeIdx]][1];
            //second node of the element- from first edge Node 2
            nodeIdx[1] = edges[elems[i][edgeIdx]][2];
            //third node of the element - from the second edge
            nodeIdx[2] = edges[elems[i][edgeIdx+1]][1];
            // if the nodes of the edges are identical swap
            if (nodeIdx[1] == nodeIdx[2] || nodeIdx[0] == nodeIdx[2])
            {
                nodeIdx[2] = edges[elems[i][edgeIdx+1]][2];
            }

            /* Check if the order of the nodes is trigonometric
            by computing the cross product
            If not then the two nodes are switched among each other*/
            Dune::FieldVector<double, 2> v(0);
            Dune::FieldVector<double, 2> w(0);
            double cross1;
            v[0] = verts[nodeIdx[0]][0] - verts[nodeIdx[1]][0];
            v[1] = verts[nodeIdx[0]][1] - verts[nodeIdx[1]][1];
            w[0] = verts[nodeIdx[0]][0] - verts[nodeIdx[2]][0];
            w[1] = verts[nodeIdx[0]][1] - verts[nodeIdx[2]][1];
            cross1 = v[0]*w[1]-v[1]*w[0];
            //If the cross product is negative switch the order of the vertices
            if (cross1 < 0)
            {
                nodeIdx[0] = edges[elems[i][edgeIdx]][2]; //node 0 is node 1
                nodeIdx[1] = edges[elems[i][edgeIdx]][1]; //node 1 is node 0
            }
            v[0] = verts[nodeIdx[0]][0] - verts[nodeIdx[1]][0];
            v[1] = verts[nodeIdx[0]][1] - verts[nodeIdx[1]][1];
            w[0] = verts[nodeIdx[0]][0] - verts[nodeIdx[2]][0];
            w[1] = verts[nodeIdx[0]][1] - verts[nodeIdx[2]][1];

            factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),
                    nodeIdx);
    #if PLOT
            std::cout << "edges of the element "<< elems[i] << std::endl;
            std::cout << "nodes of the element " << nodeIdx[0]
                      << ", " << nodeIdx[1] << ", "
                      << nodeIdx[2] << std::endl;
            std::cout << "1st "<<nodeIdx[0]<<"\t"<<"2nd "<<nodeIdx[1]
                      << "\t"<<"3rd "<<nodeIdx[2]
                      << std::endl;
    #endif
            if (nodeIdx[0] == nodeIdx[1]||nodeIdx[1] == nodeIdx[2]
                                        ||nodeIdx[0] == nodeIdx[2] )
            {
                std::cout << "Error. The node index is identical in the element"
                          << std::endl;
            }

        }

        Valgrind::CheckDefined(verts);

        return factory.createGrid();
    }

public:
    VerticesVector vertices_;
    EdgesVector edges_vector_;
    FacesVector faces_vector_ ;
    int vertexNumber;
    int edgeNumber;
    int faceNumber; //in 2D
    int elementNumber; //in 3D
};

/*!
 * \brief Maps the fractures from the Artmesh to the UG grid ones
 */
template<class GridType>
class FractureMapper
{
public:
    // mapper: one data element in every entity
    template<int dim>
    struct FaceLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim-1;
        }
    };
    template<int dim>
    struct VertexLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == 0;
        }
    };
    typedef typename GridType::LeafGridView GV;
    typedef typename GridType::ctype DT;
    enum {dim=GridType::dimension};
    typedef typename GV::IndexSet IS;
    typedef typename GridType::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<dim>::Iterator VertexLeafIterator;
    typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;
    typedef typename GridType::template Codim<0>::EntityPointer EEntityPointer;
    typedef typename GridType::Traits::GlobalIdSet IDS;
    typedef typename IDS::IdType IdType;
    typedef std::set<IdType> GIDSet;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,FaceLayout> FaceMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,VertexLayout> VertexMapper;
    typedef Dumux::ArtReader<GridType> ArtReader;

public:
    /*!
     * \brief Constructor
     * 
     * \param grid The grid
     * \param geomArt Artmesh reader
     */
    FractureMapper (const GridType& grid, ArtReader geomArt)
        : grid_(grid),
        facemapper(grid.leafView()),
        vertexmapper(grid.leafView()),
        artgeom_(geomArt)
    {}
    
    /*!
     * \brief Constructor
     * 
     * \param grid The grid
     * \param geomArt Artmesh reader
     */
    FractureMapper (const GridType& grid)
        : grid_(grid),
        facemapper(grid.leafView()),
        vertexmapper(grid.leafView())
    {}

    /*!
     * 
     */
    void fractureMapper()
    {
        //call the new_read art
        int nVertices = artgeom_.vertexNumber;
        int nEdges = artgeom_.edgeNumber;
        //The vertexes which are located on fractures
        isDuneFractureVertex_.resize(nVertices);
        std::fill(isDuneFractureVertex_.begin(), isDuneFractureVertex_.end(), false);

        //The edge which are fractures
        isDuneFractureEdge_.resize(nEdges);
        fractureEdgesIdx_.resize(nEdges);
        std::fill(isDuneFractureEdge_.begin(), isDuneFractureEdge_.end(), false);

        GV leafView = grid_.leafView();

        ElementLeafIterator eendit = leafView. template end<0>();
        for (ElementLeafIterator it = leafView .template begin<0>(); it != eendit; ++it)
        {
             Dune::GeometryType gt = it->geometry().type();
             const typename Dune::GenericReferenceElementContainer<DT,dim>::value_type&
                 refElem = Dune::GenericReferenceElements<DT,dim>::general(gt);

              // Loop over element faces
              for (int i = 0; i < refElem.size(1); i++)
              {
                  int indexFace = facemapper.map(*it, i, 1);
                  /*
                  * it maps the local element vertices "localV1Idx" -> indexVertex1
                  * then it gets the coordinates of the nodes in the ART file and
                  * by comparing them with the ones in the DUNE grid maps them too.
                  */
                  int localV1Idx = refElem.subEntity(i, 1, 0, dim);
                  int localV2Idx = refElem.subEntity(i, 1, 1, dim);
                  int indexVertex1 = vertexmapper.map(*it, localV1Idx, dim);
                  int indexVertex2 = vertexmapper.map(*it, localV2Idx, dim);
                  Dune::FieldVector<DT, dim> nodeART_from;
                  Dune::FieldVector<DT, dim> nodeART_to;
                  Dune::FieldVector<DT, dim> nodeDune_from;
                  Dune::FieldVector<DT, dim> nodeDune_to;

                  nodeDune_from = it-> geometry(). corner (localV1Idx);
                  nodeDune_to = it-> geometry(). corner (localV2Idx);

                  for (int j=0; j < nEdges; j++)
                  {
                      nodeART_from[0] = artgeom_.vertices_[artgeom_.edges_vector_[j][1]][0];
                      nodeART_from[1] = artgeom_.vertices_[artgeom_.edges_vector_[j][1]][1];
                      nodeART_to[0] = artgeom_.vertices_[artgeom_.edges_vector_[j][2]][0];
                      nodeART_to[1] = artgeom_.vertices_[artgeom_.edges_vector_[j][2]][1];
            
                      if ((nodeART_from == nodeDune_from && nodeART_to == nodeDune_to)
                          || (nodeART_from == nodeDune_to && nodeART_to == nodeDune_from))
                      {
#if PLOT
                        std::cout << " ART edge index is "<< j;
//                            std::cout << " Node ART_from: "
//                            << nodeART_from<<std::endl;
                        std::cout << " Dune edge is " << indexFace
                        <<std::endl;
#endif
                        /* assigns a value 1 for the edges
                        * which are fractures */
                        if (artgeom_.edges_vector_[j][0] < 0)
                        {
                            isDuneFractureEdge_[indexFace] = true;
                            fractureEdgesIdx_[indexFace]   = artgeom_.edges_vector_[j][0];
                            isDuneFractureVertex_[indexVertex1]=true;
                            isDuneFractureVertex_[indexVertex2]=true;

#if PLOT
                            std::cout << "isDuneFractureEdge_ "<<
                            isDuneFractureEdge_[indexFace]<<"\t";
                            std::cout << "vertex1 "<<indexVertex1<<"\t";
                            std::cout << "vertex2 "<<indexVertex2<<"" << std::endl;
#endif
                        }
                    }
                }
            }
         }

#if PLOT
        int i=0;
        for (VertexLeafIterator vIt=leafView. template begin<dim> (); 
             vIt != leafView. template end<dim> (); ++vIt)
        {
            Dune::GeometryType gt = vIt->type();
            std::cout << "visiting " << gt
                      << " at " << vIt->geometry().corner(0)
                      << "\t" << isDuneFractureVertex_[i]
                      << std::endl;
            i++;
        }
#endif

    }

private:
   const GridType& grid_;
   FaceMapper facemapper;
   VertexMapper vertexmapper;
   ArtReader artgeom_;
public:
   std::vector<bool> isDuneFractureVertex_;
   std::vector<bool> isDuneFractureEdge_;
   std::vector<int>  fractureEdgesIdx_;
};

} // end namespace

#endif // DUMUX_IO_ARTMESH_READER_HH
