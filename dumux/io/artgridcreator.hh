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

#ifndef DUMUX_ARTGRIDCREATOR_HH
#define DUMUX_ARTGRIDCREATOR_HH

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/common/math.hh>
#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(GridCreator);
}

/*!
 * \brief Reads in Artmesh files.
 */
template <class TypeTag>
class ArtGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid)  Grid;
    typedef Grid* GridPointer;
    typedef Dune::FieldVector<double, 3> Coordinates;
    typedef std::vector<Coordinates> VerticesVector;
    typedef Dune::FieldVector<int, 3> EdgePoints;
    typedef std::vector<EdgePoints> EdgesVector;
    typedef Dune::FieldVector<int, 4> Faces;
    typedef std::vector<Faces> FacesVector;

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        const std::string artFileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Grid, File);

        bool verbose = false;
        try { verbose = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Grid, Verbosity);}
        catch (Dumux::ParameterException &e) { }

        if (verbose)
            std::cout << "Opening " << artFileName << std::endl;
        std::ifstream inFile(artFileName.c_str());

        std::string jump;
        while (inFile >> jump)
        {
            //do nothing with the first lines in the code
            //start reading information only after getting to %Points
            if (jump == "VertexNumber:")
            {
                inFile >> jump;
                double dummy = atof(jump.c_str());
                vertexNumber_ = dummy;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "EdgeNumber:")
            {
                inFile >> jump;
                double dummy = atof(jump.c_str());
                edgeNumber_ = dummy;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "FaceNumber:")
            {
                inFile >> jump;
                double dummy = atof(jump.c_str());
                faceNumber_ = dummy;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "ElementNumber:")
            {
                inFile >> jump;
                double dummy = atof(jump.c_str());
                elementNumber_ = dummy;
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
                if (verbose)
                    std::cout << "Start reading the vertices" << std::endl;
                break;
            }
        }

        while (inFile >> jump)
        {
            if (jump == "$")
            {
                if (verbose)
                    std::cout << "Finished reading the vertices" << std::endl;
                break;
            }
            double dummy = atof(jump.c_str());
            Coordinates coordinate;
            coordinate[0] = dummy;
            for (int k = 1; k < 3; k++)
            {
                inFile >> jump;
                dummy = atof(jump.c_str());
                coordinate[k] = dummy;
            }
            vertices_.push_back(coordinate);
        }
    /////Reading Edges
        while (inFile >> jump)
        {
            //jump over the characters until the ones with the edges
            if (jump == "Points):")
            {
                if (verbose)
                    std::cout << std::endl << "Start reading the edges" << std::endl;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "$")
            {
                if (verbose)
                    std::cout << "Finished reading the edges" << std::endl;
                break;
            }
            double dummy = atof(jump.c_str());
            EdgePoints edgePoints;
            edgePoints[0] = dummy;
            for (int k = 1; k < 3; k++)
            {
                inFile >> jump;
                dummy = atof(jump.c_str());
                edgePoints[k] = dummy;
            }
            edges_.push_back(edgePoints);
        }

    /////Reading Faces
        while (inFile >> jump)
        {
            //jump over the characters until the ones with the edges
            if (jump == "Edges):")
            {
                if (verbose)
                    std::cout << std::endl << "Start reading the elements" << std::endl;
                break;
            }
        }
        while (inFile >> jump)
        {
            if (jump == "$")
            {
                if (verbose)
                    std::cout << "Finished reading the elements" << std::endl;
                break;
            }
            double dummy = atof(jump.c_str());
            Faces facePoints;
            facePoints[0] = dummy;
            for (int k = 1; k < 4; k++)
            {
                inFile >> jump;
                dummy = atof(jump.c_str());
                facePoints[k] = dummy;
            }
            faces_.push_back(facePoints);
        }

                // set up the grid factory
        Dune::FieldVector<double,2> position;
        Dune::GridFactory<Grid> factory;

        // Plot the vertices
        if (verbose)
        {
            std::cout << "*================*" << std::endl;
            std::cout << "* Vertices" << std::endl;
            std::cout << "*================*" << std::endl;
        }

        for (int k = 0; k < vertexNumber_; k++)
        {

            if (verbose)
            {
                // Printing the vertices vector
                std::cout << vertices_[k][0] << "\t\t";
                std::cout << vertices_[k][1] << "\t\t";
                std::cout << vertices_[k][2] << std::endl;
            }

            for (int i = 0; i < 2; i++)
            {
                position[i] = vertices_[k][i];
            }
            factory.insertVertex(position);
        }

        if (verbose)
        {
            std::cout << "*================*" << std::endl;
            std::cout << "* Edges" << std::endl;
            std::cout << "*================*" << std::endl;
            for (int k = 0; k < edgeNumber_; k++)
            {
                // Printing the Edge vector
                std::cout << edges_[k][0] << "\t\t";
                std::cout << edges_[k][1] << "\t";
                std::cout << edges_[k][2] << std::endl;
            }

            // Plot the Elements (Faces)
            std::cout << "*================*" << std::endl;
            std::cout << "* Faces" << std::endl;
            std::cout << "*================*" << std::endl;
            for (int k = 0; k < faceNumber_; k++)
            {
                std::cout << faces_[k][1] << "\t";
                std::cout << faces_[k][2] << "\t";
                std::cout << faces_[k][3] << std::endl;
            }
        }

        //***********************************************************************//
        //Create the Elements in Dune::GridFactory //
        //***********************************************************************//

        for (int i=0; i<faceNumber_; i++)
        {
            std::vector<unsigned int> vIdxGlobal(3);
            Dune::FieldVector<double,3> point(0);
            if (verbose)
            {
                std::cout << "=====================================" << std::endl;
                std::cout << "eIdxGlobal " << i << std::endl;
                std::cout << "=====================================" << std::endl;
            }

            int edgeIdx = 0;
            //first node of the element - from first edge Node 1
            vIdxGlobal[0] = edges_[faces_[i][edgeIdx+1]][1];
            //second node of the element- from first edge Node 2
            vIdxGlobal[1] = edges_[faces_[i][edgeIdx+1]][2];
            //third node of the element - from the second edge
            vIdxGlobal[2] = edges_[faces_[i][edgeIdx+2]][1];
            // if the nodes of the edges are identical swap
            if (vIdxGlobal[1] == vIdxGlobal[2] || vIdxGlobal[0] == vIdxGlobal[2])
            {
                vIdxGlobal[2] = edges_[faces_[i][edgeIdx+2]][2];
            }

            /* Check if the order of the nodes is trigonometric
            by computing the cross product
            If not then the two nodes are switched among each other*/
            Dune::FieldVector<double, 2> v(0);
            Dune::FieldVector<double, 2> w(0);
            double cross1;
            v[0] = vertices_[vIdxGlobal[0]][0] - vertices_[vIdxGlobal[1]][0];
            v[1] = vertices_[vIdxGlobal[0]][1] - vertices_[vIdxGlobal[1]][1];
            w[0] = vertices_[vIdxGlobal[0]][0] - vertices_[vIdxGlobal[2]][0];
            w[1] = vertices_[vIdxGlobal[0]][1] - vertices_[vIdxGlobal[2]][1];
            cross1 = v[0]*w[1]-v[1]*w[0];
            //If the cross product is negative switch the order of the vertices
            if (cross1 < 0)
            {
                vIdxGlobal[0] = edges_[faces_[i][edgeIdx+1]][2]; //node 0 is node 1
                vIdxGlobal[1] = edges_[faces_[i][edgeIdx+1]][1]; //node 1 is node 0
            }
            v[0] = vertices_[vIdxGlobal[0]][0] - vertices_[vIdxGlobal[1]][0];
            v[1] = vertices_[vIdxGlobal[0]][1] - vertices_[vIdxGlobal[1]][1];
            w[0] = vertices_[vIdxGlobal[0]][0] - vertices_[vIdxGlobal[2]][0];
            w[1] = vertices_[vIdxGlobal[0]][1] - vertices_[vIdxGlobal[2]][1];

            factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),
                    vIdxGlobal);

            if (verbose)
            {
                std::cout << "edges of the element "<< faces_[i] << std::endl;
                std::cout << "nodes of the element " << vIdxGlobal[0]
                          << ", " << vIdxGlobal[1] << ", "
                          << vIdxGlobal[2] << std::endl;
                std::cout << "1st " << vIdxGlobal[0]
                          << "\t" << "2nd " << vIdxGlobal[1]
                          << "\t" << "3rd " << vIdxGlobal[2]
                          << std::endl;
            }

            if (vIdxGlobal[0] == vIdxGlobal[1] || vIdxGlobal[1] == vIdxGlobal[2]
                                         || vIdxGlobal[0] == vIdxGlobal[2])
            {
                DUNE_THROW(Dune::InvalidStateException, "Two vertex indices of this element are identical");
            }
        }

        gridPtr() = factory.createGrid();
    }

    static void outputARTtoScreen()
    {
    ////////OUTPUT for verification
    //////////////////////////////////////////////////////////////////
            std::cout << std::endl << "printing VERTICES" << std::endl;
            for (int i = 0; i < vertices_.size(); i++)
            {
                for (int j = 0; j < 3; j++)
                    std::cout << vertices_[i][j] << "\t";
                std::cout << std::endl;
            }

            std::cout << std::endl << "printing EDGES" << std::endl;
            for (int i = 0; i < edges_.size(); i++)
            {
                for (int j = 0; j < 3; j++)
                    std::cout << edges_[i][j] << "\t";
                std::cout << std::endl;
            }
            ////Printing Faces
            std::cout << std::endl << "printing FACES" << std::endl;
            for (int i = 0; i < faces_.size(); i++)
            {
                for (int j = 0; j < 4; j++)
                    std::cout << faces_[i][j] << " ";
                std::cout << std::endl;
            }

            std::cout << std::endl << "Total number of vertices " << vertexNumber_ << std::endl;
            std::cout << "Total number of edges: " << edgeNumber_ << std::endl;
            std::cout << "Total number of faces: " << faceNumber_ << std::endl;
            std::cout << "Total number of elements: " << elementNumber_ << std::endl;

            if (vertices_.size() != vertexNumber_ )
            {
                std::cout << "The total given number of vertices: " << vertexNumber_
                        << " is not the same with the read number of entries: "
                        << vertices_.size() << "" << std::endl;
            }
            if (edges_.size() != edgeNumber_ )
            {
                std::cout << "The total given number of edges: " << edgeNumber_
                        << " is not the same with the read number of entries: "
                        << edges_.size() << "" << std::endl;
            }
            if (faces_.size() != faceNumber_ )
            {
                std::cout << "The total given number of faces: " << faceNumber_
                        << " is not the same with the read number of entries: "
                        << faces_.size() << "" << std::endl;
            }
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *gridPtr();
    }

    /*!
     * \brief Returns a reference to the shared pointer to the grid.
     */
    static GridPointer &gridPtr()
    {
        static GridPointer grid;
        return grid;
    }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    static void loadBalance()
    {
        gridPtr()->loadBalance();
    }

    static int vertexNumber()
    {
        return vertexNumber_;
    }

    static int edgeNumber()
    {
        return edgeNumber_;
    }

    static int faceNumber()
    {
        return faceNumber_;
    }

    static VerticesVector& vertices()
    {
        return vertices_;
    }

    static EdgesVector& edges()
    {
        return edges_;
    }

    static FacesVector& faces()
    {
        return faces_;
    }

private:
    static VerticesVector vertices_;
    static EdgesVector edges_;
    static FacesVector faces_ ;
    static int vertexNumber_;
    static int edgeNumber_;
    static int faceNumber_; //in 2D
    static int elementNumber_; //in 3D
};

template <class TypeTag>
typename Dumux::ArtGridCreator<TypeTag>::VerticesVector ArtGridCreator<TypeTag>::vertices_;
template <class TypeTag>
typename Dumux::ArtGridCreator<TypeTag>::EdgesVector ArtGridCreator<TypeTag>::edges_;
template <class TypeTag>
typename Dumux::ArtGridCreator<TypeTag>::FacesVector ArtGridCreator<TypeTag>::faces_;
template <class TypeTag>
int ArtGridCreator<TypeTag>::vertexNumber_;
template <class TypeTag>
int ArtGridCreator<TypeTag>::edgeNumber_;
template <class TypeTag>
int ArtGridCreator<TypeTag>::faceNumber_;
template <class TypeTag>
int ArtGridCreator<TypeTag>::elementNumber_;

/*!
 * \brief Maps the fractures from the Artmesh to the UG grid ones
 */
template<class TypeTag>
class FractureMapper
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView)  GridView;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator)  GridCreator;
public:
    // mapper: one data element in every entity
    template<int dim>
    struct FaceLayout
    {
        bool contains (Dune::GeometryType geomType)
        {
            return geomType.dim() == dim-1;
        }
    };
    typedef typename GridView::ctype DT;
    enum {dim = GridView::dimension};
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, FaceLayout> FaceMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGVertexLayout> VertexMapper;

public:
    /*!
     * \brief Constructor
     *
     * \param gridView The grid view
     */
    FractureMapper (const GridView& gridView)
    : gridView_(gridView),
      faceMapper_(gridView),
      vertexMapper_(gridView)
    {}

    /*! DOC ME!! and RENAME ME if you know what exactly I'm doing!
     *
     */
    void map()
    {
        //call the new_read art
        int numVertices = GridCreator::vertexNumber();
        int nEdges = GridCreator::edgeNumber();
        //The vertexes which are located on fractures
        isDuneFractureVertex_.resize(numVertices);
        std::fill(isDuneFractureVertex_.begin(), isDuneFractureVertex_.end(), false);

        //The edge which are fractures
        isDuneFractureEdge_.resize(nEdges);
        fractureEdgesIdx_.resize(nEdges);
        std::fill(isDuneFractureEdge_.begin(), isDuneFractureEdge_.end(), false);

        for (const auto& element : elements(gridView_))
        {
             Dune::GeometryType geomType = element.geometry().type();

             const auto& refElement = Dune::ReferenceElements<DT,dim>::general(geomType);

              // Loop over element faces
              for (int i = 0; i < refElement.size(1); i++)
              {
                  int indexFace = faceMapper_.subIndex(element, i, 1);

                  /*
                  * it maps the local element vertices "localV1Idx" -> indexVertex1
                  * then it gets the coordinates of the nodes in the ART file and
                  * by comparing them with the ones in the DUNE grid maps them too.
                  */
                  int localV1Idx = refElement.subEntity(i, 1, 0, dim);
                  int localV2Idx = refElement.subEntity(i, 1, 1, dim);
                  int indexVertex1 = vertexMapper_.subIndex(element, localV1Idx, dim);
                  int indexVertex2 = vertexMapper_.subIndex(element, localV2Idx, dim);

                  Dune::FieldVector<DT, dim> nodeART_from;
                  Dune::FieldVector<DT, dim> nodeART_to;
                  Dune::FieldVector<DT, dim> nodeDune_from;
                  Dune::FieldVector<DT, dim> nodeDune_to;

                  nodeDune_from = element.geometry().corner(localV1Idx);
                  nodeDune_to = element.geometry().corner(localV2Idx);

                  for (int j=0; j < nEdges; j++)
                  {
                      nodeART_from[0] = GridCreator::vertices()[GridCreator::edges()[j][1]][0];
                      nodeART_from[1] = GridCreator::vertices()[GridCreator::edges()[j][1]][1];
                      nodeART_to[0] = GridCreator::vertices()[GridCreator::edges()[j][2]][0];
                      nodeART_to[1] = GridCreator::vertices()[GridCreator::edges()[j][2]][1];

                      if ((nodeART_from == nodeDune_from && nodeART_to == nodeDune_to)
                          || (nodeART_from == nodeDune_to && nodeART_to == nodeDune_from))
                      {
                        /* assigns a value 1 for the edges
                        * which are fractures */
                        if (GridCreator::edges()[j][0] < 0)
                        {
                            isDuneFractureEdge_[indexFace] = true;
                            fractureEdgesIdx_[indexFace]   = GridCreator::edges()[j][0];
                            isDuneFractureVertex_[indexVertex1]=true;
                            isDuneFractureVertex_[indexVertex2]=true;
                        }
                    }
                }
            }
        }
    }

    bool isDuneFractureVertex(unsigned int i) const
    {
        return isDuneFractureVertex_[i];
    }

    bool isDuneFractureEdge(unsigned int i) const
    {
        return isDuneFractureEdge_[i];
    }

    int fractureEdgesIdx(unsigned int i) const
    {
        return fractureEdgesIdx_[i];
    }

private:
    const GridView gridView_;
    FaceMapper faceMapper_;
    VertexMapper vertexMapper_;
    std::vector<bool> isDuneFractureVertex_;
    std::vector<bool> isDuneFractureEdge_;
    std::vector<int>  fractureEdgesIdx_;
};

} // end namespace

#endif // DUMUX_ARTGRIDCREATOR_HH
