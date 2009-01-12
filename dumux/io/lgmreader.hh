#ifndef DUNE_STARCD_READER_HH
#define DUNE_STARCD_READER_HH

#include <fstream>

namespace Dune {

    /** @ingroup StarCD
     *    \brief File reader for the Star-CD format.
     *
     *    Reads grid data described by the Star-CD format and fills
     *    an empty grid with the data. Uses the grid creation methods
     *    described in
     *    <a href="http://www.dune-project.org/doc/devel/ugreader.html">
     *    "How to Write a File Reader for UGGrid Objects"</a>.
     *
     *    Two files \a fileName.vrt and \a fileName.cel have to be present. The file fileName.vrt contains the
     *    coordinates of the nodes, each row having the format
     *    \code
     *    idx  x-coordinate  y-coordinate  z-coordinate
     *    \endcode
     *
     *    The file fileName.cel contains the data of the volume and possibly the boundary elements, each
     *    row having the format
     *    \code
     *    idx  node1  node2  ...  node8  material/boundaryId  flag1 flag2
     *    \endcode
     *    The flags flag1 and flag2 appear to be always identical to 1 in case of a volume element
     *    and to 4 in case of a boundary element. The element types simplex, pyramid, prism, and
     *    cube are supported:
     *
     *    - For cubes, the indices node1 ... node8 are what you expect them to be.
     *    - For simplices, node3 and node4 are identical, as well as node5 ... node8.
     *    - For pyramids, node5 ... node8 are identical.
     *    - For prisms, node3 and node4 are identical, as well as node7 and node8.
     *
     *    This reader only supports three-dimensional grids.
     *
     *    Currently no boundary element data is passed to \a grid.
     */
    template <class GridType>
    class StarCDReader {

    public:

        /** \brief Read grid from a Star-CD file
         *    \param grid The grid
         *    \param fileName The base file name of the Star-CD files
         *    \param verbose Tlag to set whether information should be printed
         */
        static GridType* read(const std::string& fileName, bool verbose = true)
        {
            // extract the grid dimension
            const int dim = GridType::dimension;

            // currently only dim = 3 is implemented
            if (dim != 3)
                DUNE_THROW(Dune::NotImplemented,
                           "Reading Star-CD format is not implemented for dimension " << dim);

            // set up the grid factory
            GridFactory<GridType> factory;

            // set the name of the vertex file
            std::string vertexFileName = fileName + ".vrt";

            // set the vertex input stream
            std::ifstream vertexFile(vertexFileName.c_str());
            if (!vertexFile)
                DUNE_THROW(Dune::IOError, "Could not open " << vertexFileName);

            // read the vertices
            int dummyIdx;
            int numberOfVertices = 0;
            while (vertexFile >> dummyIdx) {
                numberOfVertices++;

                Dune::FieldVector<double,dim> position;

                for (int k = 0; k < dim; k++)
                    vertexFile >> position[k];

                factory.insertVertex(position);
            }
            if (verbose)
                std::cout << numberOfVertices << " vertices read." << std::endl;

            // set the name of the element file
            std::string elementFileName = fileName + ".cel";

            // set the element input stream
            std::ifstream elementFile(elementFileName.c_str());
            if (!elementFile)
                DUNE_THROW(Dune::IOError, "Could not open " << elementFileName);

            // read the elements
            int numberOfElements = 0;
            int numberOfSimplices = 0;
            int numberOfPyramids = 0;
            int numberOfPrisms = 0;
            int numberOfCubes = 0;;
            int maxNumberOfVertices = (int)pow(2, dim);
            int isVolume = 1;
            while (elementFile >> dummyIdx) {
                std::vector<unsigned int> vertices(maxNumberOfVertices);
                for (int k = 0; k < maxNumberOfVertices; k++)
                    elementFile >> vertices[k];

                int boundaryId;
                elementFile >> boundaryId;

                int volumeOrSurface[2];
                elementFile >> volumeOrSurface[0] >> volumeOrSurface[1];

                if (volumeOrSurface[0] == isVolume) {
                    numberOfElements++;

                    if (vertices[2] == vertices[3]) { // simplex or prism
                        if (vertices[4] == vertices[5]) { // simplex
                            numberOfSimplices++;
                            std::vector<unsigned int> simplexVertices(4);
                            for (int k = 0; k < 3; k++)
                                simplexVertices[k] = vertices[k] - 1;
                            simplexVertices[3] = vertices[4] - 1;
                            factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim), simplexVertices);
                        }
                        else { // prism
                            numberOfPrisms++;
                            std::vector<unsigned int> prismVertices(6);
                            for (int k = 0; k < 3; k++)
                                prismVertices[k] = vertices[k] - 1;
                            for (int k = 3; k < 6; k++)
                                prismVertices[k] = vertices[k+1] - 1;
                            factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism,dim), prismVertices);
                        }
                    }
                    else { // cube or pyramid
                        if (vertices[4] == vertices[5]) { // pyramid
                            numberOfPyramids++;
                            std::vector<unsigned int> pyramidVertices(5);
                            for (int k = 0; k < 5; k++)
                                pyramidVertices[k] = vertices[k] - 1;
                            factory.insertElement(Dune::GeometryType(Dune::GeometryType::pyramid,dim), pyramidVertices);
                        }
                        else { // cube
                            numberOfCubes++;
                            std::vector<unsigned int> cubeVertices(8);
                            for (int k = 0; k < 8; k++)
                                cubeVertices[k] = vertices[k] - 1;
                            std::swap(cubeVertices[2], cubeVertices[3]);
                            std::swap(cubeVertices[6], cubeVertices[7]);
                            factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim), cubeVertices);
                        }
                    }
                }
            }
            if (verbose)
                std::cout << numberOfElements << " elements read: "
                          << numberOfSimplices << " simplices, " << numberOfPyramids << " pyramids, "
                          << numberOfPrisms << " prisms, " << numberOfCubes << " cubes." << std::endl;

            // finish off the construction of the grid object
            if (verbose)
                std::cout << "Starting createGrid() ... " << std::flush;

            return factory.createGrid();

        }

    };

}

#endif
