#ifndef DUNE_LGM_READER_HH
#define DUNE_LGM_READER_HH

#include <fstream>
#include <cstdlib>


using namespace std;


namespace Dune
{

/** @ingroup IO
 *    \brief File reader for the lgm  format.
 *
 *    Reads grid data described by the lgm format and fills
 *    an empty grid with the data. Uses the grid creation methods
 *    described in
 *    <a href="http://www.dune-project.org/doc/devel/ugreader.html">
 *    "How to Write a File Reader for UGGrid Objects"</a>.
 *
 *    Two files \a fileName.lgm and \a fileName.ng have to be present.
 *
 *
 *    The file fileName.lgm contains the data of the outer nodes (points), lines and surfaces, each
 *    file having the format
 *    \code
 *    #Domain-Info
 *        name = name
 *        problemname = problemname
 *        convex = 0
 *      #Unit-Info
 *        unit 1 MATERIAL1
 *      #Line-Info
 *        line l0: points: p0, p1, ...;
 *        line l1: ...
 *      #Surface-Info
 *        surface s0: left=1; right=0; points: node 1 node2 ...; lines: l0 l1 ...; triangles: node1 node2 node3; node4 node5 node6;...;
 *        surface s1 ...: left=1; right=0; points:
 *        #Point-Info
 *        x y z for all outer nodes
 *    \endcode
 *
 *
 *    The file fileName.ng
 *    - contains boundary information following the letter B
 *    - lists inner nodes with their coordinates xyz following the letter I
 *    - gives element information following the letter E, where
 *
 *    - for simplices, there are 4 nodes.
 *    - for pyramids, there are 5 nodes.
 *    - for prisms, there are 6 nodes.
 *    - for cubes, there are 8 nodes.
 *
 *    following the letters F the nodes which are part of boundary surfaces are listed.
 *
 *    This reader only supports three-dimensional grids.
 *
 *    Currently no boundary element data is passed to \a grid.
 */
template <class GridType>
class LgmReader
{

public:

    /** \brief Read grid from a Lgm file
     *    \param fileName The base file name of the lgm files
     *    \param verbose Tlag to set whether information should be printed
     */
    static GridType* read(const std::string& fileName, bool verbose = true)
    {
        // extract the grid dimension
        const int dim = GridType::dimension;


        // currently only dim = 3 is implemented
        if (dim != 3)
            DUNE_THROW(Dune::NotImplemented,
                       "Reading lgm format is not implemented for dimension " << dim);

        // set up the grid factory
        GridFactory<GridType> factory;

        // outer vertices are placed in the .lgm file after #Point-Info, inner vertices in the .ng file following the initial I

        // set the name of the outer vertex file
        std::string vertexFileName = fileName + ".lgm";

        // set the vertex input stream
        std::ifstream vertexFile(vertexFileName.c_str());
        if (!vertexFile)
            DUNE_THROW(Dune::IOError, "Could not open " << vertexFileName);
        else if (verbose)
            std::cout << " lgm file  was opened successfully" << std::endl;


        // read the outer vertices

        std::string jump;
        int numberOfVertices = 0;
        string inputfile = vertexFileName;
        Dune::FieldVector<double,dim> position;
        int k=0;
        int anz = 0;
        std::vector<Dune:: FieldVector <double, dim> > coordinates(0);


        // jump everything before "#Point-Info"
        while (vertexFile >> jump)
        {
            if (jump == "#Point-Info")
                break;
        }

        // read the node data and insert nodes
        while (vertexFile >> position[k])
        {
            numberOfVertices++;
            for ( k = 1; k < dim ; k++)
                vertexFile >> position[k];

            factory.insertVertex(position);
            //            if (verbose)
            //                        std::cout << " insert vertex " << position << std::endl;

            vertexFile >> jump;
            k=0;
        }


        if (verbose)
            std::cout << numberOfVertices << " all vertices in lgm read." << std::endl;



        // continue with the element file:


        // set the name of the element file
        std::string elementFileName = fileName + ".ng";

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
        bool firstElement = true;
        bool followingE = true;


        if (verbose)
            std::cout << " ng file was opened successfully" << std::endl;

        int maxNumberOfVertices = 8;
        std::vector<unsigned int> vertices(maxNumberOfVertices+2);
        string vert;
        int subdomain;


        while (elementFile >> jump)
        {
            // do nothing with the information in the first lines

            //read inner vertices from the element file .ng
            if (jump == "I")
            {
                numberOfVertices++;

                // read the node data
                for ( k = 0; k < dim ; k++)
                    elementFile >> position[k];

                factory.insertVertex(position);
                //                   if (verbose)
                //                       std::cout << " insert vertex I " << position << std::endl;

                vertexFile >> jump;  // get rid of the semicolon

            }

            // first found E

            if (jump == "E")
            {
                followingE = true;
                if (verbose)
                    std::cout << numberOfVertices << " vertices read in total." << std::endl;
                break;
            }
        }

        //read element information

        while (elementFile >> jump)

        {
            if (firstElement == true) // the first time jump contains the subdomain, not an E
            {
                k=0;
                numberOfElements++;
                followingE = true;
                firstElement = false;
                elementFile >> jump;  // read first node of first element
            }

            if (jump == "F")
                followingE = false;

            else if (jump == "E")  // enters here after every E except after the first one
            {
                // insert vertices

                // simplex

                if (anz == 4)
                {
                    numberOfSimplices++;
                    std::vector<unsigned int> simplexVertices(4);

                    for (int k = 0; k < 4; k++)
                        simplexVertices[k] = vertices[k];

                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim), simplexVertices);
                }


                //  pyramid

                else if (anz == 5)
                {
                    numberOfPyramids++;
                    std::vector<unsigned int> pyramidVertices(5);

                    for (int k = 0; k < 5; k++)
                        pyramidVertices[k] = vertices[k];

                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::pyramid,dim), pyramidVertices);
                }


                // prism

                else if(anz ==6)
                {
                    numberOfPrisms++;
                    std::vector<unsigned int> prismVertices(6);

                    for (int k = 0; k < 6; k++)
                        prismVertices[k] = vertices[k];

                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism,dim), prismVertices);

                }

                else if (anz==7)
                    DUNE_THROW(Dune::IOError, "Not an Element ");


                // cube

                else if (anz==8)
                {
                    numberOfCubes++;
                    std::vector<unsigned int> cubeVertices(8);

                    for (int k = 0; k < 8; k++)
                        cubeVertices[k] = vertices[k];

                    std::swap(cubeVertices[2], cubeVertices[3]);
                    std::swap(cubeVertices[6], cubeVertices[7]);

                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim), cubeVertices);
                }

                else
                    DUNE_THROW(Dune::IOError, "Not an Element ");

                // read next element

                elementFile >> subdomain;    // the first number following an "E" tells the subdomain
                anz=0;
                numberOfElements++;
                followingE = true;
                firstElement = false;

            }

            else if (followingE == true) // read the number of a node
            {
                vertices[anz] = atoi(jump.c_str());
                anz = anz+1;
            }

        }

        //insert vertices of the last element


        // simplex

        if (anz == 4)
        {
            numberOfSimplices++;
            std::vector<unsigned int> simplexVertices(4);
            for (int k = 0; k < 4; k++)
                simplexVertices[k] = vertices[k];

            factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim), simplexVertices);
        }


        //  pyramid

        else if (anz == 5)
        {
            numberOfPyramids++;
            std::vector<unsigned int> pyramidVertices(5);
            for (int k = 0; k < 5; k++)
                pyramidVertices[k] = vertices[k];

            factory.insertElement(Dune::GeometryType(Dune::GeometryType::pyramid,dim), pyramidVertices);
        }


        // prism

        else if(anz ==6)
        {
            numberOfPrisms++;
            std::vector<unsigned int> prismVertices(6);
            for (int k = 0; k < 6; k++)
                prismVertices[k] = vertices[k];

            factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism,dim), prismVertices);
        }

        else if (anz==7)
            DUNE_THROW(Dune::IOError, "Not an Element ");


        // cube

        else if (anz==8)
        {
            numberOfCubes++;
            std::vector<unsigned int> cubeVertices(8);

            for (int k = 0; k < 8; k++)
                cubeVertices[k] = vertices[k];

            std::swap(cubeVertices[2], cubeVertices[3]);
            std::swap(cubeVertices[6], cubeVertices[7]);

            factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim), cubeVertices);
        }

        else
            DUNE_THROW(Dune::IOError, "Not an Element ");

        // final output

        if (verbose)
            std::cout << numberOfElements << " elements read: "
                      << numberOfSimplices << " simplices, " << numberOfPyramids << " pyramids, "
                      << numberOfPrisms << " prisms, " << numberOfCubes << " cubes." << std::endl;



        // finish the construction of the grid object
        if (verbose)
            std::cout << "Starting createGrid() ... " << std::flush;

        return factory.createGrid();

    }

};

}


#endif
