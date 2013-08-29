#ifndef DUMUX_INTERFACEMESHCREATOR_HH
#define DUMUX_INTERFACEMESHCREATOR_HH

#include<dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>

namespace Dumux
{

template <class Grid>
class InterfaceMeshCreator
{
public:
    enum {dim = Grid::dimension};
    typedef typename Grid::ctype Scalar;
    typedef typename Grid::LeafGridView GridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;

    static Grid* create(const std::string& dgfName, const Dune::FieldVector<int, dim>& nElements,
                        const Scalar interfaceY, const Scalar gradingFactor, const bool refineTop = false)
    {
	typedef Dune::SGrid<dim,dim> HelperGrid;
        Dune::GridPtr<HelperGrid> helperGridPtr(dgfName);
        HelperGrid& helperGrid = *helperGridPtr;
        typedef typename HelperGrid::LeafGridView HelperGridView;
        typedef typename HelperGridView::template Codim<0>::Iterator HelperElementIterator;
        typedef typename HelperGridView::Traits::template Codim<0>::Entity HelperElement;
        typedef typename HelperElement::Geometry HelperGeometry;

        HelperElementIterator helperElementIterator = helperGrid.leafView().template begin<0>();
        const HelperElement& helperElement = *helperElementIterator;
        const HelperGeometry& helperGeometry = helperElement.geometry();

        Dune::FieldVector<Scalar,dim> refinePoint(0);
	refinePoint[dim-1] = interfaceY;
        Dune::FieldVector<Scalar,dim> gradingFactors(1);
        gradingFactors[dim-1] = gradingFactor;

        Dune::FieldVector<Scalar,dim> refinePointLocal(helperGeometry.local(refinePoint));
std::cout << "rglobal = " << refinePoint << ", rlocal = " << refinePointLocal << std::endl;
        Dune::GridFactory<Grid> factory;

        int nX = nElements[0];
        int nY = nElements[1];

        std::vector<std::vector<Scalar> > localPositions(dim);
        for (int comp = 0; comp < dim; comp++)
        {
            Scalar lengthLeft = refinePointLocal[comp];
            Scalar lengthRight = 1.0 - lengthLeft;

            int nLeft, nRight;
            Scalar hLeft, hRight;
            if (lengthLeft < 1e-10)
            {
                nLeft = 0;
                nRight = nElements[comp];

                if (gradingFactors[comp] > 1.0)
                    hLeft = hRight = (1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nRight));
                else
                    hLeft = hRight = 1.0/nElements[comp];
            }
            else if (lengthLeft > 1.0 - 1e-10)
            {
                nLeft = nElements[comp];
                nRight = 0;

                if (gradingFactors[comp] > 1.0)
                    hLeft = hRight = (1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nLeft));
                else
                    hLeft = hRight = 1.0/nElements[comp];
            }
            else if (comp == dim - 1 && refineTop)
            {
                lengthLeft = refinePointLocal[comp];
                lengthRight = (1 - refinePointLocal[comp])/2;

                nLeft = nRight = nElements[comp]/3;

                if (nElements[comp]%3 == 1)
                    nLeft += 1;
                else if (nElements[comp]%3 == 2)
                    nRight += 1;

                hLeft = lengthLeft*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nLeft));
                hRight = lengthRight*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nRight));
            }
            else if (lengthLeft > 0.5)
            {
                Scalar nLeftDouble = std::ceil(-log((1.0 + sqrt(1.0 + 4.0*pow(gradingFactors[comp], nElements[comp])*lengthRight/lengthLeft))
                            /(2.0*pow(gradingFactors[comp], nElements[comp])))/log(gradingFactors[comp]));
                nLeft = std::min((int)std::ceil(nLeftDouble), nElements[comp]);

                nRight = nElements[comp] - nLeft;

                if (gradingFactors[comp] > 1.0)
                {
                    hLeft = lengthLeft*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nLeft));
                    hRight = lengthRight*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nRight));
                }
                else
                    hLeft = hRight = 1.0/nElements[comp];
            }
            else
            {
                Scalar nRightDouble = -log((1.0 + sqrt(1.0 + 4.0*pow(gradingFactors[comp], nElements[comp])*lengthLeft/lengthRight))
                            /(2.0*pow(gradingFactors[comp], nElements[comp])))/log(gradingFactors[comp]);
                nRight = std::min((int)std::ceil(nRightDouble), nElements[comp]);

                nLeft = nElements[comp] - nRight;

                if (gradingFactors[comp] > 1.0)
                {
                    hLeft = lengthLeft*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nLeft));
                    hRight = lengthRight*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nRight));
                }
                else
                    hLeft = hRight = 1.0/nElements[comp];
            }
std::cout << "lengthLeft = " << lengthLeft << ", lengthRight = " << lengthRight << ", hLeft = " << hLeft << ", hRight = " << hRight << ", nLeft = " << nLeft << ", nRight = " << nRight << std::endl;

            int nVertices = nElements[comp] + 1;
            localPositions[comp].resize(nVertices);

            localPositions[comp][0] = 0.0;
            for (int i = 0; i < nLeft; i++)
            {
                Scalar hI = hLeft*pow(gradingFactors[comp], nLeft-1-i);
                localPositions[comp][i+1] = localPositions[comp][i] + hI;
            }

            for (int i = 0; i < nRight; i++)
            {
                Scalar hI = hRight*pow(gradingFactors[comp], i);
                localPositions[comp][nLeft+i+1] = localPositions[comp][nLeft+i] + hI;
            }

            if (comp == dim - 1 && refineTop)
                for (int i = 0; i < nRight; i++)
                {
                    Scalar hI = hRight*pow(gradingFactors[comp], nRight-1-i);
                    localPositions[comp][nLeft+nRight+i+1] = localPositions[comp][nLeft+nRight+i] + hI;
                }

            if (localPositions[comp][nVertices-1] != 1.0)
            {
                for (int i = 0; i < nVertices; i++)
                    localPositions[comp][i] /= localPositions[comp][nVertices-1];
            }
        }

        Dune::FieldVector<Scalar,dim> local;
        for (int j = 0; j < nY + 1; j++)
            {
                local[1] = localPositions[1][j];
                for (int i = 0; i < nX + 1; i++)
                {
                    local[0] = localPositions[0][i];

                    Dune::FieldVector<Scalar,dim> position(helperGeometry.global(local));
                    factory.insertVertex(position);

                }
            }

        for (int j = 0; j < nY; j++)
            {
                for (int i = 0; i < nX; i++)
                {
                    std::vector<unsigned int> vertices(4);

                    vertices[0] = j*(nX+1) + i;
                    vertices[1] = j*(nX+1) + i+1;
                    vertices[2] = (j+1)*(nX+1) + i;
                    vertices[3] = (j+1)*(nX+1) + i+1;
                    
                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim), vertices);

                }
            }

        return factory.createGrid();
    }
};
}


#endif
