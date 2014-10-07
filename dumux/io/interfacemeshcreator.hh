// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A mesh creator for the coupled problems
 */

#ifndef DUMUX_INTERFACEMESHCREATOR_HH
#define DUMUX_INTERFACEMESHCREATOR_HH

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

namespace Dumux
{
/*!
 * \brief A mesh creator for the coupled problems
 *
 * A mesh creator, that creates a grid, which can be refined towards the
 * coupling interface and the top of the domain.
 */

template <class Grid>
class InterfaceMeshCreator
{
public:
    enum {dim = Grid::dimension};
    typedef typename Grid::ctype Scalar;
    typedef typename Grid::LeafGridView GridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;

    static Grid* create(const std::string& dgfName, const Dune::FieldVector<int, dim>& numElements,
                        const Scalar interfaceY, const Scalar gradingFactor, const bool refineTop = false)
    {
        typedef Dune::YaspGrid<dim> HelperGrid;
        Dune::GridPtr<HelperGrid> helperGridPtr(dgfName);
        HelperGrid& helperGrid = *helperGridPtr;
        typedef typename HelperGrid::LeafGridView HelperGridView;
        typedef typename HelperGridView::template Codim<0>::Iterator HelperElementIterator;
        typedef typename HelperGridView::Traits::template Codim<0>::Entity HelperElement;
        typedef typename HelperElement::Geometry HelperGeometry;

        HelperElementIterator helperElementIterator = helperGrid.leafGridView().template begin<0>();
        const HelperElement& helperElement = *helperElementIterator;
        const HelperGeometry& helperGeometry = helperElement.geometry();

        Dune::FieldVector<Scalar,dim> refinePoint(0);
	refinePoint[dim-1] = interfaceY;
        Dune::FieldVector<Scalar,dim> gradingFactors(1);
        gradingFactors[dim-1] = gradingFactor;

        Dune::FieldVector<Scalar,dim> refinePointLocal(helperGeometry.local(refinePoint));
std::cout << "rglobal = " << refinePoint << ", rlocal = " << refinePointLocal << std::endl;
        Dune::GridFactory<Grid> factory;

        int nX = numElements[0];
        int nY = numElements[1];

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
                nRight = numElements[comp];

                if (gradingFactors[comp] > 1.0)
                    hLeft = hRight = (1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nRight));
                else
                    hLeft = hRight = 1.0/numElements[comp];
            }
            else if (lengthLeft > 1.0 - 1e-10)
            {
                nLeft = numElements[comp];
                nRight = 0;

                if (gradingFactors[comp] > 1.0)
                    hLeft = hRight = (1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nLeft));
                else
                    hLeft = hRight = 1.0/numElements[comp];
            }
            else if (comp == dim - 1 && refineTop)
            {
                lengthLeft = refinePointLocal[comp];
                lengthRight = (1 - refinePointLocal[comp])/2;

                nLeft = nRight = numElements[comp]/3;

                if (numElements[comp]%3 == 1)
                    nLeft += 1;
                else if (numElements[comp]%3 == 2)
                    nRight += 1;

                hLeft = lengthLeft*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nLeft));
                hRight = lengthRight*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nRight));
            }
            else if (lengthLeft > 0.5)
            {
                Scalar nLeftDouble = std::ceil(-log((1.0 + sqrt(1.0 + 4.0 * pow(gradingFactors[comp], numElements[comp])
                                                                      * lengthRight/lengthLeft))
                            /(2.0*pow(gradingFactors[comp], numElements[comp])))/log(gradingFactors[comp]));
                nLeft = std::min((int)std::ceil(nLeftDouble), numElements[comp]);

                nRight = numElements[comp] - nLeft;

                if (gradingFactors[comp] > 1.0)
                {
                    hLeft = lengthLeft*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nLeft));
                    hRight = lengthRight*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nRight));
                }
                else
                    hLeft = hRight = 1.0/numElements[comp];
            }
            else
            {
                Scalar nRightDouble = -log((1.0 + sqrt(1.0 + 4.0 * pow(gradingFactors[comp], numElements[comp])
                                                             * lengthLeft/lengthRight))
                            /(2.0*pow(gradingFactors[comp], numElements[comp])))/log(gradingFactors[comp]);
                nRight = std::min((int)std::ceil(nRightDouble), numElements[comp]);

                nLeft = numElements[comp] - nRight;

                if (gradingFactors[comp] > 1.0)
                {
                    hLeft = lengthLeft*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nLeft));
                    hRight = lengthRight*(1.0 - gradingFactors[comp])/(1.0 - pow(gradingFactors[comp], nRight));
                }
                else
                    hLeft = hRight = 1.0/numElements[comp];
            }
std::cout << "lengthLeft = " << lengthLeft << ", lengthRight = " << lengthRight << ", hLeft = " << hLeft <<
             ", hRight = " << hRight << ", nLeft = " << nLeft << ", nRight = " << nRight << std::endl;

            int numVertices = numElements[comp] + 1;
            localPositions[comp].resize(numVertices);

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

            if (localPositions[comp][numVertices-1] != 1.0)
            {
                for (int i = 0; i < numVertices; i++)
                    localPositions[comp][i] /= localPositions[comp][numVertices-1];
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
