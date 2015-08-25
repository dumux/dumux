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
 * \brief A grid creator for the coupled problems, with a refined interface
 */

#ifndef DUMUX_INTERFACEGRIDCREATOR_HH
#define DUMUX_INTERFACEGRIDCREATOR_HH

#include <dune/common/deprecated.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/common/basicproperties.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
}

/*!
 * \brief A grid creator for the coupled problems, with a refined interface
 *
 * A grid creator, which can refine the grid towards the
 * coupling interface and the top of the domain.
 */
template<class TypeTag>
class InterfaceGridCreator
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef std::shared_ptr<Grid> GridPointer;
    enum {dim = Grid::dimension};

    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        if (dim != 2)
        {
            DUNE_THROW(Dune::NotImplemented, "The InterfaceGridCreator is not implemented for 1D and 3D.");
        }

        Dune::array<unsigned int, dim> numCellsDummy = {{1,1}};
        Dune::array<unsigned int, dim> numCells;
        Dune::FieldVector<Scalar, dim> lowerLeft;
        Dune::FieldVector<Scalar, dim> upperRight;
        Dune::FieldVector<Scalar, dim> refinePoint(0.0);
        Dune::FieldVector<Scalar, dim> gradingFactor(1.0);
        Dune::FieldVector<bool, dim> refineTop(false);

        // x-direction
        numCells[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, unsigned int, Grid, NumberOfCellsX);
        lowerLeft[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, LowerLeftX);
        upperRight[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightX);
        try { refinePoint[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosX); }
        catch (Dumux::ParameterException &e) { }
        try { gradingFactor[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, GradingFactorX); }
        catch (Dumux::ParameterException &e) { }
        try { refineTop[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Grid, RefineTopX); }
        catch (Dumux::ParameterException &e) { }
        // y-direction
        numCells[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, unsigned int, Grid, NumberOfCellsY);
        lowerLeft[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, LowerLeftY);
        upperRight[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightY);
        try { refinePoint[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY); }
        catch (Dumux::ParameterException &e) { }
        try { gradingFactor[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, GradingFactorY); }
        catch (Dumux::ParameterException &e) { }
        try { refineTop[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Grid, RefineTopY); }
        catch (Dumux::ParameterException &e) { }

        typedef Dune::YaspGrid<dim> HelperGrid;
        std::shared_ptr<HelperGrid> helperGrid = std::shared_ptr<HelperGrid> (
            Dune::StructuredGridFactory<HelperGrid>::createCubeGrid(lowerLeft, upperRight, numCellsDummy));

        typedef typename HelperGrid::LeafGridView HelperGridView;
        typedef typename HelperGridView::template Codim<0>::Iterator HelperElementIterator;
        typedef typename HelperGridView::Traits::template Codim<0>::Entity HelperElement;
        typedef typename HelperElement::Geometry HelperGeometry;

        HelperElementIterator helperElementIterator = helperGrid->leafGridView().template begin<0>();
        const HelperElement& helperElement = *helperElementIterator;
        const HelperGeometry& helperGeometry = helperElement.geometry();

        Dune::FieldVector<Scalar,dim> refinePointLocal(helperGeometry.local(refinePoint));
        std::cout << "refinePointGlobal = " << refinePoint
                  << ", refinePointLocal = " << refinePointLocal
                  << std::endl;
        Dune::GridFactory<Grid> factory;

        int nX = numCells[0];
        int nY = numCells[1];

        std::vector<std::vector<Scalar> > localPositions(dim);
        for (int dimIdx = 0; dimIdx < dim; dimIdx++)
        {
            Scalar lengthLeft = refinePointLocal[dimIdx];
            Scalar lengthRight = 1.0 - lengthLeft;

            int nLeft, nRight;
            Scalar hLeft, hRight;
            if (lengthLeft < 1e-10)
            {
                nLeft = 0;
                nRight = numCells[dimIdx];

                if (gradingFactor[dimIdx] > 1.0)
                    hLeft = hRight = (1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nRight));
                else
                    hLeft = hRight = 1.0/numCells[dimIdx];
            }
            else if (lengthLeft > 1.0 - 1e-10)
            {
                nLeft = numCells[dimIdx];
                nRight = 0;

                if (gradingFactor[dimIdx] > 1.0)
                    hLeft = hRight = (1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nLeft));
                else
                    hLeft = hRight = 1.0/numCells[dimIdx];
            }
            else if (refineTop[dimIdx])
            {
                lengthLeft = refinePointLocal[dimIdx];
                lengthRight = (1 - refinePointLocal[dimIdx])/2;

                nLeft = nRight = numCells[dimIdx]/3;

                if (numCells[dimIdx]%3 == 1)
                    nLeft += 1;
                else if (numCells[dimIdx]%3 == 2)
                    nRight += 1;

                hLeft = lengthLeft*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nLeft));
                hRight = lengthRight*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nRight));
            }
            else if (lengthLeft > 0.5)
            {
                Scalar nLeftDouble = std::ceil(-log((1.0 + sqrt(1.0 + 4.0 * pow(gradingFactor[dimIdx], numCells[dimIdx])
                                                                      * lengthRight/lengthLeft))
                            /(2.0*pow(gradingFactor[dimIdx], numCells[dimIdx])))/log(gradingFactor[dimIdx]));
                nLeft = std::min((unsigned int)std::ceil(nLeftDouble), numCells[dimIdx]);

                nRight = numCells[dimIdx] - nLeft;

                if (gradingFactor[dimIdx] > 1.0)
                {
                    hLeft = lengthLeft*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nLeft));
                    hRight = lengthRight*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nRight));
                }
                else
                    hLeft = hRight = 1.0/numCells[dimIdx];
            }
            else
            {
                Scalar nRightDouble = -log((1.0 + sqrt(1.0 + 4.0 * pow(gradingFactor[dimIdx], numCells[dimIdx])
                                                             * lengthLeft/lengthRight))
                            /(2.0*pow(gradingFactor[dimIdx], numCells[dimIdx])))/log(gradingFactor[dimIdx]);
                nRight = std::min((unsigned int)std::ceil(nRightDouble), numCells[dimIdx]);

                nLeft = numCells[dimIdx] - nRight;

                if (gradingFactor[dimIdx] > 1.0)
                {
                    hLeft = lengthLeft*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nLeft));
                    hRight = lengthRight*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nRight));
                }
                else
                    hLeft = hRight = 1.0/numCells[dimIdx];
            }
            std::cout << "lengthLeft = " << lengthLeft
                      << ", lengthRight = " << lengthRight
                      << ", hLeft = " << hLeft
                      << ", hRight = " << hRight
                      << ", nLeft = " << nLeft
                      << ", nRight = " << nRight
                      << std::endl;

            int numVertices = numCells[dimIdx] + 1;
            localPositions[dimIdx].resize(numVertices);

            localPositions[dimIdx][0] = 0.0;
            for (int i = 0; i < nLeft; i++)
            {
                Scalar hI = hLeft*pow(gradingFactor[dimIdx], nLeft-1-i);
                localPositions[dimIdx][i+1] = localPositions[dimIdx][i] + hI;
            }

            for (int i = 0; i < nRight; i++)
            {
                Scalar hI = hRight*pow(gradingFactor[dimIdx], i);
                localPositions[dimIdx][nLeft+i+1] = localPositions[dimIdx][nLeft+i] + hI;
            }

            if (refineTop[dimIdx])
                for (int i = 0; i < nRight; i++)
                {
                    Scalar hI = hRight*pow(gradingFactor[dimIdx], nRight-1-i);
                    localPositions[dimIdx][nLeft+nRight+i+1] = localPositions[dimIdx][nLeft+nRight+i] + hI;
                }

            if (localPositions[dimIdx][numVertices-1] != 1.0)
            {
                for (int i = 0; i < numVertices; i++)
                    localPositions[dimIdx][i] /= localPositions[dimIdx][numVertices-1];
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

        gridPtr().reset(factory.createGrid());
    }

    DUNE_DEPRECATED_MSG("create() is deprecated use makeGrid() instead")
    static Grid* create(const std::string& dgfName, const Dune::FieldVector<int, dim>& numCells,
                        const Scalar interfaceY, const Scalar grading, const bool refineTop = false)
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
        Dune::FieldVector<Scalar,dim> gradingFactor(1);
        gradingFactor[dim-1] = grading;

        Dune::FieldVector<Scalar,dim> refinePointLocal(helperGeometry.local(refinePoint));
        std::cout << "rglobal = " << refinePoint
                  << ", rlocal = " << refinePointLocal
                  << std::endl;
        Dune::GridFactory<Grid> factory;

        int nX = numCells[0];
        int nY = numCells[1];

        std::vector<std::vector<Scalar> > localPositions(dim);
        for (int dimIdx = 0; dimIdx < dim; dimIdx++)
        {
            Scalar lengthLeft = refinePointLocal[dimIdx];
            Scalar lengthRight = 1.0 - lengthLeft;

            int nLeft, nRight;
            Scalar hLeft, hRight;
            if (lengthLeft < 1e-10)
            {
                nLeft = 0;
                nRight = numCells[dimIdx];

                if (gradingFactor[dimIdx] > 1.0)
                    hLeft = hRight = (1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nRight));
                else
                    hLeft = hRight = 1.0/numCells[dimIdx];
            }
            else if (lengthLeft > 1.0 - 1e-10)
            {
                nLeft = numCells[dimIdx];
                nRight = 0;

                if (gradingFactor[dimIdx] > 1.0)
                    hLeft = hRight = (1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nLeft));
                else
                    hLeft = hRight = 1.0/numCells[dimIdx];
            }
            else if (dimIdx == dim - 1 && refineTop)
            {
                lengthLeft = refinePointLocal[dimIdx];
                lengthRight = (1 - refinePointLocal[dimIdx])/2;

                nLeft = nRight = numCells[dimIdx]/3;

                if (numCells[dimIdx]%3 == 1)
                    nLeft += 1;
                else if (numCells[dimIdx]%3 == 2)
                    nRight += 1;

                hLeft = lengthLeft*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nLeft));
                hRight = lengthRight*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nRight));
            }
            else if (lengthLeft > 0.5)
            {
                Scalar nLeftDouble = std::ceil(-log((1.0 + sqrt(1.0 + 4.0 * pow(gradingFactor[dimIdx], numCells[dimIdx])
                                                                      * lengthRight/lengthLeft))
                            /(2.0*pow(gradingFactor[dimIdx], numCells[dimIdx])))/log(gradingFactor[dimIdx]));
                nLeft = std::min((int)std::ceil(nLeftDouble), numCells[dimIdx]);

                nRight = numCells[dimIdx] - nLeft;

                if (gradingFactor[dimIdx] > 1.0)
                {
                    hLeft = lengthLeft*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nLeft));
                    hRight = lengthRight*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nRight));
                }
                else
                    hLeft = hRight = 1.0/numCells[dimIdx];
            }
            else
            {
                Scalar nRightDouble = -log((1.0 + sqrt(1.0 + 4.0 * pow(gradingFactor[dimIdx], numCells[dimIdx])
                                                             * lengthLeft/lengthRight))
                            /(2.0*pow(gradingFactor[dimIdx], numCells[dimIdx])))/log(gradingFactor[dimIdx]);
                nRight = std::min((int)std::ceil(nRightDouble), numCells[dimIdx]);

                nLeft = numCells[dimIdx] - nRight;

                if (gradingFactor[dimIdx] > 1.0)
                {
                    hLeft = lengthLeft*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nLeft));
                    hRight = lengthRight*(1.0 - gradingFactor[dimIdx])/(1.0 - pow(gradingFactor[dimIdx], nRight));
                }
                else
                    hLeft = hRight = 1.0/numCells[dimIdx];
            }
            std::cout << "lengthLeft = " << lengthLeft
                      << ", lengthRight = " << lengthRight
                      << ", hLeft = " << hLeft
                      << ", hRight = " << hRight
                      << ", nLeft = " << nLeft
                      << ", nRight = " << nRight
                      << std::endl;

            int numVertices = numCells[dimIdx] + 1;
            localPositions[dimIdx].resize(numVertices);

            localPositions[dimIdx][0] = 0.0;
            for (int i = 0; i < nLeft; i++)
            {
                Scalar hI = hLeft*pow(gradingFactor[dimIdx], nLeft-1-i);
                localPositions[dimIdx][i+1] = localPositions[dimIdx][i] + hI;
            }

            for (int i = 0; i < nRight; i++)
            {
                Scalar hI = hRight*pow(gradingFactor[dimIdx], i);
                localPositions[dimIdx][nLeft+i+1] = localPositions[dimIdx][nLeft+i] + hI;
            }

            if (dimIdx == dim - 1 && refineTop)
                for (int i = 0; i < nRight; i++)
                {
                    Scalar hI = hRight*pow(gradingFactor[dimIdx], nRight-1-i);
                    localPositions[dimIdx][nLeft+nRight+i+1] = localPositions[dimIdx][nLeft+nRight+i] + hI;
                }

            if (localPositions[dimIdx][numVertices-1] != 1.0)
            {
                for (int i = 0; i < numVertices; i++)
                    localPositions[dimIdx][i] /= localPositions[dimIdx][numVertices-1];
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

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *gridPtr();
    }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    static void loadBalance()
    {
        gridPtr()->loadBalance();
    }

    /*!
     * \brief Returns a reference to the shared pointer to the grid.
     */
    static GridPointer &gridPtr()
    {
        static GridPointer interfaceGridCreator;
        return interfaceGridCreator;
    }
};
}


#endif
