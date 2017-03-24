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
 * \brief Provides a grid creator for a piece of cake grid
 */

#ifndef DUMUX_CAKE_GRID_CREATOR_HH
#define DUMUX_CAKE_GRID_CREATOR_HH

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <dune/common/dynvector.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dumux/common/basicproperties.hh>
#include <dumux/common/propertysystem.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
}

/*!
 * \brief Provides a grid creator with a method for creating creating vectors
 *        with polar Coordinates and one for creating a cartesian grid from
 *        these polar coordinates.
 */
template <class TypeTag>
class CakeGridCreator
{
    using ThisType = CakeGridCreator<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridFactory = Dune::GridFactory<Grid>;
    using GridPointer = std::shared_ptr<Grid>;

    enum { dim = Grid::dimension,
           dimWorld = Grid::dimensionworld
         };

    using ScalarVector = std::vector<Scalar>;
    using CellArray = std::array<unsigned int, dim>;
    using GlobalPosition = Dune::FieldVector<typename Grid::ctype, dimWorld>;

public:
    /*!
     * \brief Make the grid.
     */
    static void makeGrid()
    {
        static_assert(dim == 2 || dim == 3, "The CakeGridCreator is only implemented for 2D and 3D.");

        bool verbose = false;
        try { verbose = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, bool, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), Verbosity);}
        catch (Dumux::ParameterException &e) {  }

        std::array<ScalarVector, dim> polarCoordinates;
        // Indices specifing in which direction the piece of cake is oriented
        Dune::FieldVector<int, dim> indices(-1);

        createVectors(polarCoordinates, indices, verbose);
        gridPtr() = ThisType::createCakeGrid(polarCoordinates, indices, verbose);
    }

    /*!
    * \brief Create vectors containing polar coordinates of all points.
    *
    * All keys are expected to be in group GridParameterGroup.
    * The following keys are recognized:
    * - Radial : min/max value for radial coordinate
    * - Angular : min/max value for angular coordinate
    * - Axial : min/max value for axial coordinate
    *   Adding 0, 1 (or 2 in 3D) specifies in which direction (x, y and z, respectively)
    *   the radial, angular and axial direction are oriented
    * - Cells : number of cells array for x-coordinate (Again, an added 0, 1 or 3 specifies x, y or z
    * - Grading : grading factor array for x-coordinate (Same here)
    * - Verbosity : whether the grid construction should output to standard out
    *
    * The grading factor \f$ g \f$ specifies the ratio between the next and the current cell size:
    * \f$ g = \frac{h_{i+1}}{h_i} \f$.
    * Negative grading factors are converted to
    * \f$ g = -\frac{1}{g_\textrm{negative}} \f$
    * to avoid issues with imprecise fraction numbers.
    */
    static void createVectors(std::array<ScalarVector, dim> &polarCoordinates,
                              Dune::FieldVector<int, dim> &indices,
                              bool verbose = false)
    {
        // Only construction from the input file is possible
        // Look for the necessary keys to construct from the input file
        try
        {
            // The positions
            std::array<ScalarVector, dim> positions;

            for (int i = 0; i < dim; ++i)
            {
                const auto paramNameRadial = GET_PROP_VALUE(TypeTag, GridParameterGroup)
                                             + ".Radial" + std::to_string(i);
                try
                {
                    positions[i] = GET_RUNTIME_PARAM_CSTRING(TypeTag, ScalarVector, paramNameRadial.c_str());
                    indices[0] = i; // Index specifing radial direction
                }
                catch (Dumux::ParameterException &e) { }

                const auto paramNameAngular = GET_PROP_VALUE(TypeTag, GridParameterGroup)
                                              + ".Angular" + std::to_string(i);
                try
                {
                    positions[i] = GET_RUNTIME_PARAM_CSTRING(TypeTag, ScalarVector, paramNameAngular.c_str());
                    indices[1] = i; // Index specifing angular direction
                }
                catch (Dumux::ParameterException &e) { }

                if (dim == 3)
                {
                    const auto paramNameAxial = GET_PROP_VALUE(TypeTag, GridParameterGroup)
                                                + ".Axial" + std::to_string(i);
                    try
                    {
                        positions[i] = GET_RUNTIME_PARAM_CSTRING(TypeTag, ScalarVector, paramNameAxial.c_str());
                        indices[2] = i; // Index specifing axial direction
                    }
                    catch (Dumux::ParameterException &e) { }
                }
            }

            if (dim == 3)
            {
                if ( !( ((3 > indices[0]) && (indices[0] >= 0)
                    &&   (3 > indices[1]) && (indices[1] >= 0)
                    &&   (3 > indices[2]) && (indices[2] >= 0))
                    && ((indices[2] != indices[1]) && (indices[2] != indices[0])
                    && (indices[1] != indices[0]))
                    && (indices[2] + indices[1] + indices[0] == 3) ) )
                {
                    if (verbose)
                        std::cout << " indices[0] " << indices[0]
                                  << " indices[1] " << indices[1]
                                  << " indices[2] " << indices[2] << std::endl;
                    DUNE_THROW(Dune::RangeError, "Please specify Positions Angular and Radial correctly and unambiguously!" << std::endl);
                }
            }
            else
            {
                if ( !( ((2 > indices[0]) && (indices[0] >= 0)
                    &&   (2 > indices[1]) && (indices[1] >= 0))
                    &&   (indices[1] != indices[0])
                    &&   (indices[1] + indices[0] == 1) ) )
                {
                    if (verbose)
                        std::cout << " indices[0] " << indices[0]
                                  << " indices[1] " << indices[1] << std::endl;
                    DUNE_THROW(Dune::RangeError, "Please specify Positions Angular and Radial correctly and unambiguously!" << std::endl);
                }
            }

            // the number of cells (has a default)
            std::array<std::vector<int>, dim> cells;
            for (int i = 0; i < dim; ++i)
            {
                const auto paramName = GET_PROP_VALUE(TypeTag, GridParameterGroup)
                                       + ".Cells" + std::to_string(i);
                try { cells[i] = GET_RUNTIME_PARAM_CSTRING(TypeTag, std::vector<int>, paramName.c_str()); }
                catch (Dumux::ParameterException &e) { cells[i].resize(positions[i].size()-1, 1.0); }
            }

            // grading factor (has a default)
            std::array<ScalarVector, dim> grading;
            for (int i = 0; i < dim; ++i)
            {
                const auto paramName = GET_PROP_VALUE(TypeTag, GridParameterGroup)
                                       + ".Grading" + std::to_string(i);
                try { grading[i] = GET_RUNTIME_PARAM_CSTRING(TypeTag, ScalarVector, paramName.c_str()); }
                catch (Dumux::ParameterException &e) { grading[i].resize(positions[i].size()-1, 1.0); }
            }

            // make the grid
            // sanity check of the input parameters
            for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                // check that cells and position have same size
                if (cells[dimIdx].size() + 1 != positions[dimIdx].size())
                    DUNE_THROW(Dune::RangeError, "Make sure to specify correct \"Cells\" and \"Positions\" arrays");

                // check that grading and position have same size
                else if (grading[dimIdx].size() + 1 != positions[dimIdx].size())
                    DUNE_THROW(Dune::RangeError, "Make sure to specify correct \"Grading\" and \"Positions\" arrays");

                // make sure we have monotonically increasing positions
                for (unsigned int posIdx = 0; posIdx < positions[dimIdx].size(); ++posIdx)
                    if (!std::is_sorted(positions[dimIdx].begin(), positions[dimIdx].end()))
                        DUNE_THROW(Dune::GridError, "Make sure to specify a monotone increasing \"Positions\" array");
            }

            std::array<ScalarVector, dim> globalPositions;
            using std::pow;
            for (int dimIdx = 0; dimIdx < dim; dimIdx++)
            {
                std::size_t numGlobalPositions = 0;
                for (int zoneIdx = 0; zoneIdx < cells[dimIdx].size(); ++zoneIdx)
                    numGlobalPositions += cells[dimIdx][zoneIdx] + 1;

                globalPositions[dimIdx].resize(numGlobalPositions);
                std::size_t posIdx = 0;
                for (int zoneIdx = 0; zoneIdx < cells[dimIdx].size(); ++zoneIdx)
                {
                    const Scalar lower = positions[dimIdx][zoneIdx];
                    const Scalar upper = positions[dimIdx][zoneIdx+1];
                    const int numCells = cells[dimIdx][zoneIdx];
                    Scalar gradingFactor = grading[dimIdx][zoneIdx];
                    const Scalar length = upper - lower;
                    Scalar height = 1.0;
                    bool increasingCellSize = false;

                    if (verbose)
                    {
                        std::cout << "dim " << dimIdx
                                  << " lower "  << lower
                                  << " upper "  << upper
                                  << " numCells "  << numCells
                                  << " grading "  << gradingFactor;
                    }

                    if (gradingFactor > 1.0)
                        increasingCellSize = true;

                    // take absolute values and reverse cell size increment to achieve
                    // reverse behavior for negative values
                    if (gradingFactor < 0.0)
                    {
                        using std::abs;
                        gradingFactor = abs(gradingFactor);

                        if (gradingFactor < 1.0)
                            increasingCellSize = true;
                    }

                    // if the grading factor is exactly 1.0 do equal spacing
                    if (gradingFactor > 1.0 - 1e-7 && gradingFactor < 1.0 + 1e-7)
                    {
                        height = 1.0 / numCells;
                        if (verbose)
                            std::cout << " -> h "  << height * length << std::endl;
                    }

                    // if grading factor is not 1.0, do power law spacing
                    else
                    {
                        height = (1.0 - gradingFactor) / (1.0 - pow(gradingFactor, numCells));

                        if (verbose)
                        {
                            std::cout << " -> grading_eff "  << gradingFactor
                                      << " h_min "  << height * pow(gradingFactor, 0) * length
                                      << " h_max "  << height * pow(gradingFactor, numCells-1) * length
                                      << std::endl;
                        }
                    }

                    // the positions inside a global segment
                    Dune::DynamicVector<Scalar> localPositions(numCells, 0.0);
                    for (int i = 1; i < numCells; i++)
                    {
                        Scalar hI = height;
                        if (!(gradingFactor < 1.0 + 1e-7 && gradingFactor > 1.0 - 1e-7))
                        {
                            if (increasingCellSize)
                                hI *= pow(gradingFactor, i-1);

                            else
                                hI *= pow(gradingFactor, numCells-i);
                        }
                        localPositions[i] = localPositions[i-1] + hI;
                    }

                    localPositions *= length;
                    localPositions += lower;

                    for (int i = 0; i < numCells; ++i)
                        globalPositions[dimIdx][posIdx++] = localPositions[i];
                }

                globalPositions[dimIdx][posIdx] = positions[dimIdx].back();
            }

            polarCoordinates[0] = globalPositions[indices[0]];
            polarCoordinates[1] = globalPositions[indices[1]];
            if (dim == 3)
                polarCoordinates[2] = globalPositions[dim - indices[0] - indices[1]];

            // convert angular coordinates into radians
            std::transform(polarCoordinates[1].begin(), polarCoordinates[1].end(),
                           polarCoordinates[1].begin(),
                           [](Scalar s){ return s*M_PI/180; });

        }
        catch (Dumux::ParameterException &e)
        {
            DUNE_THROW(Dumux::ParameterException, "Please supply all mandatory parameters:" << std::endl
                                                  << GET_PROP_VALUE(TypeTag, GridParameterGroup)
                                                  << ".Positions0, ..." << std::endl);
        }
        catch (...) { throw; }
    }

     /*!
     * \brief Creates cartesian grid from polar coordinates.
     *
     * \param polarCoordinates Vector containing radial, angular and axial coordinates (in this order)
     * \param indices Indices specifing the radial, angular and axial direction (in this order)
     */
    static std::shared_ptr<Grid> createCakeGrid(std::array<ScalarVector, dim> &polarCoordinates,
                                                Dune::FieldVector<int, dim> &indices,
                                                bool verbose = false)
    {
        ScalarVector dR = polarCoordinates[0];
        ScalarVector dA = polarCoordinates[1];

        GridFactory gridFactory;
        Dune::GeometryType type; type.makeCube(dim);

        // create nodes
        if (dim == 3)
        {
            ScalarVector dZ = polarCoordinates[2];
            for (int j = 0; j <= dA.size() - 1; ++j)
            {
                for (int l = 0; l <= dZ.size() - 1; ++l)
                {
                    for (int i = 0; i <= dR.size()- 1; ++i)
                    {
                        // Get radius for the well (= a hole) in the center
                        const Scalar wellRadius = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, WellRadius);

                        // transform into cartesian Coordinates
                        using std::cos;
                        using std::sin;
                        Dune::FieldVector <double, dim> v(0.0);
                        v[indices[2]] = dZ[l];
                        v[indices[0]] = cos(dA[j])*wellRadius + cos(dA[j])*dR[i];
                        v[indices[1]] = sin(dA[j])*wellRadius + sin(dA[j])*dR[i];
                        if (verbose)
                        {
                            std::cout << "Coordinates of : "
                                      << v[0] << " " << v[1] << " " << v[2] << std::endl;
                        }
                        gridFactory.insertVertex(v);
                    }
                }
            }

            if (verbose)
                std::cout << "Filled node vector" << std::endl;

            // assign nodes
            unsigned int z = 0;
            unsigned int t = 0;
            for (int j = 0; j < dA.size() - 1; ++j)
            {
                for (int l = 0; l < dZ.size() - 1; ++l)
                {
                    if (j < dA.size() - 2)
                    {
                        for (int i = 0; i < dR.size() - 1; ++i)
                        {
                            unsigned int rSize = dR.size();
                            unsigned int zSize = dZ.size();
                            std::vector<unsigned int> vid({z, z+1, z+rSize*zSize,
                                                           z+rSize*zSize+1, z+rSize, z+rSize+1,
                                                           z+rSize*zSize+rSize, z+rSize*zSize+rSize+1});

                            if (verbose)
                            {
                                std::cout << "element vertex ids: ";
                                for (int k = 0; k < vid.size(); ++k)
                                    std::cout << vid[k] << " ";
                                std::cout << std::endl;
                            }

                            gridFactory.insertElement(type, vid);

                            z = z+1;
                        }
                        z = z+1;
                    }
                    else
                    {
                        // assign nodes for 360°-cake
                        for (int i = 0; i < dR.size() - 1; ++i)
                        {
                            // z = z + 1;
                            unsigned int rSize = dR.size();
                            std::vector<unsigned int> vid({z, z+1, t,
                                                           t+1, z+rSize, z+rSize+1,
                                                           t+rSize, t+rSize+1});

                            if (verbose)
                            {
                                std::cout << "element vertex ids: ";
                                for (int k = 0; k < vid.size(); ++k)
                                    std::cout << vid[k] << " ";
                                std::cout << std::endl;
                            }

                            gridFactory.insertElement(type, vid);
                            t = t + 1;
                            z = z+1;
                        }
                        t = t + 1;
                        z = z+1;
                    }
                    if (verbose)
                        std::cout << "assign nodes 360° ends..." << std::endl;
                }

                z = z + dR.size();
            }
        }

        // for dim = 2
        else
        {
            for (int j = 0; j <= dA.size() - 1; ++j)
            {
                for (int i = 0; i <= dR.size()- 1; ++i)
                {
                    // Get radius for the well (= a hole) in the center
                    const Scalar wellRadius = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, WellRadius);

                    // transform into cartesian Coordinates
                    Dune::FieldVector <double, dim> v(0.0);

                    v[indices[0]] = cos(dA[j])*wellRadius + cos(dA[j])*dR[i];
                    v[indices[1]] = sin(dA[j])*wellRadius + sin(dA[j])*dR[i];
                    if(verbose) std::cout << "Coordinates of : " << v[0] << " " << v[1] << std::endl;
                    gridFactory.insertVertex(v);
                }
            }
            std::cout << "Filled node vector" << std::endl;

            // assign nodes
            unsigned int z = 0;
            unsigned int t = 0;
            for (int j = 0; j < dA.size() - 1; ++j)
            {
                if (j < dA.size() - 2)
                {
                    for (int i = 0; i < dR.size() - 1; ++i)
                    {
                        unsigned int rSize = dR.size();
                        std::vector<unsigned int> vid({z, z+1, z+rSize, z+rSize+1});

                        if (verbose)
                        {
                            std::cout << "element vertex ids: ";
                            for (int k = 0; k < vid.size(); ++k)
                                std::cout << vid[k] << " ";
                            std::cout << std::endl;
                        }

                        gridFactory.insertElement(type, vid);
                        z = z+1;
                    }
                    z = z+1;
                }
                else
                {
                    // assign nodes for 360°-cake
                    for (int i = 0; i < dR.size() - 1; ++i)
                    {
                        // z = z + 1;
                        std::vector<unsigned int> vid({z, z+1, t, t+1});

                        if (verbose)
                        {
                            std::cout << "element vertex ids: ";
                            for (int k = 0; k < vid.size(); ++k)
                                std::cout << vid[k] << " ";
                            std::cout << std::endl;
                        }

                        gridFactory.insertElement(type, vid);
                        t = t + 1;
                        z = z+1;
                    }
                    t = t + 1;
                    z = z+1;
                }
                std::cout << "assign nodes 360 ends..." << std::endl;
            }
        }
        // return the grid pointer
        return std::shared_ptr<Grid>(gridFactory.createGrid());

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
        static GridPointer cakeGrid;
        return cakeGrid;
    }
};

} // end namespace Dumux

#endif
