// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup InputOutput
 * \brief Provides a grid manager for a 3D landfill grid
 */
#ifndef DUMUX_LANDFILL_GRID_MANAGER_HH
#define DUMUX_LANDFILL_GRID_MANAGER_HH

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <dune/common/dynvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/math.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for creating a three-dimensional
 *        grid for a landfill.
 */
template <class Grid>
class LandfillGridManager
{
    using Scalar = typename Grid::ctype;

    using GridFactory = Dune::GridFactory<Grid>;
    using GridPointer = std::shared_ptr<Grid>;

    enum { dim = Grid::dimension,
           dimWorld = Grid::dimensionworld };

public:
    /*!
     * \brief Make the grid.
     */
    void init(const std::string& modelParamGroup = "")
    {
        static_assert(dim == 3, "The Landfill Gridmanager is only implemented for 3D.");

        const bool verbose = getParamFromGroup<bool>(modelParamGroup, "Grid.Verbosity", false);

        gridPtr() = createLandfillGrid(modelParamGroup, verbose);

        loadBalance();
    }

    /*!
     * \brief Creates a grid consisting a landfill and the underground below.
     *
     * \param modelParamGroup name of the model parameter group
     * \param verbose if the output should be verbose
     */
    std::unique_ptr<Grid> createLandfillGrid(const std::string& modelParamGroup,
                                           bool verbose = false)
    {
        const auto sizeXLF = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.SizeXLF");
        const auto sizeYLF = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.SizeYLF");
        const auto sizeZLF = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.SizeZLF");
        const auto cellsLF = getParamFromGroup<std::vector<unsigned>>(modelParamGroup, "Grid.CellsLandfill");

        const auto sizeXAQ = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.SizeXAQ");
        const auto sizeYAQ = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.SizeYAQ");
        const auto sizeZAQ = getParamFromGroup<std::vector<Scalar>>(modelParamGroup, "Grid.SizeZAQ");
        const auto cellsAQ = getParamFromGroup<std::vector<unsigned>>(modelParamGroup, "Grid.CellsAquifer");

        //auto slopeLF = getParamFromGroup<Scalar>(modelParamGroup, "Grid.SlopeLF", 0.0);

        // grid size in each coordinate direction
        const auto hXLF = (sizeXLF[1] - sizeXLF[0])/cellsLF[0];
        const auto hYLF = (sizeYLF[1] - sizeYLF[0])/cellsLF[1];
        const auto hZLF = (sizeZLF[1] - sizeZLF[0])/cellsLF[2];

        // grid size in each coordinate direction
        const auto hXAQ = (sizeXAQ[1] - sizeXAQ[0])/cellsAQ[0];
        const auto hYAQ = (sizeYAQ[1] - sizeYAQ[0])/cellsAQ[1];
        const auto hZAQ = (sizeZAQ[1] - sizeZAQ[0])/cellsAQ[2];

        GridFactory gridFactory;
        const auto noPts = (cellsAQ[0]+1)*(cellsAQ[1]+1)*(cellsAQ[2]+cellsLF[2]);
        // vector to hold all points
        std::vector<Dune::FieldVector <double, dim>> vertex(noPts, Dune::FieldVector <double, dim>(0.0));

        // Landfill:
        // create all landfill points
        auto n = 0u;
        for (auto k = 0u; k <= cellsLF[2]; ++k)
        {
            for (auto j = 0u; j <= cellsLF[1]+k*2; ++j)
            {
                for (auto i = 0u; i <= cellsLF[0]+k*2; ++i)
                {
                    Dune::FieldVector <double, dim> v(0.0);

                    // determine coordinates
                    v[0] = sizeXLF[0] + i*hXLF - k*hXLF;
                    v[1] = sizeYLF[0] + j*hYLF - k*hYLF;
                    v[2] = sizeZLF[1] - k*hZLF; //k*(hZLF+(1-v[0])/(sizeXLF[1]-sizeXLF[0]));
                    vertex[n] = v;

                    if(verbose)
                        printCoordinate(n, v);
                    gridFactory.insertVertex(v);

                    ++n;
                }
            }
        }

        auto noPtsLF = 0u;
        for (auto k = 0u; k < cellsLF[2]; ++k)
        {
            const auto noPtsX = cellsLF[0]+1+k*2;
            const auto noPtsX2 = cellsLF[0]+1+(k+1)*2;
            const auto noPtsY = cellsLF[1]+1+k*2;

            for (auto j = 0u; j < noPtsY-1; ++j)
            {
                for (auto i = 0u; i < noPtsX-1; ++i)
                {
                    const std::vector<unsigned> vid({i+j*noPtsX+noPtsLF, i+j*noPtsX+noPtsLF+1,
                                                    i+(j+1)*noPtsX+noPtsLF, i+(j+1)*noPtsX+noPtsLF+1,
                                                    i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+1, i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+2,
                                                    i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+1, i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+2});
                    if (verbose)
                        printIndices(vid);
                    constexpr auto type = Dune::GeometryTypes::cube(dim);
                    gridFactory.insertElement(type, vid);

                    if(i == 0)
                    {
                        const std::vector<unsigned> vid({i+j*noPtsX+noPtsLF,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+1,
                                                        i+(j+1)*noPtsX+noPtsLF,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+1
                                                        });
                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::prism;
                        gridFactory.insertElement(type, vid);
                    }
                    if(i == noPtsX-2)
                    {
                        const std::vector<unsigned> vid({i+j*noPtsX+noPtsLF+1,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+2,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+3,
                                                        i+(j+1)*noPtsX+noPtsLF+1,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+2,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+3
                                                        });
                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::prism;
                        gridFactory.insertElement(type, vid);

                    }

                    if(j == 0)
                    {
                        const std::vector<unsigned> vid({i+j*noPtsX+noPtsLF,
                                                        i+noPtsX*noPtsY+noPtsLF+1,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+1,
                                                        i+j*noPtsX+noPtsLF+1,
                                                        i+noPtsX*noPtsY+noPtsLF+2,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+2
                                                        });
                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::prism;
                        gridFactory.insertElement(type, vid);

                    }
                    if(j == noPtsY-2)
                    {
                        const std::vector<unsigned> vid({i+(j+1)*noPtsX+noPtsLF,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+1,
                                                        i+(j+3)*noPtsX2+noPtsX*noPtsY+noPtsLF+1,
                                                        i+(j+1)*noPtsX+noPtsLF+1,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+2,
                                                        i+(j+3)*noPtsX2+noPtsX*noPtsY+noPtsLF+2
                                                        });
                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::prism;
                        gridFactory.insertElement(type, vid);

                    }
                    if(i == 0 && j == 0)
                    {
                        const std::vector<unsigned> vid({i+j*noPtsX+noPtsLF,
                                                        i+j*noPtsX2+noPtsX*noPtsY+noPtsLF+1,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF,
                                                        i+j*noPtsX2+noPtsX*noPtsY+noPtsLF,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+1
                                                        });
                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::pyramid;
                        gridFactory.insertElement(type, vid);
                    }
                    if(i == noPtsX-2 && j == 0)
                    {
                        const std::vector<unsigned> vid({i+j*noPtsX+noPtsLF+1,
                                                        i+j*noPtsX2+noPtsX*noPtsY+noPtsLF+3,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+2,
                                                        i+j*noPtsX2+noPtsX*noPtsY+noPtsLF+2,
                                                        i+(j+1)*noPtsX2+noPtsX*noPtsY+noPtsLF+3
                                                        });
                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::pyramid;
                        gridFactory.insertElement(type, vid);
                    }
                    if(i == 0 && j == noPtsY-2)
                    {
                        const std::vector<unsigned> vid({i+(j+1)*noPtsX+noPtsLF,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+1,
                                                        i+(j+3)*noPtsX2+noPtsX*noPtsY+noPtsLF,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF,
                                                        i+(j+3)*noPtsX2+noPtsX*noPtsY+noPtsLF+1
                                                        });
                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::pyramid;
                        gridFactory.insertElement(type, vid);

                    }
                    if(i == noPtsX-2 && j == noPtsY-2)
                    {
                        const std::vector<unsigned> vid({i+(j+1)*noPtsX+noPtsLF+1,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+3,
                                                        i+(j+3)*noPtsX2+noPtsX*noPtsY+noPtsLF+2,
                                                        i+(j+2)*noPtsX2+noPtsX*noPtsY+noPtsLF+2,
                                                        i+(j+3)*noPtsX2+noPtsX*noPtsY+noPtsLF+3
                                                        });
                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::pyramid;
                        gridFactory.insertElement(type, vid);

                    }

                }
            }
            noPtsLF += noPtsX*noPtsY;
        }
        noPtsLF = n;

        // Aquifer:
        // create all aquifer points
        for (auto k = 0u; k <= cellsAQ[2]; ++k)
        {
            for (auto j = 0u; j <= cellsAQ[1]; ++j)
            {
                for (auto i = 0u; i <= cellsAQ[0]; ++i)
                {
                    Dune::FieldVector <double, dim> v(0.0);

                    // determine coordinates
                    v[0] = sizeXAQ[0] + i*hXAQ;
                    v[1] = sizeYAQ[0] + j*hYAQ;
                    v[2] = sizeZAQ[1] - k*hZAQ;
                    vertex[n] = v;

                    if(verbose)
                        printCoordinate(n, v);
                    gridFactory.insertVertex(v);

                    ++n;
                }
            }
        }

        auto noPtsTotal = noPtsLF;
        for (auto k = 0u; k < cellsAQ[2]; ++k)
        {

            const auto noPtsX = cellsAQ[0]+1;
            const auto noPtsY = cellsAQ[1]+1;
            for (auto j = 0u; j < cellsAQ[1]; ++j)
            {
                for (auto i = 0u; i < cellsAQ[0]; ++i)
                {
                    const std::vector<unsigned> vid({i+j*noPtsX+noPtsTotal, i+j*noPtsX+noPtsTotal+1,
                                                    i+(j+1)*noPtsX+noPtsTotal, i+(j+1)*noPtsX+noPtsTotal+1,
                                                    i+j*noPtsX+noPtsX*noPtsY+noPtsTotal, i+j*noPtsX+noPtsX*noPtsY+noPtsTotal+1,
                                                    i+(j+1)*noPtsX+noPtsX*noPtsY+noPtsTotal, i+(j+1)*noPtsX+noPtsX*noPtsY+noPtsTotal+1});
                    if (verbose)
                        printIndices(vid);
                    constexpr auto type = Dune::GeometryTypes::cube(dim);
                    gridFactory.insertElement(type, vid);

                }
            }
            noPtsTotal += noPtsX*noPtsY;
        }

        // connect the landfill and aquifer with elements
        const auto noPtsXLF = cellsLF[0]+1+2*cellsLF[2];
        const auto noPtsYLF = cellsLF[0]+1+2*cellsLF[2];
        const auto noPtsXAQ = (cellsAQ[0]+1);
        const auto noPtsYAQ = (cellsAQ[1]+1);

        const auto connectionStart = noPtsLF-noPtsXLF*noPtsYLF;

        for (auto j = 0u; j < noPtsYLF-1; ++j)
        {
            for (auto i = connectionStart; i < connectionStart+noPtsXLF-1; ++i)
            {
                for (auto k = noPtsLF; k < noPtsLF+noPtsXAQ*noPtsYAQ; ++k)
                {
                    if(vertex[i+j*noPtsXLF][0] == vertex[k][0] && vertex[i+j*noPtsXLF][1] == vertex[k][1])
                    {
                        const std::vector<unsigned> vid({i+j*noPtsXLF, i+j*noPtsXLF+1,
                                                        i+(j+1)*noPtsXLF, i+(j+1)*noPtsXLF+1,
                                                        k, k+1,
                                                        k+noPtsXAQ, k+noPtsXAQ+1});

                        if (verbose)
                            printIndices(vid);
                        constexpr auto type = Dune::GeometryTypes::cube(dim);
                        gridFactory.insertElement(type, vid);
                    }
                }
            }
        }



        // return the grid pointer
        return std::unique_ptr<Grid>(gridFactory.createGrid());
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    {
        return *gridPtr();
    }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    void loadBalance()
    {
        gridPtr()->loadBalance();
    }

protected:
    static void printCoordinate(const int i, const Dune::FieldVector <double, dim>& v)
    {
        std::cout << "Coordinates of Point index " << i << " : ";
        for (int k = 0; k < v.size(); ++k)
            std::cout << v[k] << " ";
        std::cout << std::endl;
    }

    static void printIndices(const std::vector<unsigned>& vid)
    {
        std::cout << "element vertex indices: ";
        for (int k = 0; k < vid.size(); ++k)
            std::cout << vid[k] << " ";
        std::cout << std::endl;
    }

    /*!
     * \brief Returns a reference to the shared pointer to the grid.
     */
    GridPointer& gridPtr()
    {
        return landfillGrid_;
    }

private:
    GridPointer landfillGrid_;
};

} // end namespace Dumux

#endif
