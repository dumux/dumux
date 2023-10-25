// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Draw Grid function for face-centered staggered. Modified from Dune-Grid's printgrid method
 */

#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_DRAWGRID_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_DRAWGRID_HH

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dumux/discretization/facecentered/staggered/normalaxis.hh>

namespace Dumux {

namespace {

// Shift a point closer to basegeo's center by factor scale (used for drawing relative to the element)
template<class GlobalPosition>
GlobalPosition shift(const GlobalPosition& centerCoords, const GlobalPosition& coords, const double shiftScale, const int direction)
{
    double shift = (coords[direction] - centerCoords[direction]) * shiftScale;
    GlobalPosition ret = coords;
    ret[direction] += shift;
    return ret;
}

// Scale a point closer to basegeo's center by factor scale (used for drawing relative to the element)
template<class GlobalPosition>
GlobalPosition centrify(const GlobalPosition& centerCoords, const GlobalPosition& coords, const double scale)
{
    GlobalPosition ret = coords;
    ret -= centerCoords;
    ret *= scale;
    ret += centerCoords;
    return ret;
}

// Add a line to the plotfile from p1 to p2
template<class GlobalPosition>
void drawLine (std::ofstream &plotfile, const GlobalPosition &p1, const GlobalPosition &p2, const std::string& options) {
    plotfile << "set object poly from ";
    plotfile << p1[0] << "," << p1[1] << " to ";
    plotfile << p2[0] << "," << p2[1] << " to ";
    plotfile << p1[0] << "," << p1[1];
    plotfile << " " << options << std::endl;
}

template<class GlobalPosition>
void drawCircle(std::ofstream& plotfile, const GlobalPosition& center, const double size, std::string options = "fc rgb \"navy\"")
{
    plotfile << "set object circle at ";
    plotfile << center[0] << "," << center[1];
    plotfile << " size " << size;
    plotfile << " " << options << std::endl;
}

template<class GlobalPosition>
void drawArrow(std::ofstream& plotfile, const GlobalPosition& center, const double elementRadius, const double scale, const int direction)
{
    const double size = elementRadius * scale;
    GlobalPosition start = center;
    if (direction == 0)
        start[0] -= size/2;
    else
        start[1] -= size/2;

    const std::string shift = direction ? " 0.0," + std::to_string(size) : std::to_string(size) + ",0.0";
    const std::string options = direction ? "lw 1 lc rgb \"red\"" : "lw 1 lc rgb \"blue\"";

    plotfile << "set arrow from ";
    plotfile << start[0] << "," << start[1];
    plotfile << " rto " << shift;
    plotfile << " " << options << std::endl;
}

} // end anonymous namespace

/** \brief Print a grid geometry as a gnuplot for testing and development
*  \tparam GridGeometry the grid geometry used
*  \param gridGeometry the grid geometry to print
*  \param outputFileName the base of the output filename
*  \param size size of the plot in pixels; increase if plot is too cramped
*  \param executePlot whether to execute gnuplot automatically
*  \param png whether to use PNG or SVG as the output format
*  Creates a gnuplot (one per process if parallel) showing the grid structure with indices, intersection types etc.
*/
template<class GridGeometry>
void drawGridGeometry(const GridGeometry& gridGeometry,
                    const std::string& outputFileName = "printgrid",
                    int size = 2000,
                    bool executePlot = true,
                    bool png = true)
{
    static_assert(GridGeometry::Grid::dimension == 2, "drawGridGeometry only works for 2D grids");

    // Create output file
    std::string plotFileName = outputFileName + ".gnuplot";
    std::ofstream plotFile (plotFileName, std::ios::out | std::ios::trunc);

    if (!plotFile.is_open())
    {
        DUNE_THROW(Dune::IOError, "Could not create plot file " << outputFileName << "!");
        return;
    }

    // Basic plot output settings
    plotFile << "set size ratio -1" << std::endl;
    if (png)
    {
        plotFile << "set terminal png size " << size << "," << size << std::endl;
        plotFile << "set output '" << outputFileName << ".png'" << std::endl;
    }
    else
    {
        plotFile << "set terminal svg size " << size << "," << size << " enhanced background rgb 'white'" << std::endl;
        plotFile << "set output '" << outputFileName << ".svg'" << std::endl;
    }

    plotFile << "set style textbox opaque noborder" << std::endl;

    // Get GridView
    const auto gridView = gridGeometry.gridView();

    // Will contain min/max coordinates. Needed for scaling of the plot
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    GlobalPosition maxCoord;
    GlobalPosition minCoord;

    // Iterate over elements
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        fvGeometry.bind(element);

        // Plot element index at element center
        const GlobalPosition elementCenter = element.geometry().center();
        plotFile << "set label at " << elementCenter[0] << "," << elementCenter[1] << " '" << eIdx << "' center" << std::endl;

        // draw element boundaries
        for (const auto& intersection : intersections(gridView, element))
        {
            // draw grid intersections
            auto intersectionGeometry = intersection.geometry();
            auto elementRadius = (elementCenter - intersectionGeometry.center()).two_norm();
            drawLine(plotFile, intersectionGeometry.corner(0), intersectionGeometry.corner(1), "fs empty border 1");

            // draw arrows at dofs
            GlobalPosition shiftedIntersectionCenter = intersectionGeometry.center();
            if (Dumux::normalAxis(intersection.centerUnitOuterNormal()) == 0)
                shiftedIntersectionCenter[1] -= 0.05;
            else
                shiftedIntersectionCenter[0] -= 0.05;
            drawArrow(plotFile, shiftedIntersectionCenter, elementRadius, 0.15, Dumux::normalAxis(intersection.centerUnitOuterNormal()));

            // draw an element Box around element index
            // Thin line with color according to neighbor()
            auto innerCorner1 = centrify(elementCenter, intersectionGeometry.corner(0), 0.1);
            auto innerCorner2 = centrify(elementCenter, intersectionGeometry.corner(1), 0.1);

            if (intersection.neighbor())
                drawLine(plotFile, innerCorner1, innerCorner2, "fs empty border 2");
            else
                drawLine(plotFile, innerCorner1, innerCorner2, "fs empty border 1");

            // Thick line in case of boundary()
            if (intersection.boundary())
                drawLine(plotFile, innerCorner1, innerCorner2, "fs empty border 3 lw 4");

            // set max and min grid coordinates
            for (int i = 0; i < intersectionGeometry.corners(); ++i)
            {
                // Adapt min / max coordinates
                for (int dim = 0; dim < 2; ++dim)
                {
                    if (intersectionGeometry.corner(i)[dim] < minCoord[dim])
                        minCoord[dim] = intersectionGeometry.corner(i)[dim];
                    else if (intersectionGeometry.corner(i)[dim] > maxCoord[dim])
                        maxCoord[dim] = intersectionGeometry.corner(i)[dim];
                }
            }

            // draw shifted frontal scv boundaries
            int dirIdx = Dumux::normalAxis(intersection.centerUnitOuterNormal());
            std::array<GlobalPosition, 4> pseudoScvCorners;
            pseudoScvCorners[0] = shift(elementCenter,
                                        centrify(elementCenter, intersectionGeometry.corner(0), 0.45),
                                        0.15, dirIdx);
            pseudoScvCorners[1] = shift(elementCenter,
                                        centrify(elementCenter, intersectionGeometry.corner(1), 0.45),
                                        0.15, dirIdx);
            pseudoScvCorners[2] = shift(elementCenter,
                                        centrify(elementCenter, intersectionGeometry.corner(1), 0.45),
                                        1.20, dirIdx);
            pseudoScvCorners[3] = shift(elementCenter,
                                        centrify(elementCenter, intersectionGeometry.corner(0), 0.45),
                                        1.20, dirIdx);

            for (int i = 0; i < pseudoScvCorners.size(); i++)
            {
                if (dirIdx == 0)
                {
                    if(i == 3)
                        drawLine(plotFile, pseudoScvCorners[3], pseudoScvCorners[0], "fs empty border 6");
                    else
                        drawLine(plotFile, pseudoScvCorners[i], pseudoScvCorners[i+1], "fs empty border 6");
                }
                else
                {
                    if(i == 3)
                        drawLine(plotFile, pseudoScvCorners[3], pseudoScvCorners[0], "fs empty border 7");
                    else
                        drawLine(plotFile, pseudoScvCorners[i], pseudoScvCorners[i+1], "fs empty border 7");
                }
            }
        }

        // Iterate through the scvs (SCVs and SCVFs only have center positions stored?)
        for (auto&& scv : scvs(fvGeometry))
        {
            // Plot scv index at the scv center
            // shift scv center in it's direction
            GlobalPosition shiftedScvCenter = shift(elementCenter, scv.center(), 0.5, scv.dofAxis());
            plotFile << "set label at " << shiftedScvCenter[0] << "," << shiftedScvCenter[1] << " '" << scv.index() << "' center" << std::endl;
            drawCircle(plotFile, shiftedScvCenter, 0.05, "fs empty border 0");

            // Plot the Dof index adjacent to the scv center
            const std::string suffix = scv.dofAxis() == 0 ? "' boxed center" : "' boxed center rotate by 90";
            const double shiftFactor = scv.boundary() ? 1.2 : 1.0;
            const GlobalPosition shiftedDofLocation = shift(elementCenter, scv.center(), shiftFactor, scv.dofAxis());
            plotFile << "set label at " << shiftedDofLocation[0] << "," << shiftedDofLocation[1] << " 'dof " << int(scv.dofIndex()) << suffix << std::endl;

            // Iterate through each scvf in each scv
            for (auto&& scvf : scvfs(fvGeometry, scv))
            {
                //scale the scvf centers
                auto shiftedScvfCenter = centrify(scv.center(), scvf.center(), 0.25);
                shiftedScvfCenter = shift(elementCenter, shiftedScvfCenter, 0.5, scv.dofAxis());

                // Plot scvf index at the scvf center, (offset)
                if (scvf.boundary())
                    plotFile << "set label at " << shiftedScvfCenter[0] << "," << shiftedScvfCenter[1] << " '" << scvf.index() << "(b)' center" << std::endl;
                else
                    plotFile << "set label at " << shiftedScvfCenter[0] << "," << shiftedScvfCenter[1] << " '" << scvf.index() << "' center" << std::endl;
            }
        }
    }

    // Finish plot, pass extend of the grid
    GlobalPosition extend(maxCoord - minCoord);

    extend *= 0.075;
    minCoord -= extend;
    maxCoord += extend;
    plotFile << "plot [" << minCoord[0] << ":" << maxCoord[0] << "] [" << minCoord[1]
                << ":" << maxCoord[1] << "] NaN notitle" << std::endl;
    plotFile.close();

#if DUMUX_HAVE_GNUPLOT
    if (executePlot)
    {
        std::string cmd = "gnuplot -p '" + plotFileName + "'";
        if (std::system (cmd.c_str()) != 0)
            DUNE_THROW(Dune::Exception, "Error running GNUPlot: " << cmd);
    }
#else
    std::cout << "No gnuplot found. Cannot visualize result. \n";
#endif

}

} // end namespace Dumux

#endif
