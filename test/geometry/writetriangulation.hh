//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Common
 * \brief A function to write a triangulation to vtk
 */
#ifndef DUMUX_TEST_WRITE_TRIANGULATION_HH
#define DUMUX_TEST_WRITE_TRIANGULATION_HH

#include <fstream>
#include <string>
#include <iterator>

namespace Dumux {

//! TriangleVector has to be a nested vector of triangles in 3d
template<class TriangleVector>
void writeVTKPolyDataTriangle(const TriangleVector& triangles,
                              const std::string& filename)
{
    std::ofstream fout(filename + ".vtp");
    fout << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PolyData\" version=\"0.1\">\n"
         << "  <PolyData>\n"
         << "    <Piece NumberOfPoints=\"" << triangles.size()*3 << "\" NumberOfLines=\"" << triangles.size()*3 << "\">\n"
         << "      <Points>\n"
         << "        <DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto& t : triangles)
    {
        for (const auto& p : t)
        {
            fout << p << " ";
            if (p.size() == 1)
                fout << "0.0 0.0 ";
            else if (p.size() == 2)
                fout << "0.0 ";
        }

    }
    fout << '\n';

    fout << "        </DataArray>\n"
         << "      </Points>\n"
         << "      <Lines>\n"
         << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    int offset = 0;
    for (const auto& t : triangles)
    {
        for (std::size_t i = 0; i < t.size()-1; ++i)
            fout << i + offset*3 << " " << i+1 + offset*3 << " ";
        fout << t.size()-1 + offset*3 << " " << offset*3 << " ";
        ++offset;
    }

    fout << "        </DataArray>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    offset = 0;
    for (const auto& t : triangles)
    {
        for (std::size_t i = 1; i <= t.size(); ++i)
            fout << i*2  + offset*6 << " ";
        ++offset;
    }
    fout << '\n';

    fout << "        </DataArray>\n"
         << "      </Lines>\n"
         << "    </Piece>\n"
         << "</PolyData>\n"
         << "</VTKFile>\n";
}

//! TetrahedronVector has to be a nested vector of tets in 3d
template<class TetrahedronVector>
void writeVTUTetrahedron(const TetrahedronVector& tets,
                         const std::string& filename)
{
    std::ofstream fout(filename + ".vtu");
    fout << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
         << "  <UnstructuredGrid>\n"
         << "    <Piece NumberOfPoints=\"" << tets.size()*4 << "\" NumberOfCells=\"" << tets.size() << "\">\n"
         << "      <Points>\n"
         << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto& t : tets)
    {
        for (const auto& p : t)
            fout << p << " ";

        fout << '\n';
    }

    fout << "        </DataArray>\n"
         << "      </Points>\n"
         << "      <Cells>\n"
         << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    int offset = 0;
    for (const auto& t : tets)
    {
        fout << offset << " " << offset+1 << " " << offset+2 << " " << offset+3  << "\n";
        offset += t.size();
    }

    fout << "        </DataArray>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    offset = 4;
    for (const auto& t : tets)
    {
        fout << offset << '\n';
        offset += t.size();
    }

    fout << "        </DataArray>\n";
    fout << "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    for (int i = 0; i < tets.size(); ++i)
        fout << "10\n";

    fout << "        </DataArray>\n"
         << "      </Cells>\n"
         << "    </Piece>\n"
         << "</UnstructuredGrid>\n"
         << "</VTKFile>\n";
}

template<class LineVector>
void writeVTKPolyDataLines(const LineVector& lines,
                           const std::string& filename)
{
    std::ofstream fout(filename + ".vtp");
    fout << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PolyData\" version=\"0.1\">\n"
         << "  <PolyData>\n"
         << "    <Piece NumberOfPoints=\"" << lines.size()*2 << "\" NumberOfLines=\"" << lines.size() << "\">\n"
         << "      <Points>\n"
         << "        <DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto& l : lines)
    {
        for (const auto& p : l)
        {
            fout << p << " ";
            if (p.size() == 1)
                fout << "0.0 0.0 ";
            else if (p.size() == 2)
                fout << "0.0 ";
        }
    }
    fout << '\n';

    fout << "        </DataArray>\n"
         << "      </Points>\n"
         << "      <Lines>\n"
         << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    int offset = 0;
    for (int i = 0; i < lines.size(); ++i)
    {
        fout << offset*2 << " " << offset*2 + 1 << "\n";
        ++offset;
    }

    fout << "        </DataArray>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    for (int i = 0; i < lines.size(); ++i)
        fout << (i+1)*2 << "\n";

    fout << "        </DataArray>\n"
         << "      </Lines>\n"
         << "    </Piece>\n"
         << "</PolyData>\n"
         << "</VTKFile>\n";
}

} // end namespace Dumux

#endif
