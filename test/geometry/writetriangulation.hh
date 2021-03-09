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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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
         << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
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

template<class LineVector>
void writeVTKPolyDataLines(const LineVector& lines,
                           const std::string& filename)
{
    std::ofstream fout(filename + ".vtp");
    fout << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
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
