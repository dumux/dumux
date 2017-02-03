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
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 */
#ifndef STAGGERED_VTK_WRITER_HH
#define STAGGGERED_VTK_WRITER_HH

#include <dune/common/fvector.hh>

#include <dumux/io/vtkoutputmodulebase.hh>

namespace Properties
{
NEW_PROP_TAG(VtkAddVelocity);
NEW_PROP_TAG(VtkAddProcessRank);
}

namespace Dumux
{

/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 *
 * Handles the output of scalar and vector fields to VTK formatted file for multiple
 * variables and timesteps. Certain predefined fields can be registered on problem / model
 * initialization and/or be turned on/off using the designated properties. Additionally
 * non-standardized scalar and vector fields can be added to the writer manually.
 */
template<typename TypeTag>
class StaggeredVtkWriter
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    static constexpr unsigned int precision = 6;
    static constexpr unsigned int numBeforeLineBreak = 15;

public:
    StaggeredVtkWriter(const Problem& problem) : problem_(problem)
    {}


    template<class T>
    void write(const T& info)
    {
        file_.open(fileName_());
        write_(info);
        file_.close();
        ++curWriterNum_;
    }


private:

    void writeHeader_(const int numPoints)
    {
        std::string header = "<?xml version=\"1.0\"?>\n";
                    header += "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
                    header += "<PolyData>\n";
                    header += "<Piece NumberOfLines=\"0\" NumberOfPoints=\"" + std::to_string(numPoints) + "\">\n";
        file_ << header;
    }

    template<class T>
    void write_(const T& info)
    {
        const int numPoints = problem_.model().numFaceDofs();
        writeHeader_(numPoints);


        // get fields for all primary variables
        std::vector<Scalar> positions(numPoints*3);
        std::vector<std::vector<Scalar>> priVarScalarData(info.size(), std::vector<Scalar>(numPoints));

        for(auto&& facet : facets(problem_.gridView()))
        {
            const int dofIdxGlobal = problem_.gridView().indexSet().index(facet);

            // first, get the coordinates
            auto pos = facet.geometry().center();

            if(dim == 1)
            {
                positions[dofIdxGlobal*3] = pos[0];
                positions[dofIdxGlobal*3 + 1] = 0.0;
                positions[dofIdxGlobal*3 + 2] = 0.0;
            }
            else if(dim == 2)
            {
                positions[dofIdxGlobal*3] = pos[0];
                positions[dofIdxGlobal*3 + 1] = pos[1];
                positions[dofIdxGlobal*3 + 2] = 0.0;
            }
            else
            {
                positions[dofIdxGlobal*3] = pos[0];
                positions[dofIdxGlobal*3 + 1] = pos[1] ;
                positions[dofIdxGlobal*3 + 2] = pos[2] ;
            }

            // now, get the actual values
            for(int pvIdx = 0; pvIdx < info.size(); ++pvIdx)
            {
                priVarScalarData[pvIdx][dofIdxGlobal] = problem_.model().curSol()[faceIdx][dofIdxGlobal][pvIdx];
            }
        }

        writeCoordinates_(positions);
        writeFacePriVars_(priVarScalarData, info);

        file_ << "</Piece>\n";
        file_ << "</PolyData>\n";
        file_ << "</VTKFile>";
    }

    void writeCoordinates_(const std::vector<Scalar>& positions)
    {
        // write the positions to the file
        file_ << "<Points>\n";
        file_ << "<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";

        int counter = 0;
        for(auto&& x : positions)
        {
            file_ << x << " ";
            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file_ << std::endl;
                counter = 0;
            }
        }

        file_ << "\n</DataArray>\n";
        file_ << "</Points>\n";
    }

    template<class T>
    void writeFacePriVars_(const std::vector<std::vector<Scalar>>& priVarScalarData, const T& info)
    {
        // write the priVars to the file
        file_ << "<PointData Scalars=\"" << info[0].name << "\">";

        for(int pvIdx = 0; pvIdx < info.size(); ++pvIdx)
        {
            file_ << "<DataArray type=\"Float32\" Name=\"" << info[pvIdx].name << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";

            int counter = 0;
            for(auto&& x : priVarScalarData[pvIdx])
            {
                file_ << x << " " ;
                // introduce a line break after a certain time
                if((++counter)  > numBeforeLineBreak)
                {
                    file_ << std::endl;
                    counter = 0;
                }
            }

            file_ << "\n</DataArray>\n";
        }

        file_ << "</PointData>\n";
    }

    std::string fileName_() const
    {
        std::ostringstream oss;
        oss << problem_.name() << "-face-" << std::setw(5) << std::setfill('0') << curWriterNum_ << ".vtp";
        return oss.str();
    }

    const Problem& problem_;
    std::ofstream file_;
    int curWriterNum_{0};
};
} // end namespace Dumux

#endif
