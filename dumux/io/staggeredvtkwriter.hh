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
 * \brief A VTK writer specialized for staggered grid implementations with dofs on the faces
 */
#ifndef STAGGERED_VTK_WRITER_HH
#define STAGGERED_VTK_WRITER_HH

#include <dune/common/fvector.hh>

#include <dumux/io/vtkoutputmodulebase.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>

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
template<class TypeTag, class WriterData>
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

    using ScalarInfo = std::vector<typename WriterData::PriVarScalarDataInfo>;
    using VectorInfo = std::vector<typename WriterData::PriVarVectorDataInfo>;
    using Positions = typename WriterData::Positions;
    using Data = typename WriterData::Data;

public:
    StaggeredVtkWriter(const Problem& problem) : problem_(problem)
    {}

    void write(const WriterData& data)
    {
        if(data.priVarScalarDataInfo.empty() &&
           data.priVarVectorDataInfo.empty() &&
           data.secondVarScalarDataInfo.empty() &&
           data.secondVarVectorDataInfo.empty())
            return;

        file_.open(fileName_());
        write_(data);
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

     /*!
     * \brief Writes the coordinates and the simulation data to the file
     *
     * \param priVarScalarData Container to store the scalar-valued data
     * \param priVarVectorData Container to store the vector-valued data
     * \param scalarInfo Information concerning the data (e.g. names)
     * \param vectorInfo Information concerning the data (e.g. names)
     * \param positions Container to store the positions
     */
    void write_(const WriterData& data)
    {
        const int numPoints = problem_.model().numFaceDofs();
        writeHeader_(numPoints);
        writeCoordinates_(data.positions);

        std::string scalarName;
        std::string vectorName;

        bool scalarValuesPresent = false;
        bool vectorValuesPresent = false;


        if(!(data.priVarScalarDataInfo.empty() && data.secondVarScalarDataInfo.empty()))
        {
            scalarValuesPresent = true;
            scalarName = data.priVarScalarDataInfo.empty() ? data.secondVarScalarDataInfo[0].name :
                                                                                  data.priVarScalarDataInfo[0].name;
        }
        if(!(data.priVarVectorDataInfo.empty() && data.secondVarVectorDataInfo.empty()))
        {
            vectorValuesPresent = true;
            vectorName = data.priVarVectorDataInfo.empty() ? data.secondVarVectorDataInfo[0].name :
                                                                                  data.priVarVectorDataInfo[0].name;
        }

        if(scalarValuesPresent)
            if(vectorValuesPresent)
                file_ << "<PointData Scalars=\"" << scalarName << "\" Vectors=\"" << vectorName <<"\">";
            else
                file_ << "<PointData Scalars=\"" << scalarName << "\">";
        else if(vectorValuesPresent)
            file_ << "<PointData Vectors=\"" << vectorName << "\">";
        else
            return;

        file_ << std::endl;
        writeAllData_(data);

        file_ << "</PointData>\n";
        file_ << "</Piece>\n";
        file_ << "</PolyData>\n";
        file_ << "</VTKFile>";
    }

     /*!
     * \brief Writes the coordinates to the file
     *
     * \param positions Container to store the positions
     */
    void writeCoordinates_(const Positions& positions)
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

    void writeAllData_(const WriterData& data)
    {
        // write scalar-valued primary variable data
        if(!data.priVarScalarDataInfo.empty())
            writeScalarData_(data.priVarScalarData, data.priVarScalarDataInfo);

        // write vector-valued primary variable data
        if(!data.priVarVectorDataInfo.empty())
            writeVectorData_(data.priVarVectorData, data.priVarVectorDataInfo);

        // write scalar-valued secondary variable data
        if(!data.secondVarScalarDataInfo.empty())
            writeScalarData_(data.secondVarScalarData, data.secondVarScalarDataInfo);

        // write vector-valued primary variable data
        if(!data.secondVarVectorDataInfo.empty())
            writeVectorData_(data.secondVarVectorData, data.secondVarVectorDataInfo);
    }

     /*!
     * \brief Writes scalar-valued data to the file
     *
     * \param priVarScalarData Container to store the data
     * \param info Information concerning the data (e.g. names)
     */
    void writeScalarData_(const Data& priVarScalarData, const ScalarInfo& info)
    {
        // write the priVars to the file
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
    }

     /*!
     * \brief Writes vector-valued data to the file
     *
     * \param priVarVectorData Container to store the data
     * \param info Information concerning the data (e.g. names)
     */
    void writeVectorData_(const Data& priVarVectorData, const VectorInfo& info)
    {
        // write the priVars to the file
        for(int i = 0; i < info.size(); ++i)
        {
            file_ << "<DataArray type=\"Float32\" Name=\"" << info[i].name << "\" NumberOfComponents=\"" << info[i].pvIdx.size() << "\" format=\"ascii\">\n";
            int counter = 0;
            for(auto&& x : priVarVectorData[i])
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
    }

     /*!
     * \brief Returns the file name for each timestep
     */
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
