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
#include <dune/grid/io/file/vtk/common.hh>

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
template<int dim>
class StaggeredVtkWriter
{
    using Scalar = double;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    static constexpr unsigned int precision = 6;
    static constexpr unsigned int numBeforeLineBreak = 15;

    template<class ContainerType>
    class VTKLocalFunction
    {
    public:
        VTKLocalFunction(const ContainerType& data, const std::string& name, const int numComponents) : data_(data), name_(name), numComponents_(numComponents)
        {}

        const std::string& name() const
        {
            return name_;
        }

        int numComponents() const
        {
            return numComponents_;
        }

        auto& operator() (const int dofIdx) const  { return data_[dofIdx]; }

        decltype(auto) begin() const
        {
            return data_.begin();
        }

        decltype(auto) end() const
        {
            return data_.end();
        }

        int size() const
        {
            return data_.size();
        }

    private:
        const ContainerType& data_;
        const std::string name_;
        const int numComponents_;
    };


public:
    using ScalarLocalFunction = VTKLocalFunction<std::vector<Scalar>>;
    using VectorLocalFunction = VTKLocalFunction<std::vector<GlobalPosition>>;


    StaggeredVtkWriter(const std::vector<GlobalPosition>& coordinates) : coordinates_(coordinates)
    {}

    void write(/*const WriterData& data*/)
    {

        file_.open(fileName_());
        writeHeader_();
        writeCoordinates_(coordinates_);
        writeDataInfo_();

        for(auto&& data : scalarPointData_)
        {
            writeData_(data);
        }
        for(auto&& data :vectorPointData_)
        {
            writeData_(data);
        }

        file_ << "</PointData>\n";
        file_ << "</Piece>\n";
        file_ << "</PolyData>\n";
        file_ << "</VTKFile>";

        clear();

        file_.close();
        ++curWriterNum_;
    }


    std::string write ( const std::string &name,
                        Dune::VTK::OutputType type = Dune::VTK::ascii )
    {
        return "dummy";
    }


    void addPointData(const std::vector<Scalar>& v, const std::string &name, int ncomps = 1)
    {
        assert(v.size() == ncomps * coordinates_.size());
        scalarPointData_.push_back(ScalarLocalFunction(v, name, ncomps));
    }

    void addPointData(const std::vector<GlobalPosition>& v, const std::string &name, int ncomps = 1)
    {
        assert(v.size() == coordinates_.size());
        vectorPointData_.push_back(VectorLocalFunction(v, name, ncomps));
    }


    void clear()
    {
        scalarPointData_.clear();
        vectorPointData_.clear();
    }


private:
    void writeHeader_()
    {
        std::string header = "<?xml version=\"1.0\"?>\n";
                    header += "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
                    header += "<PolyData>\n";
                    header += "<Piece NumberOfLines=\"0\" NumberOfPoints=\"" + std::to_string(coordinates_.size()) + "\">\n";
        file_ << header;
    }

     /*!
     * \brief Writes information about the data to the file
     */
    void writeDataInfo_()
    {
        std::string scalarName;
        std::string vectorName;
        bool foundScalar = false;
        bool foundVector = false;

        for(auto&& data : scalarPointData_)
        {
            if(data.numComponents() == 1 && !foundScalar)
            {
                scalarName = data.name();
                foundScalar = true;
                continue;
            }

            if(data.numComponents() > 1 && !foundVector)
            {
                vectorName = data.name();
                foundVector = true;
            }
        }

        for(auto&& data : vectorPointData_)
        {
            if(data.numComponents() > 1 && !foundVector)
            {
                vectorName = data.name();
                foundVector = true;
            }
        }

        if(foundScalar)
            if(foundVector)
                file_ << "<PointData Scalars=\"" << scalarName << "\" Vectors=\"" << vectorName <<"\">\n";
            else
                file_ << "<PointData Scalars=\"" << scalarName << "\">\n";
        else if(foundVector)
            file_ << "<PointData Vectors=\"" << vectorName << "\">\n";
        else
            return;
    }

     /*!
     * \brief Writes the coordinates to the file
     *
     * \param positions Container to store the positions
     */
    void writeCoordinates_(const std::vector<GlobalPosition>& positions)
    {
        // write the positions to the file
        file_ << "<Points>\n";
        file_ << "<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        int counter = 0;
        for(auto&& x : positions)
        {
            file_ << x ;

            if(x.size() == 1)
                file_ << " 0 0 ";
            if(x.size() == 2)
                file_ << " 0 ";

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

     /*!
     * \brief Writes data to the file
     *
     * \param data The data container which hold the data itself, as well as the name of the data set and the number of its components
     */
    template<class T>
    void writeData_(const T& data)
    {
        file_ << "<DataArray type=\"Float32\" Name=\"" << data.name() << "\" NumberOfComponents=\"" << data.numComponents() << "\" format=\"ascii\">\n";
        int counter = 0;
        for(auto&& value : data)
        {
            // forward to specialized function
            writeToFile_(value);

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file_ << std::endl;
                counter = 0;
            }
        }
        file_ << "\n</DataArray>\n";
    }

     /*!
     * \brief Writes a scalar to the file
     *
     * \param s The scalar
     */
    void writeToFile_(const Scalar& s)
    {
        file_ << s << " ";
    }

     /*!
     * \brief Writes a vector to the file
     *
     * \param g The vector
     */
    void writeToFile_(const GlobalPosition& g)
    {
        assert(g.size() > 1 && g.size() < 4);
        if(g.size() < 3)
            file_ << g << " 0 ";
        else
            file_ << g;
    }


     /*!
     * \brief Returns the file name for each timestep
     */
    std::string fileName_() const
    {
        std::ostringstream oss;
        oss << /*problem_.name() <<*/ "face-" << std::setw(5) << std::setfill('0') << curWriterNum_ << ".vtp";
        return oss.str();
    }

    const std::vector<GlobalPosition>& coordinates_;
    std::ofstream file_;
    int curWriterNum_{0};

    std::list<ScalarLocalFunction> scalarPointData_;
    std::list<VectorLocalFunction> vectorPointData_;
};
} // end namespace Dumux

#endif
