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
 * \brief A VTK writer specialized for staggered grid implementations with dofs on the faces
 */
#ifndef DUMUX_POINTCLOUD_VTK_WRITER_HH
#define DUMUX_POINTCLOUD_VTK_WRITER_HH

#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <iomanip>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/path.hh>
#include <dune/grid/io/file/vtk/common.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 *
 * Handles the output of scalar and vector fields to VTK formatted file for multiple
 * variables and timesteps. Certain predefined fields can be registered on problem / model
 * initialization and/or be turned on/off using the designated properties. Additionally
 * non-standardized scalar and vector fields can be added to the writer manually.
 */
template<class Scalar, class GlobalPosition>
class PointCloudVtkWriter
{
    // GlobalPosition is used for the point coordinates, DimWorldVector for the actual data.
    // GlobalPosition's ctype does not necessarily equal Scalar.
    using DimWorldVector = Dune::FieldVector<Scalar, GlobalPosition::size()>;

    static constexpr unsigned int precision = 6;
    static constexpr unsigned int numBeforeLineBreak = 15;

     /*!
     * \brief A class holding a data container and additional information
     */
    template<class ContainerType>
    class VTKFunction
    {
    public:
        /*!
        * \brief A class holding a data container and additional information
        *
        * \param data The data container
        * \param name The name of the data set
        * \param numComponents The number of components of the data set
        */
        VTKFunction(const ContainerType& data, const std::string& name, const int numComponents) : data_(data), name_(name), numComponents_(numComponents)
        {}

        /*!
        * \brief Returns the name of the data set
        */
        const std::string& name() const
        {
            return name_;
        }

        /*!
        * \brief Returns the number of components of the data set
        */
        int numComponents() const
        {
            return numComponents_;
        }

        /*!
        * \brief Allows random acces to data
        *
        * \param idx The index
        */
        auto& operator() (const int idx) const  { return data_[idx]; }

        decltype(auto) begin() const
        {
            return data_.begin();
        }

        decltype(auto) end() const
        {
            return data_.end();
        }

        /*!
        * \brief Returns the size of the data container
        */
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
    using ScalarFunction = VTKFunction<std::vector<Scalar>>;
    using VectorFunction = VTKFunction<std::vector<DimWorldVector>>;


    PointCloudVtkWriter(const std::vector<GlobalPosition>& coordinates) : coordinates_(coordinates)
    {}

     /*!
     * \brief Create an output file
     *
     * \param name The base name
     * \param type The output type
     */
    void write(const std::string& name, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        auto filename = getSerialPieceName(name, "");

        std::ofstream file;
        file.open(filename);
        writeHeader_(file);
        writeCoordinates_(file, coordinates_);
        writeDataInfo_(file);

        for (auto&& data : scalarPointData_)
            writeData_(file, data);

        for (auto&& data :vectorPointData_)
            writeData_(file, data);

        if (!scalarPointData_.empty() || !vectorPointData_.empty())
            file << "</PointData>\n";

        file << "</Piece>\n";
        file << "</PolyData>\n";
        file << "</VTKFile>";

        clear();
        file.close();
    }

    /*!
     * \brief Create an output file in parallel
     *
     * \param name 			Base name of the output files.  This should not
     *                   	contain any directory part and not filename
     *                   	extensions.  It will be used both for each processes
     *                   	piece as well as the parallel collection file
     * \param path  		Directory where to put the parallel collection
     *                      (.pvtu/.pvtp) file.  If it is relative, it is taken
     *                      relative to the current directory
     * \param extendpath 	Directory where to put the piece file (.vtu/.vtp) of
     *                   	this process.  If it is relative, it is taken
     *                   	relative to the directory denoted by path
     * \param type			How to encode the data in the file
     */
    void pwrite(const std::string & name,  const std::string & path, const std::string & extendpath,
                Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        DUNE_THROW(Dune::NotImplemented, "parallel point cloud vtk output not supported yet");
    }

    /*!
     * \brief Add a vector of scalar data that live on arbitrary points to the visualization.
     *
     * \param v The vector containing the data
     * \param name The name of the data set
     */
    void addPointData(const std::vector<Scalar>& v, const std::string &name)
    {
        assert(v.size() == coordinates_.size());
        scalarPointData_.push_back(ScalarFunction(v, name, 1));
    }

    /*!
     * \brief Add a vector of vector data that live on arbitrary points to the visualization.
     *
     * \param v The vector containing the data
     * \param name The name of the data set
     */
    void addPointData(const std::vector<DimWorldVector>& v, const std::string &name)
    {
        assert(v.size() == coordinates_.size());
        vectorPointData_.push_back(VectorFunction(v, name, 3));
    }

    /*!
     * \brief Clears the data lists
     */
    void clear()
    {
        scalarPointData_.clear();
        vectorPointData_.clear();
    }

    /*!
     * \brief Return name of a serial header file
     *
     * \param name     Base name of the VTK output.  This should be without
     *                 any directory parts and without a filename extension.
     * \param path     Directory part of the resulting header name.  May be
     *                 empty, in which case the resulting name will not have a
     *                 directory part.  If non-empty, may or may not have a
     *                 trailing '/'.  If a trailing slash is missing, one is
     *                 appended implicitly.
     */
    std::string getSerialPieceName(const std::string& name,
                                   const std::string& path) const
    {
      static const std::string extension = ".vtp";

      return Dune::concatPaths(path, name+extension);
    }

    /*!
     * \brief Return name of a parallel header file
     *
     * \param name     Base name of the VTK output.  This should be without
     *                 any directory parts and without a filename extension.
     * \param path     Directory part of the resulting header name.  May be
     *                 empty, in which case the resulting name will not have a
     *                 directory part.  If non-empty, may or may not have a
     *                 trailing '/'.  If a trailing slash is missing, one is
     *                 appended implicitly.
     * \param commSize Number of processes writing a parallel vtk output.
     */
    std::string getParallelHeaderName(const std::string& name,
                                      const std::string& path,
                                      int commSize) const
    {
      std::ostringstream s;
      if(path.size() > 0) {
        s << path;
        if(path[path.size()-1] != '/')
          s << '/';
      }
      s << 's' << std::setw(4) << std::setfill('0') << commSize << '-';
      s << name;
      s << ".pvtp";
      return s.str();
    }


private:
    /*!
     * \brief Writes the header to the file
     */
    void writeHeader_(std::ostream& file)
    {
        std::string header = "<?xml version=\"1.0\"?>\n";
                    header += "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
                    header += "<PolyData>\n";
                    header += "<Piece NumberOfLines=\"0\" NumberOfPoints=\"" + std::to_string(coordinates_.size()) + "\">\n";
        file << header;
    }

    /*!
     * \brief Writes information about the data to the file
     */
    void writeDataInfo_(std::ostream& file)
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
                file << "<PointData Scalars=\"" << scalarName << "\" Vectors=\"" << vectorName <<"\">\n";
            else
                file << "<PointData Scalars=\"" << scalarName << "\">\n";
        else if(foundVector)
            file << "<PointData Vectors=\"" << vectorName << "\">\n";
        else
            return;
    }

    /*!
     * \brief Writes the coordinates to the file
     *
     * \param file The output file
     * \param positions Container to store the positions
     */
    void writeCoordinates_(std::ostream& file, const std::vector<GlobalPosition>& positions)
    {
        // write the positions to the file
        file << "<Points>\n";
        file << "<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        int counter = 0;
        for(auto&& x : positions)
        {
            file << x ;

            if(x.size() == 1)
                file << " 0 0 ";
            if(x.size() == 2)
                file << " 0 ";
            if(x.size() == 3)
                file << " ";

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << std::endl;
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
        file << "</Points>\n";
    }

    /*!
     * \brief Writes data to the file
     *
     * \param file The output file
     * \param data The data container which hold the data itself, as well as the name of the data set and the number of its components
     */
    template<class T>
    void writeData_(std::ostream& file, const T& data)
    {
        file << "<DataArray type=\"Float32\" Name=\"" << data.name() << "\" NumberOfComponents=\"" << data.numComponents() << "\" format=\"ascii\">\n";
        int counter = 0;
        for(auto&& value : data)
        {
            // forward to specialized function
            writeToFile_(file, value);

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << std::endl;
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
    }

    /*!
     * \brief Writes a scalar to the file
     *
     * \param file The output file
     * \param s The scalar
     */
    void writeToFile_(std::ostream& file, const Scalar& s)
    {
        file << s << " ";
    }

    /*!
     * \brief Writes a vector to the file
     *
     * \param file The output file
     * \param g The vector
     */
    void writeToFile_(std::ostream& file, const DimWorldVector& g)
    {
        assert(g.size() > 1 && g.size() < 4);
        if(g.size() < 3)
            file << g << " 0 ";
        else
            file << g << " ";
    }

    const std::vector<GlobalPosition>& coordinates_;
    std::list<ScalarFunction> scalarPointData_;
    std::list<VectorFunction> vectorPointData_;
};
} // end namespace Dumux

#endif
