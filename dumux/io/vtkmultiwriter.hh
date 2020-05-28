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
 * \brief Simplifies writing multi-file VTK datasets.
 */
#ifndef VTK_MULTI_WRITER_HH
#define VTK_MULTI_WRITER_HH

#warning "This header is deprecated. Use the new vtkoutputmodule."

#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <string>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/function.hh>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Provides a vector-valued function using Dune::FieldVectors
 *        as elements.
 * DEPRECATED will be removed once this header is removed
 */
template <class GridView, class Mapper, class Buffer>
class VtkNestedFunction : public Dune::VTKFunction<GridView>
{
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    VtkNestedFunction(std::string name,
                      const GridView &gridView,
                      const Mapper &mapper,
                      const Buffer &buf,
                      int codim,
                      int numComp)
        : name_(name)
        , gridView_(gridView)
        , mapper_(mapper)
        , buf_(buf)
        , codim_(codim)
        , numComp_(numComp)
    {
        assert(std::size_t(buf_.size()) == std::size_t(mapper_.size()));
    }

    virtual std::string name () const
    { return name_; }

    virtual int ncomps() const
    { return numComp_; }

    virtual double evaluate(int mycomp,
                            const Element &element,
                            const Dune::FieldVector< ctype, dim > &xi) const
    {
        int idx;
        if (codim_ == 0) {
            // cells. map element to the index
            idx = mapper_.index(element);
        }
        else if (codim_ == dim) {
            // find vertex which is closest to xi in local
            // coordinates. This code is based on Dune::P1VTKFunction
            double min=1e100;
            int imin=-1;
            int n = element.subEntities(dim);

            for (int i=0; i < n; ++i)
            {
                Dune::FieldVector<ctype,dim> local = referenceElement(element).position(i,dim);
                local -= xi;
                if (local.infinity_norm()<min)
                {
                    min = local.infinity_norm();
                    imin = i;
                }
            }

            // map vertex to an index
            idx = mapper_.subIndex(element, imin, codim_);
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Only element and vertex based vector "
                       " fields are supported so far.");

        double val = buf_[idx][mycomp];
        using std::abs;
        if (abs(val) < std::numeric_limits<float>::min())
            val = 0;

        return val;
    }

private:
    const std::string name_;
    const GridView gridView_;
    const Mapper &mapper_;
    const Buffer &buf_;
    int codim_;
    int numComp_;
};

/*!
 * \ingroup InputOutput
 * \brief Simplifies writing multi-file VTK datasets.
 *
 * This class automatically keeps the meta file up to date and
 * simplifies writing datasets consisting of multiple files. (i.e.
 * multiple time steps or grid refinements within a time step.)
 * \todo This class can most likely be replaced by Dune::VTKSequenceWriter
 */
template<class GridView, Dune::VTK::OutputType OutputValue = Dune::VTK::ascii >
class [[deprecated("Use VtkOutputModule instead!")]] VtkMultiWriter
{
    enum { dim = GridView::dimension };
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
public:
    using VtkWriter = Dune::VTKWriter<GridView>;
    VtkMultiWriter(const GridView &gridView,
                   const std::string &simName = "",
                   std::string multiFileName = "")
        : gridView_(gridView)
        , elementMapper_(gridView, Dune::mcmgElementLayout())
        , vertexMapper_(gridView, Dune::mcmgVertexLayout())
    {
        simName_ = (simName.empty())?"sim":simName;
        multiFileName_ = multiFileName;
        if (multiFileName_.empty()) {
            multiFileName_ = simName_;
            multiFileName_ += ".pvd";
        }

        curWriterNum_ = 0;

        commRank_ = 0;
        commSize_ = 1;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &commRank_);
        MPI_Comm_size(MPI_COMM_WORLD, &commSize_);
#endif
    }

    ~VtkMultiWriter()
    {
        finishMultiFile_();

        if (commRank_ == 0)
            multiFile_.close();
    }

    /*!
     * \brief Returns the number of the current VTK file.
     */
    int curWriterNum() const
    { return curWriterNum_; }

    /*!
     * \brief Updates the internal data structures after mesh
     *        refinement.
     *
     * If the grid changes between two calls of beginWrite(), this
     * method _must_ be called before the second beginWrite()!
     */
    void gridChanged()
    {
        elementMapper_.update();
        vertexMapper_.update();
    }

    /*!
     * \brief Called when ever a new time step or a new grid
     *        must be written.
     */
    void beginWrite(double t)
    {
        if (!multiFile_.is_open()) {
            startMultiFile_(multiFileName_);
        }

        curWriter_ = std::make_shared<VtkWriter>(gridView_, Dune::VTK::conforming);
        curTime_ = t;
        curOutFileName_ = fileName_();
    }

    /*!
     * \brief Allocate a managed buffer for a scalar field
     *
     * The buffer will be deleted automatically after the data has
     * been written by to disk.
     */
    template <class Scalar, int nComp>
    Dune::BlockVector<Dune::FieldVector<Scalar, nComp> > *allocateManagedBuffer(int nEntities)
    {
        using VectorField = Dune::BlockVector<Dune::FieldVector<Scalar, nComp> >;

        std::shared_ptr<ManagedVectorField_<VectorField> > vfs =
            std::make_shared<ManagedVectorField_<VectorField> >(nEntities);
        managedObjects_.push_back(vfs);
        return &(vfs->vf);
    }

    template <class Scalar>
    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > *allocateManagedBuffer(int nEntities)
    { return allocateManagedBuffer<Scalar, 1>(nEntities); }
    Dune::BlockVector<Dune::FieldVector<double, 1> > *allocateManagedBuffer(int nEntities)
    { return allocateManagedBuffer<double, 1>(nEntities); }

    /*!
     * \brief Add a finished vertex centered vector field to the
     *        output.
     * \brief Add a vertex-centered quantity to the output.
     *
     * If the buffer is managed by the VtkMultiWriter, it must have
     * been created using createField() and may not be used by
     * anywhere after calling this method. After the data is written
     * to disk, it will be deleted automatically.
     *
     * If the buffer is not managed by the MultiWriter, the buffer
     * must exist at least until the call to endWrite()
     * finishes.
     *
     * In both cases, modifying the buffer between the call to this
     * method and endWrite() results in _undefined behaviour_.
     */
    template <class DataBuffer>
    void attachVertexData(DataBuffer &buf, std::string name, int nComps=1)
    {
        sanitizeBuffer_(buf, nComps);

        using FunctionPtr = std::shared_ptr<const typename VtkWriter::VTKFunction>;
        using VtkFn = VtkNestedFunction<GridView, VertexMapper, DataBuffer>;
        FunctionPtr fnPtr(new VtkFn(name,
                                    gridView_,
                                    vertexMapper_,
                                    buf,
                                    /*codim=*/dim,
                                    nComps));
        curWriter_->addVertexData(fnPtr);
    }

    /*!
     * \brief Add a cell centered quantity to the output.
     *
     * If the buffer is managed by the VtkMultiWriter, it must have
     * been created using createField() and may not be used by
     * anywhere after calling this method. After the data is written
     * to disk, it will be deleted automatically.
     *
     * If the buffer is not managed by the MultiWriter, the buffer
     * must exist at least until the call to endWrite()
     * finishes.
     *
     * In both cases, modifying the buffer between the call to this
     * method and endWrite() results in _undefined behaviour_.
     */
    template <class DataBuffer>
    void attachCellData(DataBuffer &buf, std::string name, int nComps = 1)
    {
        sanitizeBuffer_(buf, nComps);

        using FunctionPtr = std::shared_ptr<const typename VtkWriter::VTKFunction>;
        using VtkFn = VtkNestedFunction<GridView, ElementMapper, DataBuffer>;
        FunctionPtr fnPtr(new VtkFn(name,
                                    gridView_,
                                    elementMapper_,
                                    buf,
                                    /*codim=*/0,
                                    nComps));
        curWriter_->addCellData(fnPtr);
    }

    /*!
     * \brief Add data associated with degrees of freedom to the output.
     *
     * This is a wrapper for the functions attachVertexData and attachCelldata.
     * Depending on the value of \a isVertexData, the appropriate function
     * is selected.
     */
    template <class DataBuffer>
    void attachDofData(DataBuffer &buf, std::string name, bool isVertexData, int nComps = 1)
    {
        if (isVertexData)
            attachVertexData(buf, name, nComps);
        else
            attachCellData(buf, name, nComps);
    }

    /*!
     * \brief Finalizes the current writer.
     *
     * This means that everything will be written to disk, except if
     * the onlyDiscard argument is true. In this case only all managed
     * buffers are deleted, but no output is written.
     */
    void endWrite(bool onlyDiscard = false)
    {
        if (!onlyDiscard) {
            ++curWriterNum_;
            curWriter_->write(curOutFileName_.c_str(), OutputValue);

            // determine name to write into the multi-file for the
            // current time step
            std::string fileName;
            std::string suffix = fileSuffix_();
            if (commSize_ == 1) {
                fileName = curOutFileName_;
                multiFile_.precision(16);
                multiFile_ << "   <DataSet timestep=\""
                           << curTime_
                           << "\" file=\""
                           << fileName << "." << suffix << "\"/>\n";
            }
            if (commSize_ > 1 && commRank_ == 0) {
                // only the first process updates the multi-file
                for (int part=0; part < commSize_; ++part) {
                    fileName = fileName_(part);
                    multiFile_.precision(16);
                    multiFile_ << "   <DataSet "
                               << " part=\"" << part << "\""
                               << " timestep=\"" << curTime_ << "\""
                               << " file=\"" << fileName << "." << suffix << "\"/>\n";
                }
            }
        }

        // discard managed objects
        while (managedObjects_.begin() != managedObjects_.end()) {
            managedObjects_.pop_front();
        }

        // temporarily write the closing XML mumbo-jumbo to the mashup
        // file so that the data set can be loaded even if the
        // simulation is aborted (or not yet finished)
        finishMultiFile_();
    }

    /*!
     * \brief Write the multi-writer's state to a restart file.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        res.serializeSectionBegin("VTKMultiWriter");
        res.serializeStream() << curWriterNum_ << "\n";

        if (commRank_ == 0) {
            size_t fileLen = 0;
            size_t filePos = 0;
            if (multiFile_.is_open()) {
                // write the meta file into the restart file
                filePos = multiFile_.tellp();
                multiFile_.seekp(0, std::ios::end);
                fileLen = multiFile_.tellp();
                multiFile_.seekp(filePos);
            }

            res.serializeStream() << fileLen << "  " << filePos << "\n";

            if (fileLen > 0) {
                std::ifstream multiFileIn(multiFileName_.c_str());
                char *tmp = new char[fileLen];
                multiFileIn.read(tmp, fileLen);
                res.serializeStream().write(tmp, fileLen);
                delete[] tmp;
            }
        }

        res.serializeSectionEnd();
    }

    /*!
     * \brief Read the multi-writer's state from a restart file.
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.deserializeSectionBegin("VTKMultiWriter");
        res.deserializeStream() >> curWriterNum_;

        if (commRank_ == 0) {
            std::string dummy;
            std::getline(res.deserializeStream(), dummy);

            // recreate the meta file from the restart file
            size_t filePos, fileLen;
            res.deserializeStream() >> fileLen >> filePos;
            std::getline(res.deserializeStream(), dummy);
            if (multiFile_.is_open())
                multiFile_.close();

            if (fileLen > 0) {
                multiFile_.open(multiFileName_.c_str());

                char *tmp = new char[fileLen];
                res.deserializeStream().read(tmp, fileLen);
                multiFile_.write(tmp, fileLen);
                delete[] tmp;
            }

            multiFile_.seekp(filePos);
        }
        else {
            std::string tmp;
            std::getline(res.deserializeStream(), tmp);
        }
        res.deserializeSectionEnd();
    }

private:
    std::string fileName_()
    {
        std::ostringstream oss;
        oss << simName_ << "-" << std::setw(5) << std::setfill('0') << curWriterNum_;
        return oss.str();
    }

    std::string fileName_(int rank)
    {
        if (commSize_ > 1) {
            std::ostringstream oss;
            oss << "s" << std::setw(4) << std::setfill('0') << commSize_
                << "-p" << std::setw(4) << std::setfill('0') << rank
                << "-" << simName_ << "-"
                << std::setw(5) << curWriterNum_;
            return oss.str();

            return oss.str();
        }
        else {
            std::ostringstream oss;
            oss << simName_ << "-" << std::setw(5) << std::setfill('0') << curWriterNum_;
            return oss.str();
        }
    }

    std::string fileSuffix_()
    {
        return (GridView::dimension == 1)?"vtp":"vtu";
    }


    void startMultiFile_(const std::string &multiFileName)
    {
        // only the first process writes to the multi-file
        if (commRank_ == 0) {
            // generate one meta vtk-file holding the individual time steps
            multiFile_.open(multiFileName.c_str());
            multiFile_ <<
                "<?xml version=\"1.0\"?>\n"
                "<VTKFile type=\"Collection\"\n"
                "         version=\"0.1\"\n"
                "         byte_order=\"LittleEndian\"\n"
                "         compressor=\"vtkZLibDataCompressor\">\n"
                " <Collection>\n";
        }
    }

    void finishMultiFile_()
    {
        // only the first process writes to the multi-file
        if (commRank_ == 0) {
            // make sure that we always have a working meta file
            std::ofstream::pos_type pos = multiFile_.tellp();
            multiFile_ <<
                " </Collection>\n"
                "</VTKFile>\n";
            multiFile_.seekp(pos);
            multiFile_.flush();
        }
    }

    // make sure that all values can be displayed by paraview
    template <class DataBuffer>
    void sanitizeBuffer_(DataBuffer &b, int nComps)
    {
        for (unsigned int i = 0; i < b.size(); ++i) {
            for (int j = 0; j < nComps; ++j) {
                // set values which are too small to 0 to avoid
                // problems with paraview
                if (std::abs(b[i][j])
                    < std::numeric_limits<float>::min())
                {
                    b[i][j] = 0.0;
                }
            }
        }
    }

    //////////////////////////////
    // Trick: when ever we attach some data we need to copy the
    //       vector field (that's because Dune::VTKWriter is not
    //       able to write fields one at a time and using
    //       VTKWriter::add*Data doesn't copy the data's
    //       representation so that once we arrive at the point
    //       where we want to write the data to disk, it might
    //       exist anymore). The problem we encounter there is
    //       that add*Data accepts arbitrary types as vector
    //       fields, but we need a single type for the linked list
    //       which keeps track of the data added. The trick we use
    //       here is to define a non-template base class with a
    //       virtual destructor for the type given to the linked
    //       list and a derived template class which actually
    //       knows the type of the vector field it must delete.

    class ManagedObject_
    {
    public:
        virtual ~ManagedObject_()
        {}
    };

    template <class VF>
    class ManagedVectorField_ : public ManagedObject_
    {
    public:
        ManagedVectorField_(int size)
            : vf(size)
        { }
        VF vf;
    };
    // end hack
    ////////////////////////////////////

    const GridView gridView_;
    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    std::string simName_;
    std::ofstream multiFile_;
    std::string multiFileName_;

    int commSize_; // number of processes in the communicator
    int commRank_; // rank of the current process in the communicator

    std::shared_ptr<VtkWriter> curWriter_;
    double curTime_;
    std::string curOutFileName_;
    int curWriterNum_;

    std::list<std::shared_ptr<ManagedObject_> > managedObjects_;
};
}

#endif
