/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Simplyfies writing multi-file VTK datasets.
 */
#ifndef VTK_MULTI_WRITER_HH
#define VTK_MULTI_WRITER_HH

#include "vtknestedfunction.hh"

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dumux/common/valgrind.hh>

#include <boost/format.hpp>

#if HAVE_MPI
#include <mpi.h>
#endif

#include <list>
#include <iostream>
#include <string>

#include <limits>

namespace Dumux {
/*!
 * \brief Simplyfies writing multi-file VTK datasets.
 *
 * This class automatically keeps the meta file up to date and
 * simplifies writing datasets consisting of multiple files. (i.e.
 * multiple timesteps or grid refinements within a timestep.)
 */
template<class GridView>
class VtkMultiWriter
{
    typedef typename GridView::Grid Grid;
    enum { dim = GridView::dimension };

#if DUNE_VERSION_NEWER_REV(DUNE_GRID, 2, 1, 0)
    // DUNE 2.1 and above
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGVertexLayout> VertexMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;
#else
    // DUNE 2.0 and below
    template<int dim>
    struct VertexLayout {
        bool contains (Dune::GeometryType gt) const
        { return gt.dim() == 0; } };
    template<int dim>
    struct ElementLayout {
        bool contains (Dune::GeometryType gt) const
        { return gt.dim() == dim; } };

    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, VertexLayout> VertexMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, ElementLayout> ElementMapper;
#endif

    // this constructor won't work anymore. Please use the variant
    // below which also includes the GridView as an argument.
    DUNE_DEPRECATED
    VtkMultiWriter(const std::string &simName = "",
                   std::string multiFileName = "")
    {}

public:
    typedef Dune::VTKWriter<GridView> VtkWriter;
    VtkMultiWriter(const GridView &gridView,
                   const std::string &simName = "",
                   std::string multiFileName = "")
        : gridView_(gridView)
        , elementMapper_(gridView)
        , vertexMapper_(gridView)
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
     * \brief Called when ever a new timestep or a new grid
     *        must be written.
     */
    void beginWrite(double t)
    {
        if (!multiFile_.is_open()) {
            startMultiFile_(multiFileName_);
        }


        curWriter_ = new VtkWriter(gridView_,
                                   Dune::VTKOptions::conforming);
        ++curWriterNum_;

        curTime_ = t;
        curOutFileName_ = fileName_();
    }

    void beginTimestep(double t, const GridView &gridView)
        DUNE_DEPRECATED // use beginWrite()
    { gridChanged(); beginWrite(t); }

    /*!
     * \brief Allocate a managed buffer for a scalar field
     *
     * The buffer will be deleted automatically after the data has
     * been written by to disk.
     */
    template <class Scalar = double, int nComp = 1>
    Dune::BlockVector<Dune::FieldVector<Scalar, nComp> > *allocateManagedBuffer(int nEntities)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, nComp> > VectorField;

        ManagedVectorField_<VectorField> *vfs =
            new ManagedVectorField_<VectorField>(nEntities);
        managedObjects_.push_back(vfs);
        return &(vfs->vf);
    }

    template <class Scalar, int nComp>
    DUNE_DEPRECATED // use allocateManagedBuffer() instead!
    Dune::BlockVector<Dune::FieldVector<Scalar, nComp> > *createField(int nEntities)
    {  return allocateManagedBuffer<Scalar, nComp>(nEntities);   }

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
    void attachVertexData(DataBuffer &buf, const char *name, int nComps=1)
    {
        typedef typename VtkWriter::VTKFunctionPtr FunctionPtr;
        typedef Dumux::VtkNestedFunction<Grid, VertexMapper, DataBuffer> VtkFn;
        FunctionPtr fnPtr(new VtkFn(name,
                                    gridView_.grid(),
                                    vertexMapper_,
                                    buf,
                                    /*codim=*/dim,
                                    nComps));
        curWriter_->addVertexData(fnPtr);
    }

    template <class DataArray>
    DUNE_DEPRECATED // use attachVertexData() instead!
    void addVertexData(DataArray *field, const char *name, int nComps=1)
    {
      curWriter_->addVertexData(*field, name);
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
    void attachCellData(DataBuffer &buf, const char *name, int nComps = 1)
    {
        typedef typename VtkWriter::VTKFunctionPtr FunctionPtr;
        typedef Dumux::VtkNestedFunction<Grid, ElementMapper, DataBuffer> VtkFn;
        FunctionPtr fnPtr(new VtkFn(name,
                                    gridView_.grid(),
                                    elementMapper_,
                                    buf,
                                    /*codim=*/0,
                                    nComps));
        curWriter_->addCellData(fnPtr);
    }

    template <class VectorField>
    DUNE_DEPRECATED // use attachCellData() instead
    void addCellData(VectorField *field, const char *name, int nComps = 1)
    {
        curWriter_->addCellData(*field, name);
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
            curWriter_->write(curOutFileName_.c_str(),
                              Dune::VTKOptions::ascii);

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
        else
            -- curWriterNum_;

        // discard managed objects and the current VTK writer
        delete curWriter_;
        while (managedObjects_.begin() != managedObjects_.end()) {
            delete managedObjects_.front();
            managedObjects_.pop_front();
        }

        // temporarily write the closing XML mumbo-jumbo to the mashup
        // file so that the data set can be loaded even if the
        // simulation is aborted (or not yet finished)
        finishMultiFile_();
    };

    void endTimestep(bool onlyDiscard=false)
        DUNE_DEPRECATED // use endWrite()
    { endWrite(onlyDiscard); }


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
        };
        res.deserializeSectionEnd();
    }


private:
    std::string fileName_()
    {
        return (boost::format("%s-%05d")
                %simName_%curWriterNum_).str();
    }

    std::string fileName_(int rank)
    {
        if (commSize_ > 1) {
            return (boost::format("s%04d:p%04d:%s-%05d")
                    %commSize_
                    %rank
                    %simName_
                    %curWriterNum_).str();
        }
        else {
            return (boost::format("%s-%05d")
                    %simName_%curWriterNum_).str();
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
            // generate one meta vtk-file holding the individual timesteps
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

    // make sure the field is well defined if running under valgrind
    // and make sure that all values can be displayed by paraview
    template <class DataBuffer>
    void sanitizeBuffer_(DataBuffer &b, int nComps)
    {
        for (int i = 0; i < b.size(); ++i) {
            for (int j = 0; j < nComps; ++j) {
                Valgrind::CheckDefined(b[i][j]);

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
    // HACK: when ever we attach some data we need to copy the
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

    /** \todo Please doc me! */

    class ManagedObject_
    {
    public:
        virtual ~ManagedObject_()
        {}
    };

    /** \todo Please doc me! */

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

    VtkWriter *curWriter_;
    double curTime_;
    std::string curOutFileName_;
    int curWriterNum_;

    std::list<ManagedObject_*> managedObjects_;
};
}

#endif
