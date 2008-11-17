// $Id$

/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Simplyfies writing multi-file VTK datasets.
 */
#ifndef VTK_MULTI_WRITER_HH
#define VTK_MULTI_WRITER_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <boost/format.hpp>

#include <list>
#include <iostream>
#include <string>

namespace Dune {
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
    public:

        typedef Dune::VTKWriter<GridView> VtkWriter;
        VtkMultiWriter(const std::string &simName = "", std::string multiFileName = "")
            {
                simName_ = (simName.empty())?"sim":simName;

                if (multiFileName.empty())
                    multiFileName = (boost::format("%s.pvd")%simName).str();

                writerNum_ = 0;

                beginMultiFile_(multiFileName);
            }

        ~VtkMultiWriter()
            {
                endMultiFile_();
                multiFile_.close();
            }

        /*!
         * \brief Called when ever a new timestep or a new grid
         *        must be written.
         */
        void beginTimestep(double t, const GridView &gridView)
            {
                curWriter_ = new VtkWriter(gridView);
                ++writerNum_;
                curTime_ = t;
                curGridView_ = &gridView;

                curOutFileName_ = (boost::format("%s-%05d")
                                   %simName_%writerNum_).str();
                const char *suffix = (GridView::dimension == 1)?"vtp":"vtu";
                multiFile_ << boost::format("   <DataSet timestep=\"%lf\" file=\"%s.%s\"/>\n")
                    %curTime_%curOutFileName_%suffix;
            };

        /*!
         * \brief Write a vertex centered vector field to disk.
         */
        template <class Scalar, int nComp>
        Dune::BlockVector<Dune::FieldVector<Scalar, nComp> > *createField(int nEntities)
            {
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, nComp> > VectorField;

                VtkVectorFieldStoreImpl_<VectorField> *vfs =
                    new VtkVectorFieldStoreImpl_<VectorField>(nEntities);
                vectorFields_.push_back(vfs);
                return &(vfs->vf);
            }

        /*!
         * \brief Add a finished vertex centered vector field to the
         *        output. The field must have been created using
         *        createField() and may not be modified after calling
         *        this method.
         */
        template <class VectorField>
        void addVertexData(VectorField *field, const char *name)
            {
                curWriter_->addVertexData(*field, name);
            }

        /*!
         * \brief Add a finished cell centered vector field to the
         *        output. The field must have been created using
         *        createField() and may not be modified after calling
         *        this method.
         */
        template <class VectorField>
        void addCellData(VectorField *field, const char *name)
            {
                curWriter_->addCellData(*field, name);
            }

        /*!
         * \brief Evaluates a single component of a function defined on the grid at the
         *        vertices and appends it to the writer.
         */
        template <class Function>
        void addScalarVertexFunction(const char *name,
                                     const Function &fn,
                                     int comp)
            {
                /*
                // useful typedefs
                typedef typename Function::RangeFieldType                      Scalar;
                typedef typename GridView::Traits::template Codim<GridView::dimension> VertexTraits;
                typedef typename VertexTraits::Entity                          Vertex;
                typedef typename VertexTraits::LeafIterator                    VertexIterator;
                typedef Dune::ReferenceElement<typename GridView::ctype, 0>        VertexReferenceElement;
                typedef Dune::ReferenceElements<typename GridView::ctype, 0>       VertexReferenceElements;
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> >       ScalarField;

                // create a vertex based scalar field.
                ScalarField *field = createField<Scalar, 1>(vertexMap.size());
                std::vector<bool> vertexVisited(vertexMap.size(), false);

                // fill the Scalar field
                VertexIterator it = curGridView_->template leafbegin<GridView::dimension>();
                VertexIterator endIt = curGridView_->template leafend<GridView::dimension>();
                for (; it != endIt; ++it) {
                    // extract the current solution's Sn component
                    const VertexReferenceElement &refElem =
                        VertexReferenceElements::general(it->geometry().type());
                    Scalar compValue = fn.evallocal(comp,
                                                    *it,
                                                    refElem.position(0,0));

                    // find out the cell's index
                    unsigned vertexIndex = vertexMap.map(*it);
                    (*field)[vertexIndex] = compValue;
                }

                addVertexData(field, name);
                */

                // this is pretty hacky as it assumes that the mapping
                // to the vertices is the same for the function a and
                // the vertex mapper
                typedef typename Function::RangeFieldType                Scalar;
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

                unsigned nVerts = (*fn).size();
                ScalarField *field = createField<Scalar, 1>(nVerts);
                for (int i = 0; i < (*fn).size(); i++) {
                    (*field)[i] = (*fn)[i][comp];
                }

                addVertexData(field, name);
            }

        /*!
         * \brief Evaluates a single component of a function defined on the grid at the
         *        cell centers and appends it to the writer.
         */
        template <class Function, class CellMap>
        void addScalarCellFunction(const char *name,
                                   const Function &fn,
                                   const CellMap &cellMap,
                                   int comp)
            {
                // some typedefs

                typedef typename Function::RT                                        Scalar;
                typedef typename GridView::template Codim<0>::Entity                 Cell;
                typedef typename GridView::template Codim<0>::Iterator               CellIterator;
                typedef Dune::ReferenceElement<typename GridView::ctype, GridView::dimgrid>  CellReferenceElement;
                typedef Dune::ReferenceElements<typename GridView::ctype, GridView::dimgrid> CellReferenceElements;
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> >             ScalarField;

                // create a cell based scalar field.
                ScalarField *field = createField<Scalar, 1>(cellMap.size());

                // fill the Scalar field
                CellIterator it = curGridView_->template begin<0>();
                CellIterator endIt = curGridView_->template end<0>();
                for (; it != endIt; ++it) {
                    // extract the current solution's Sn component
                    const CellReferenceElement &refElem =
                        CellReferenceElements::general(it->geometry().type());
                    Scalar compValue = fn.evallocal(comp,
                                                    *it,
                                                    refElem.position(0,0));

                    // find out the cell's index
                    unsigned cellIndex = cellMap.map(*it);
                    (*field)[cellIndex] = compValue;
                }

                addCellData(field, name);
            }

        /*!
         * \brief Finalizes the current writer.
         *
         * This means that everything will be written to disk.
         */
        void endTimestep()
            {
                curWriter_->write(curOutFileName_.c_str(),
                                  Dune::VTKOptions::ascii);
                delete curWriter_;
                while (vectorFields_.begin() != vectorFields_.end()) {
                    delete vectorFields_.front();
                    vectorFields_.pop_front();
                }

                // temporarily write the closing XML mumbo-jumbo to
                // the mashup file so that the data set can be loaded
                // even if the programm is aborted
                endMultiFile_();
            };


    private:
        void beginMultiFile_(const std::string &multiFileName)
            {
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

        void endMultiFile_()
            {
                // make sure that we always have a working meta file
                std::ofstream::pos_type pos = multiFile_.tellp();
                multiFile_ <<
                    " </Collection>\n"
                    "</VTKFile>\n";
                multiFile_.seekp(pos);
                multiFile_.flush();
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
        class VtkVectorFieldStoreBase_
        {
        public:
            virtual ~VtkVectorFieldStoreBase_()
                {}
        };

        template <class VF>
        class VtkVectorFieldStoreImpl_ : public VtkVectorFieldStoreBase_
        {
        public:
            VtkVectorFieldStoreImpl_(int size)
                : vf(size)
                { }
            VF vf;
        };
        // end hack
        ////////////////////////////////////

        std::ofstream   multiFile_;
        std::string     simName_;

        double          curTime_;
        VtkWriter     * curWriter_;
        const GridView* curGridView_;
        std::string     curOutFileName_;
        int             writerNum_;

        std::list<VtkVectorFieldStoreBase_*> vectorFields_;
    };
}

#endif
