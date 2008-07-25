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
#ifndef STUPID_VTK_MULTI_WRITER_HH
#define STUPID_VTK_MULTI_WRITER_HH

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
    template<class Grid>
    class VtkMultiWriter
    {
    public:
        typedef Dune::VTKWriter<Grid> VtkWriter;

        VtkMultiWriter(const std::string &simName = "", std::string multiFileName = "")
            {
                _simName = (simName.empty())?"sim":simName;

                if (multiFileName.empty())
                    multiFileName = (boost::format("%s.pvd")%simName).str();

                _writerNum = 0;

                _beginMultiFile(multiFileName);
            }

        ~VtkMultiWriter()
            {
                _endMultiFile();
                _multiFile.close();
            }

        /*!
         * \brief Called when ever a new timestep or a new grid
         *        must be written.
         */
        void beginTimestep(double t, const Grid &grid)
            {
                _curWriter = new VtkWriter(grid);
                ++_writerNum;
                _curTime = t;
                _curGrid = &grid;

                _curOutFileName = (boost::format("%s-%05d")
                                   %_simName%_writerNum).str();
                const char *suffix = (Grid::dimension == 1)?"vtp":"vtu";
                _multiFile << boost::format("   <DataSet timestep=\"%lf\" file=\"%s.%s\"/>\n")
                    %_curTime%_curOutFileName%suffix;
            };

        /*!
         * \brief Write a vertex centered vector field to disk.
         */
        template <class Scalar, int nComp>
        Dune::BlockVector<Dune::FieldVector<Scalar, nComp> > *createField(int nEntities)
            {
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, nComp> > VectorField;

                _VtkVectorFieldStoreImpl<VectorField> *vfs =
                    new _VtkVectorFieldStoreImpl<VectorField>(nEntities);
                _vectorFields.push_back(vfs);
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
                _curWriter->addVertexData(*field, name);
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
                _curWriter->addCellData(*field, name);
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
                typedef typename Grid::Traits::template Codim<Grid::dimension> VertexTraits;
                typedef typename VertexTraits::Entity                          Vertex;
                typedef typename VertexTraits::LeafIterator                    VertexIterator;
                typedef Dune::ReferenceElement<typename Grid::ctype, 0>        VertexReferenceElement;
                typedef Dune::ReferenceElements<typename Grid::ctype, 0>       VertexReferenceElements;
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> >       ScalarField;

                // create a vertex based scalar field.
                ScalarField *field = createField<Scalar, 1>(vertexMap.size());
                std::vector<bool> vertexVisited(vertexMap.size(), false);

                // fill the Scalar field
                VertexIterator it = _curGrid->template leafbegin<Grid::dimension>();
                VertexIterator endIt = _curGrid->template leafend<Grid::dimension>();
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
                typedef typename Grid::Traits::template Codim<0>                     CellTraits;
                typedef typename CellTraits::Entity                                  Cell;
                typedef typename CellTraits::LeafIterator                            CellIterator;
                typedef Dune::ReferenceElement<typename Grid::ctype, Grid::dimgrid>  CellReferenceElement;
                typedef Dune::ReferenceElements<typename Grid::ctype, Grid::dimgrid> CellReferenceElements;
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> >             ScalarField;

                // create a cell based scalar field.
                ScalarField *field = createField<Scalar, 1>(cellMap.size());

                // fill the Scalar field
                CellIterator it = _curGrid->template leafbegin<0>();
                CellIterator endIt = _curGrid->template leafend<0>();
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
                _curWriter->write(_curOutFileName.c_str(),
                                  Dune::VTKOptions::ascii);
                delete _curWriter;
                while (_vectorFields.begin() != _vectorFields.end()) {
                    delete _vectorFields.front();
                    _vectorFields.pop_front();
                }
                
                // temporarily write the closing XML mumbo-jumbo to
                // the mashup file so that the data set can be loaded
                // even if the programm is aborted
                _endMultiFile();
            };


    private:
        void _beginMultiFile(const std::string &multiFileName)
            {
                // generate one meta vtk-file holding the individual timesteps
                _multiFile.open(multiFileName.c_str());
                _multiFile <<
                    "<?xml version=\"1.0\"?>\n"
                    "<VTKFile type=\"Collection\"\n"
                    "         version=\"0.1\"\n"
                    "         byte_order=\"LittleEndian\"\n"
                    "         compressor=\"vtkZLibDataCompressor\">\n"
                    " <Collection>\n";
            }

        void _endMultiFile()
            {
                // make sure that we always have a working meta file
                std::ofstream::pos_type pos = _multiFile.tellp();
                _multiFile <<
                    " </Collection>\n"
                    "</VTKFile>\n";
                _multiFile.seekp(pos);
                _multiFile.flush();
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
        class _VtkVectorFieldStoreBase
        {
        public:
            virtual ~_VtkVectorFieldStoreBase()
                {}
        };

        template <class VF>
        class _VtkVectorFieldStoreImpl : public _VtkVectorFieldStoreBase
        {
        public:
            _VtkVectorFieldStoreImpl(int size)
                : vf(size)
                { }
            VF vf;
        };
        // end hack
        ////////////////////////////////////

        std::ofstream  _multiFile;
        std::string    _simName;

        double         _curTime;
        VtkWriter     *_curWriter;
        const Grid    *_curGrid;
        std::string    _curOutFileName;
        int            _writerNum;

        std::list<_VtkVectorFieldStoreBase*> _vectorFields;
    };
}

#endif
