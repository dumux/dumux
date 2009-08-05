// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief Generic class to read/write restart files.
 */
#ifndef DUMUX_RESTART_HH
#define DUMUX_RESTART_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <boost/format.hpp>

#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>


namespace Dune {

/*!
 * \brief Load or save a state of a model to/from the harddisk.
 */
template <class GridView>
class Restart {
    //! \brief Create a magic cookie for restart files, so that it is
    //!        unlikely to load a restart file for an incorrectly.
    static const std::string magicRestartCookie_(const GridView &gridView)
    {
        const std::string gridName = gridView.grid().name();
        const int dim = GridView::dimension;

        int numVertices = gridView.template size(dim);
        int numElements = gridView.template size(0);
        int numEdges = gridView.template size(dim-1);
        int numCPUs = gridView.comm().size();
        int rank = gridView.comm().rank();

        return
            (boost::format("DuMux restart file: "
                           "gridName='%s' "
                           "numCPUs=%d "
                           "myRank=%d "
                           "numElements=%d "
                           "numEdges=%d "
                           "numVertices=%d")
             %gridName
             %numCPUs
             %rank
             %numElements
             %numEdges
             %numVertices).str();
    }

    //! \brief Return the restart file name.
    static const std::string restartFileName_(const GridView &gridView,
                                              const std::string &simName,
                                              double t)
    {
        int rank = gridView.comm().rank();
        return (boost::format("%s_time=%3.1d_rank=%05d.drs")
                %simName
                %t
                %rank).str();
    };


public:
    /*!
     * \brief Write the current state of the model to disk.
     */
    void serializeBegin(const GridView &gridView,
                        const std::string &simName,
                        double t)
    {
        const std::string magicCookie = magicRestartCookie_(gridView);
        const std::string fileName = restartFileName_(gridView, simName, t);

        // open output file and write magic cookie
        outStream_.open(fileName.c_str());
        outStream_.precision(20);

        serializeSection(magicCookie);
    }

    /*!
     * \brief The output stream to write the serialized data.
     */
    std::ostream &serializeStream()
    {
        return outStream_;
    }

    /*!
     * \brief Start a new section in the serialized output.
     */
    void serializeSection(const std::string &cookie)
    {
        outStream_ << cookie << "\n";
    }

    /*!
     * \brief Serialize all leaf entities of a codim in a gridView.
     *
     * The actual work is done by Serializer::serialize(Entity)
     */
    template <int codim, class Serializer>
    void serializeEntities(Serializer &serializer,
                           const GridView &gridView)
    {
        std::string cookie = (boost::format("Entities: Codim %d")%codim).str();
        serializeSection(cookie);

        // write element data
        typedef typename GridView::template Codim<codim>::Iterator Iterator;

        Iterator it = gridView.template begin<codim>();
        const Iterator &endIt = gridView.template end<codim>();
        for (; it != endIt; ++it) {
            serializer.serializeEntity(outStream_, *it);
            outStream_ << "\n";
        };
    }

    /*!
     * \brief Finish the restart file.
     */
    void serializeEnd()
    {
        outStream_.close();
    };


    /*!
     * \brief Start reading a restart file at a certain simulated
     *        time.
     */
    void deserializeBegin(const GridView &grid,
                          const std::string &simName,
                          double t)
    {
        const std::string fileName = restartFileName_(grid, simName, t);

        // open input file and read magic cookie
        inStream_.open(fileName.c_str());
        if (!inStream_.good()) {
            DUNE_THROW(Dune::IOError,
                       "Restart file '"
                       << fileName
                       << "' could not be opened properly");
        }

        // make sure that we don't open an empty file
        inStream_.seekg(0, std::ios::end);
        int pos = inStream_.tellg();
        if (pos == 0) {
            DUNE_THROW(IOError,
                       "Restart file '"
                       << fileName
                       << "' is empty");
        }
        inStream_.seekg(0, std::ios::beg);

        const std::string magicCookie = magicRestartCookie_(grid);
        deserializeSection(magicCookie);
    }

    /*!
     * \brief The input stream to read the data which ought to be
     *        deserialized.
     */
    std::istream &deserializeStream()
    {
        return inStream_;
    }

    /*!
     * \brief Start reading a new section of the restart file.
     */
    void deserializeSection(const std::string &cookie)
    {
        if (!inStream_.good())
            DUNE_THROW(IOError,
                       "Encountered unexpected EOF in restart file.");
        std::string buf;
        std::getline(inStream_, buf);
        if (buf != cookie)
            DUNE_THROW(IOError,
                       "Could not start section '" << cookie << "'");
    };


    /*!
     * \brief Deserialize all leaf entities of a codim in a grid.
     *
     * The actual work is done by Deserializer::deserialize(Entity)
     */
    template <int codim, class Deserializer>
    void deserializeEntities(Deserializer &deserializer,
                             const GridView &gridView)
    {
        std::string cookie = (boost::format("Entities: Codim %d")%codim).str();
        deserializeSection(cookie);

        std::string curLine;

        // read entity data
        typedef typename GridView::template Codim<codim>::Iterator Iterator;
        Iterator it = gridView.template begin<codim>();
        const Iterator &endIt = gridView.template end<codim>();
        for (; it != endIt; ++it) {
            if (!inStream_.good()) {
                DUNE_THROW(IOError,
                           "Restart file is corrupted");
            }

            std::getline(inStream_, curLine);
            std::istringstream curLineStream(curLine);
            deserializer.deserializeEntity(curLineStream, *it);
        };
    }

    /*!
     * \brief Stop reading the restart file.
     */
    void deserializeEnd()
    {
        inStream_.close();
    }

    /*!
     * \brief Returns the list of restart files in the current directory.
     */
    static void restartFileList(std::list<std::string> &fileList,
                                const std::string directory=".")
    {
        DUNE_THROW(NotImplemented,
                   "Dune::Restart::restartFileList()");
    };


private:
    std::ifstream inStream_;
    std::ofstream outStream_;
};
}

#endif
