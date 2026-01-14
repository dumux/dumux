// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_PVDREADER_HH
#define DUMUX_PVDREADER_HH

#include <iostream>
#include <iterator>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <dumux/io/xml/tinyxml2.h>
#include <dune/common/exceptions.hh>

#include <gridformat/gridformat.hpp>
#include <gridformat/grid/discontinuous.hpp>
#include <gridformat/traits/dune.hpp>

/*
 * Description: This class is responsible to store the acess nodes
 * for the vtk-files based on a given pvd file
 * \note: This should be done with use of gridformat
 */

template<class TypeTag>
class SolutionLoader
{
    public:
    // It seems not as if this class needs a ctor doing something.
    // It knows everything based on the TypeTag.
    SolutionLoader()
    { }

    // function to create the nodeslisting
    void makeSolutionVectorsListing(const std::string& problemName)
    {
        std::string fileName = problemName+".pvd";
        using namespace tinyxml2;
        XMLDocument pvdDoc;
        // try to open the given file
        const auto eResult = pvdDoc.LoadFile(fileName.c_str() );
        if( eResult != XML_SUCCESS)
            DUNE_THROW(Dune::IOError, "Couldn't open XML file " << fileName << ".");
        const XMLElement* collectionNode = getCollectionNode_(pvdDoc, fileName);
        if (collectionNode==nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'Collection' node in " << fileName << ".");

        // fill up the timesFilesList
        const tinyxml2::XMLElement* dataSet = collectionNode->FirstChildElement("DataSet");
        for (; dataSet != nullptr; dataSet = dataSet->NextSiblingElement("DataSet"))
        {
            const char *timestepText = dataSet->Attribute("timestep");
            const char *vtkfilename = dataSet->Attribute("file");
            std::cout <<"inserting node time = " << timestepText << " file " << vtkfilename << std::endl;
            nodes_.push_back(  std::make_pair(timestepText,vtkfilename) );
        }
        // update the number of Entries.
        numberOfEntries_=nodes_.size();
    }

    // This function needs a timeLoop to keep track of the time.
    template<class TimeLoop>
    std::string nextSolution(TimeLoop& timeLoop)
    {
        ++currentEntryIndex_;
        return loadFromEntry(timeLoop, currentEntryIndex_);
    }

    // A function to get a specific solution, indentified by its entry index
    template<class TimeLoop>
    std::string loadFromEntry(TimeLoop& timeLoop, int idx)
    {
        if(currentEntryIndex_>numberOfEntries_)
            return std::string{"Error: no valid entry has been found."};
        else
        {
            auto nodeData = nodes_[currentEntryIndex_];
            // As the pvd files times are not exact enough, we need to track time with use of a timeLoop.
            // This only as we have a reporting that has constant time steps
            auto checkthisval  = std::stold(nodeData.first);
            timeLoop.setTime(simTime_ );

            std::cout << " Loading node : " << checkthisval << " which is equal to: "<< simTime_ << std::endl;
            simTime_+=dt_;
            // return the file name
            return nodeData.second;
        }

    }

    // getters
    int currentEntryIndex()
    {
        return currentEntryIndex_;
    }

    // Allows to set a time stepping difference. This
    // is used as the pvd time step could be not exact enough.
    void setDt(double dt)
    {
        dt_ = dt;
    }
    double simTime()
    {
        return simTime_;
    }
    // Also a setter, just in case
    void setSimTime(double st)
    {
        simTime_=st;
    }

    private:
    // some private helpers

    // to get the collection node
    const tinyxml2::XMLElement* getCollectionNode_(const tinyxml2::XMLDocument& doc, const std::string& fileName)
    {
        using namespace tinyxml2;

        const XMLElement* collectionNode = doc.FirstChildElement("VTKFile");

        if (collectionNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'VTKFile' node in " << fileName << ".");

        return collectionNode->FirstChildElement("Collection");
    }
    std::vector<std::pair<const std::string,const std::string> > nodes_;

    int currentEntryIndex_{-1}; // this seems to be convenient in some way...
    int numberOfEntries_{0};
    // double check variables
    double simTime_{0.0};
    double dt_{0.0};
};
#endif // DUMUX_PVDREADER_HH
