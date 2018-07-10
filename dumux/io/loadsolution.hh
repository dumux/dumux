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
 * \ingroup InputOutput
 * \brief read from a file into a solution vector
 */
#ifndef DUMUX_LOADSOLUTION_HH
#define DUMUX_LOADSOLUTION_HH

#include <string>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/io/tinyxml2/tinyxml2.h>

namespace Dumux {

//! helper struct detecting if a PrimaryVariables object has a state() function
struct hasState
{
    template<class PrimaryVariables>
    auto operator()(PrimaryVariables&& priVars)
    -> decltype(priVars.state())
    {}
};

/*!
 * \ingroup InputOutput
 * \brief helper class to read from a file into a solution vector
 */
class LoadSolutionHelper {
public:
    template <class FVGridGeometry, class SolutionVector, class PvNamesFunc>
    static void loadFromVtkFile(const std::string fileName,
                                const FVGridGeometry& fvGridGeometry,
                                PvNamesFunc pvNamesFunc,
                                SolutionVector& sol)
    {
        using namespace tinyxml2;

        XMLDocument xmlDoc;
        auto eResult = xmlDoc.LoadFile(fileName.c_str());
        if (eResult != XML_SUCCESS)
            DUNE_THROW(Dune::IOError, "Couldn't open XML file " << fileName << ".");

        XMLElement *pieceNode = xmlDoc.FirstChildElement("VTKFile")->FirstChildElement("UnstructuredGrid")->FirstChildElement("Piece");
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get Piece node in " << fileName << ".");

        XMLElement *dataNode = nullptr;
        if (FVGridGeometry::discMethod == DiscretizationMethod::box)
            dataNode = pieceNode->FirstChildElement("PointData");
        else
            dataNode = pieceNode->FirstChildElement("CellData");

        if (dataNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get data node in " << fileName << ".");

        setPrimaryVariables_(dataNode, fvGridGeometry, pvNamesFunc, sol);
    }

private:

    template <class Scalar, class FVGridGeometry>
    static std::vector<Scalar> extractDataToVector_(tinyxml2::XMLElement *xmlNode,
                                                    const std::string& name,
                                                    const FVGridGeometry& fvGridGeometry)
    {
        // loop over XML node siblings to find the correct data array
        tinyxml2::XMLElement *dataArray = xmlNode->FirstChildElement("DataArray");
        for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
        {
            const char *attributeText = dataArray->Attribute("Name");

            if (attributeText == nullptr)
                DUNE_THROW(Dune::IOError, "Couldn't get Name attribute.");

            if (std::string(attributeText) == name)
                break;
        }
        if (dataArray == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't find the data array " << name << ".");

        std::stringstream dataStream(dataArray->GetText());
        std::vector<Scalar> vec(fvGridGeometry.numDofs());
        for (auto& val : vec)
        {
            dataStream >> val;
        }

        return vec;
    }

    template<class SolutionVector, class FVGridGeometry, class PvNamesFunc>
    static auto setPrimaryVariables_(tinyxml2::XMLElement *node,
                                     const FVGridGeometry& fvGridGeometry,
                                     PvNamesFunc pvNamesFunc,
                                     SolutionVector& sol)
    -> typename std::enable_if_t<!decltype(isValid(hasState())(sol[0]))::value, void>
    {
        using PrimaryVariables = typename SolutionVector::block_type;
        using Scalar = typename PrimaryVariables::field_type;

        for (size_t pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
        {
            auto pvName = pvNamesFunc(pvIdx);

            auto vec = extractDataToVector_<Scalar>(node, pvName, fvGridGeometry);

            for (size_t i = 0; i < sol.size(); ++i)
                sol[i][pvIdx] = vec[i];
        }
    }

    template<class SolutionVector, class FVGridGeometry, class PvNamesFunc>
    static auto setPrimaryVariables_(tinyxml2::XMLElement *node,
                                     const FVGridGeometry& fvGridGeometry,
                                     PvNamesFunc pvNamesFunc,
                                     SolutionVector& sol)
    -> typename std::enable_if_t<decltype(isValid(hasState())(sol[0]))::value, void>
    {
        auto vec = extractDataToVector_<int>(node, "phase presence", fvGridGeometry);

        int minState = std::numeric_limits<int>::max();
        int maxState = std::numeric_limits<int>::min();
        for (size_t i = 0; i < sol.size(); ++i)
        {
            int state = vec[i];
            sol[i].setState(state);
            using std::min;
            minState = min(minState, state);
            using std::max;
            maxState = max(maxState, state);
        }

        using PrimaryVariables = typename SolutionVector::block_type;
        using Scalar = typename PrimaryVariables::field_type;
        for (size_t pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
        {
            if (pvNamesFunc(pvIdx, 1) == pvNamesFunc(pvIdx, 2))
            {
                auto vec = extractDataToVector_<Scalar>(node, pvNamesFunc(pvIdx, 1), fvGridGeometry);

                for (size_t i = 0; i < sol.size(); ++i)
                    sol[i][pvIdx] = vec[i];
            }
            else
            {
                std::vector<std::vector<Scalar>> switchedPvsSol;
                for (int state = minState; state <= maxState; ++state)
                    switchedPvsSol.push_back(extractDataToVector_<Scalar>(node,
                                                                          pvNamesFunc(pvIdx, state),
                                                                          fvGridGeometry));

                for (size_t i = 0; i < sol.size(); ++i)
                    sol[i][pvIdx] = switchedPvsSol[sol[i].state() - minState][i];
            }
        }
    }
};

template<class ModelTraits, class FluidSystem>
std::string primaryVariableName(int pvIdx, int state)
{
    if (haveParam("LoadSolution.PrimaryVariableNames"))
        DUNE_THROW(Dune::NotImplemented, "reading PrimaryVariableNames from input file");
    else
        return ModelTraits::template primaryVariableName<FluidSystem>(pvIdx, state);
}

template<class ModelTraits>
std::string primaryVariableName(int pvIdx)
{
    if (haveParam("LoadSolution.PrimaryVariableNames"))
        DUNE_THROW(Dune::NotImplemented, "reading PrimaryVariableNames from input file");
    else
        return ModelTraits::primaryVariableName(pvIdx);
}

template <class FVGridGeometry, class SolutionVector, class PvNamesFunc>
void loadSolution(const std::string& fileName,
                  const FVGridGeometry& fvGridGeometry,
                  PvNamesFunc pvNamesFunc,
                  SolutionVector& sol)
{
    auto extension = fileName.substr(fileName.find_last_of(".") + 1);

    if (extension == "vtu")
        LoadSolutionHelper::loadFromVtkFile(fileName, fvGridGeometry, pvNamesFunc, sol);
    else
        DUNE_THROW(Dune::NotImplemented, "loadSolution for extension " << extension);
}

} // namespace Dumux

#endif
