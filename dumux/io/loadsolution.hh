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
#ifndef DUMUX_IO_LOADSOLUTION_HH
#define DUMUX_IO_LOADSOLUTION_HH

#include <string>
#include <iostream>
#include <vector>
#include <unordered_set>

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/io/vtk/vtureader.hh>

namespace Dumux {

namespace Detail {
//! helper struct detecting if a PrimaryVariables object has a state() function
struct hasState
{
    template<class PrimaryVariables>
    auto operator()(PrimaryVariables&& priVars)
    -> decltype(priVars.state())
    {}
};

} // end namespace Detail

/*!
 * \ingroup InputOutput
 * \brief helper function to read from a file into a solution vector
 */
template <class SolutionVector, class PvNamesFunc>
auto loadSolutionFromVtuFile(const std::string fileName,
                             const VTUReader::DataType& dataType,
                             PvNamesFunc&& pvNamesFunc,
                             SolutionVector& sol)
-> typename std::enable_if_t<!decltype(isValid(Detail::hasState())(sol[0]))::value, void>
{
    VTUReader vtu(fileName);
    using PrimaryVariables = typename SolutionVector::block_type;
    using Scalar = typename PrimaryVariables::field_type;

    for (size_t pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
    {
        const auto pvName = pvNamesFunc(pvIdx);
        const auto vec = vtu.readData<std::vector<Scalar>>(pvName, dataType);
        if (vec.size() != sol.size())
            DUNE_THROW(Dune::IOError, "Size mismatch between solution vector and read data (" << sol.size() << " != " << vec.size() << ")");

        for (std::size_t i = 0; i < sol.size(); ++i)
            sol[i][pvIdx] = vec[i];
    }
}

/*!
 * \ingroup InputOutput
 * \brief helper function to read from a file into a solution vector
 */
template <class SolutionVector, class PvNamesFunc>
auto loadSolutionFromVtuFile(const std::string fileName,
                             const VTUReader::DataType& dataType,
                             PvNamesFunc&& pvNamesFunc,
                             SolutionVector& sol)
-> typename std::enable_if_t<decltype(isValid(Detail::hasState())(sol[0]))::value, void>
{
    VTUReader vtu(fileName);
    const auto vec = vtu.readData<std::vector<int>>("phase presence", dataType);
    std::unordered_set<int> states;
    for (size_t i = 0; i < sol.size(); ++i)
    {
        const int state = vec[i];
        sol[i].setState(state);
        states.insert(state);
    }

    using PrimaryVariables = typename SolutionVector::block_type;
    using Scalar = typename PrimaryVariables::field_type;
    for (size_t pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
    {
        if (pvNamesFunc(pvIdx, 1) == pvNamesFunc(pvIdx, 2))
        {
            const auto vec = vtu.readData<std::vector<Scalar>>(pvNamesFunc(pvIdx, 1), dataType);
            for (size_t i = 0; i < sol.size(); ++i)
                sol[i][pvIdx] = vec[i];
        }
        else
        {
            std::unordered_map<int, std::vector<Scalar>> switchedPvsSol;
            for (const auto& state : states)
                switchedPvsSol[state] = vtu.readData<std::vector<Scalar>>(pvNamesFunc(pvIdx, state), dataType);

            for (size_t i = 0; i < sol.size(); ++i)
                sol[i][pvIdx] = switchedPvsSol[sol[i].state()][i];
        }
    }
}

/*!
 * \ingroup InputOutput
 * \brief helper function to determine the primray variable names of a model with privar state
 */
template<class ModelTraits, class FluidSystem>
std::string primaryVariableName(int pvIdx, int state)
{
    using std::pow;
    static auto numStates = pow(2, ModelTraits::numPhases()) - 1;
    const auto paramNameWithState = "LoadSolution.PriVarNamesState" + std::to_string(state);

    if (haveParam("LoadSolution.PriVarNames"))
    {
        DUNE_THROW(Dune::NotImplemented, "please provide LoadSolution.PriVarNamesState1..." << numStates
                   << " or remove LoadSolution.PriVarNames to use default names");
    }
    else if (haveParam(paramNameWithState))
    {
        const auto pvNames = getParam<std::vector<std::string>>(paramNameWithState);
        return pvNames[pvIdx];
    }
    else
        return ModelTraits::template primaryVariableName<FluidSystem>(pvIdx, state);
}

/*!
 * \ingroup InputOutput
 * \brief helper function to determine the primray variable names of a model
 */
template<class ModelTraits>
std::string primaryVariableName(int pvIdx)
{
    if (haveParam("LoadSolution.PriVarNames"))
    {
        static auto pvNames = getParam<std::vector<std::string>>("LoadSolution.PriVarNames");
        return pvNames[pvIdx];
    }
    else
        return ModelTraits::primaryVariableName(pvIdx);
}


/*!
 * \ingroup InputOutput
 * \brief load a solution vector from file
 * \note Supports the following file extensions: *.vtu
 */
template <class SolutionVector, class PvNamesFunc>
void loadSolution(const std::string& fileName,
                  DiscretizationMethod discMethod,
                  PvNamesFunc&& pvNamesFunc,
                  SolutionVector& sol)
{
    const auto extension = fileName.substr(fileName.find_last_of(".") + 1);

    if (extension == "vtu")
    {
        const auto dataType = discMethod == DiscretizationMethod::box
                              ? VTUReader::DataType::pointData : VTUReader::DataType::cellData;
        loadSolutionFromVtuFile(fileName, dataType, pvNamesFunc, sol);
    }
    else
        DUNE_THROW(Dune::NotImplemented, "loadSolution for extension " << extension);
}

} // namespace Dumux

#endif
