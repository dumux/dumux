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
 * \brief read from a file into a solution vector
 */
#ifndef DUMUX_IO_LOADSOLUTION_HH
#define DUMUX_IO_LOADSOLUTION_HH

#include <string>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <type_traits>
#include <functional>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/common/typetraits/state.hh>
#include <dumux/io/vtk/vtkreader.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief a data handle to communicate the solution on ghosts and overlaps
 *        when reading from vtk file in parallel
 */
template <class Container, class EntityMapper, int codim>
class LoadSolutionDataHandle
: public Dune::CommDataHandleIF< LoadSolutionDataHandle<Container, EntityMapper, codim>,
                                 std::decay_t<decltype(std::declval<Container>()[0])> >
{
    using FieldType = std::decay_t<decltype(std::declval<Container>()[0])>;
public:
    LoadSolutionDataHandle(Container& container,
                           const EntityMapper& mapper)
    : mapper_(mapper)
    , container_(container)
    {}

    bool contains(int dim, int cd) const
    { return cd == codim; }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedSize(int dim, int cd) const
    { return true; }

    template<class EntityType>
    std::size_t size (const EntityType &e) const
    { return 1; }

    template<class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp& buff, const EntityType& e) const
    {
        const auto vIdx = mapper_.index(e);
        buff.write(container_[vIdx]);
    }

    template<class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp& buff, const EntityType& e, std::size_t n)
    {
        const auto vIdx = mapper_.index(e);
        FieldType tmp;
        buff.read(tmp);
        container_[vIdx] = tmp;
    }

private:
    EntityMapper mapper_;
    Container& container_;
};

/*!
 * \ingroup InputOutput
 * \brief read from a vtk file into a solution vector with primary variables without state
 */
template <class SolutionVector, class PvNameFunc, class GridGeometry>
auto loadSolutionFromVtkFile(SolutionVector& sol,
                             const std::string fileName,
                             PvNameFunc&& targetPvNameFunc,
                             const GridGeometry& gridGeometry,
                             const VTKReader::DataType& dataType)
-> typename std::enable_if_t<!decltype(isValid(Detail::hasState())(sol[0]))::value, void>
{
    VTKReader vtu(fileName);

    using PrimaryVariables = typename SolutionVector::block_type;
    using Scalar = typename PrimaryVariables::field_type;
    constexpr auto dim = GridGeometry::GridView::dimension;
    const std::size_t targetSolutionSize = PrimaryVariables::dimension;

    std::size_t matchingLoadedArrays = 0;
    for (std::size_t i = 0; i < targetSolutionSize; i++)
        if (vtu.hasData(targetPvNameFunc(i,0), dataType))
            matchingLoadedArrays++;

    if (matchingLoadedArrays < targetSolutionSize)
        std::cout << "The loaded solution does not provide a data array for each of the primary variables. \n"
                  << "The target solution has "<< targetSolutionSize << " entries, "
                  << "whereas the loaded solution provides only " << matchingLoadedArrays << " data array(s). \n"
                  << "Make sure that the model concepts are compatible, "
                  << "and be sure to provide initial conditions for the missing primary variables. \n";

    for (std::size_t targetPvIdx = 0; targetPvIdx < targetSolutionSize; ++targetPvIdx)
    {
        std::vector<Scalar> vec;
        const auto targetPvName = targetPvNameFunc(targetPvIdx, 0);

        if (vtu.hasData(targetPvName, dataType))
            vec = vtu.readData<std::vector<Scalar>>(targetPvName, dataType);
        else
        {
            std::cout << "The loaded solution does not have a field named \"" << targetPvName << "\". "
                      << "Make sure this field is filled using the initial method in the problem definition. \n";
            continue;
        }

        if (dataType == VTKReader::DataType::cellData)
        {
            std::size_t i = 0;
            for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
            {
                const auto eIdx = gridGeometry.elementMapper().index(element);
                sol[eIdx][targetPvIdx] = vec[i++];
            }
        }
        // for staggered face data (which is written out as VTK point data) we just read in the vector
        else if (dataType == VTKReader::DataType::pointData && GridGeometry::discMethod == DiscretizationMethod::staggered)
        {
            if (sol.size() != vec.size())
                DUNE_THROW(Dune::InvalidStateException, "Solution size (" << sol.size() << ") does not match input size (" << vec.size() << ")!");

            for (std::size_t i = 0; i < sol.size(); ++i)
                sol[i][targetPvIdx] = vec[i];
        }
        else
        {
            std::size_t i = 0;
            std::vector<bool> visited(gridGeometry.gridView().size(dim), false);
            for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
            {
                for (int vIdxLocal = 0; vIdxLocal < element.subEntities(dim); ++vIdxLocal)
                {
                    const auto vIdxGlobal = gridGeometry.vertexMapper().subIndex(element, vIdxLocal, dim);
                    if (!visited[vIdxGlobal])
                    {
                        sol[vIdxGlobal][targetPvIdx] = vec[i++];
                        visited[vIdxGlobal] = true;
                    }
                }
            }
        }
    }
}

/*!
 * \ingroup InputOutput
 * \brief read from a sequential file into a solution vector with primary variables with state
 */
template <class SolutionVector, class PvNameFunc, class GridGeometry>
auto loadSolutionFromVtkFile(SolutionVector& sol,
                             const std::string fileName,
                             PvNameFunc&& targetPvNameFunc,
                             const GridGeometry& gridGeometry,
                             const VTKReader::DataType& dataType)
-> typename std::enable_if_t<decltype(isValid(Detail::hasState())(sol[0]))::value, void>
{
    VTKReader vtu(fileName);
    // get states at each dof location
    const auto stateAtDof = vtu.readData<std::vector<int>>("phase presence", dataType);

    // determine all states that are present
    std::unordered_set<int> states;
    for (std::size_t i = 0; i < stateAtDof.size(); ++i)
        states.insert(stateAtDof[i]);

    using PrimaryVariables = typename SolutionVector::block_type;
    using Scalar = typename PrimaryVariables::field_type;
    const std::size_t targetSolutionSize = PrimaryVariables::dimension;

    std::unordered_set<std::string> matchingNames;
    for (std::size_t i = 0; i < targetSolutionSize; i++)
        for (const auto& state : states)
            if ( vtu.hasData(targetPvNameFunc(i,state), dataType))
                matchingNames.insert(targetPvNameFunc(i,state));

    const std::size_t matchingLoadedArrays = matchingNames.size() - (states.size()-1);

    if (matchingLoadedArrays < targetSolutionSize)
        std::cout << "The loaded solution does not provide a data array for each of the primary variables. \n"
                  << "The target solution has "<< targetSolutionSize << " entries, "
                  << "whereas the loaded solution provides only " << matchingLoadedArrays << " data array(s). \n"
                  << "Make sure that the model concepts are compatible, "
                  << "and be sure to provide initial conditions for the missing primary variables. \n";

    for (std::size_t targetPvIdx = 0; targetPvIdx < targetSolutionSize; ++targetPvIdx)
    {
        std::unordered_map<int, std::vector<Scalar>> data;
        for (const auto& state : states)
        {
            const auto targetPvName = targetPvNameFunc(targetPvIdx, state);

            if (vtu.hasData(targetPvName, dataType))
                data[state] = vtu.readData<std::vector<Scalar>>(targetPvName, dataType);
            else
            {
                std::cout << "Loaded Solution does not have a field named \"" << targetPvName << "\". "
                          << "Make sure this field is filled using the initial method in the problem definition. \n";
                continue;
            }
        }

        if (dataType == VTKReader::DataType::cellData)
        {
            std::size_t i = 0;
            for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
            {
                const auto eIdx = gridGeometry.elementMapper().index(element);
                const auto state = stateAtDof[i];
                sol[eIdx][targetPvIdx] = data[state][i++];
                sol[eIdx].setState(state);
            }
        }
        else
        {
            std::size_t i = 0;
            constexpr int dim = GridGeometry::GridView::dimension;
            std::vector<bool> visited(gridGeometry.gridView().size(dim), false);
            for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
            {
                for (int vIdxLocal = 0; vIdxLocal < element.subEntities(dim); ++vIdxLocal)
                {
                    const auto vIdxGlobal = gridGeometry.vertexMapper().subIndex(element, vIdxLocal, dim);
                    if (!visited[vIdxGlobal])
                    {
                        const auto state = stateAtDof[i];
                        sol[vIdxGlobal][targetPvIdx] = data[state][i++];
                        sol[vIdxGlobal].setState(state);
                        visited[vIdxGlobal] = true;
                    }
                }
            }
        }
    }
}

/*!
 * \ingroup InputOutput
 * \brief helper function to determine the primary variable names of a model with privar state
 * \note use this as input for the load solution function
 */
template<class IOFields, class PrimaryVariables, class ModelTraits = void, class FluidSystem = void, class SolidSystem = void>
auto createPVNameFunction(const std::string& paramGroup = "")
-> typename std::enable_if_t<decltype(isValid(Detail::hasState())(PrimaryVariables(0)))::value, std::function<std::string(int,int)>>
{
    return  [paramGroup](int pvIdx, int state = 0)
            {
                static auto numStates = (1 << ModelTraits::numFluidPhases()) - 1;
                const auto paramNameWithState = "LoadSolution.PriVarNamesState" + std::to_string(state);

                if (hasParamInGroup(paramGroup, "LoadSolution.PriVarNames") && !hasParamInGroup(paramGroup, paramNameWithState))
                    DUNE_THROW(Dune::NotImplemented, "please provide LoadSolution.PriVarNamesState1..." << numStates
                              << " or remove LoadSolution.PriVarNames to use the model's default primary variable names");

                else if (hasParamInGroup(paramGroup, paramNameWithState))
                {
                    const auto pvName = getParamFromGroup<std::vector<std::string>>(paramGroup, paramNameWithState);
                    return pvName[pvIdx];
                }
                else
                    return IOFields::template primaryVariableName<ModelTraits, FluidSystem, SolidSystem>(pvIdx, state);
            };
}

/*!
 * \ingroup InputOutput
 * \brief helper function to determine the primary variable names of a model without state
 * \note use this as input for the load solution function
 */
template<class IOFields, class PrimaryVariables, class ModelTraits = void, class FluidSystem = void, class SolidSystem = void>
auto createPVNameFunction(const std::string& paramGroup = "")
-> typename std::enable_if_t<!decltype(isValid(Detail::hasState())(PrimaryVariables(0)))::value, std::function<std::string(int,int)>>
{
    if (hasParamInGroup(paramGroup, "LoadSolution.PriVarNames"))
    {
        const auto pvName = getParamFromGroup<std::vector<std::string>>(paramGroup, "LoadSolution.PriVarNames");
        return [n = std::move(pvName)](int pvIdx, int state = 0){ return n[pvIdx]; };
    }
    else
        return [](int pvIdx, int state = 0){ return IOFields::template primaryVariableName<ModelTraits, FluidSystem, SolidSystem>(pvIdx, state); };
}

/*!
 * \ingroup InputOutput
 * \brief load a solution vector from file
 * \note Supports the following file extensions: *.vtu *.vtp *.pvtu, *.pvtp
 * \param sol the solution vector to read from file
 * \param fileName the file name of the file to read from
 * \param targetPvNameFunc a function with the signature std::string(int pvIdx)
 *        in case the primary variables have a state the signature is std::string(int pvIdx, int state)
 * \param gridGeometry the grid geometry of the discretization method used
 */
template <class SolutionVector, class PvNameFunc, class GridGeometry>
void loadSolution(SolutionVector& sol,
                  const std::string& fileName,
                  PvNameFunc&& targetPvNameFunc,
                  const GridGeometry& gridGeometry)
{
    const auto extension = fileName.substr(fileName.find_last_of(".") + 1);
    auto dataType = GridGeometry::discMethod == DiscretizationMethod::box ?
                    VTKReader::DataType::pointData : VTKReader::DataType::cellData;

    if (extension == "vtu" || extension == "vtp")
    {
        if (GridGeometry::discMethod == DiscretizationMethod::staggered && extension == "vtp")
            dataType = VTKReader::DataType::pointData;

        loadSolutionFromVtkFile(sol, fileName, targetPvNameFunc, gridGeometry, dataType);
    }
    else if (extension == "pvtu" || extension == "pvtp")
    {
        if (GridGeometry::discMethod == DiscretizationMethod::staggered)
            DUNE_THROW(Dune::NotImplemented, "reading staggered solution from a parallel vtk file");

        loadSolutionFromVtkFile(sol, fileName, targetPvNameFunc, gridGeometry, dataType);
    }
    else
        DUNE_THROW(Dune::NotImplemented, "loadSolution for file with extension " << extension);

    // communicate solution on ghost and overlap dofs
    if (gridGeometry.gridView().comm().size() > 1)
    {
        using GridView = typename GridGeometry::GridView;
        if (dataType == VTKReader::DataType::cellData)
        {
            LoadSolutionDataHandle<SolutionVector, typename GridGeometry::ElementMapper, 0>
                dataHandle(sol, gridGeometry.elementMapper());
            gridGeometry.gridView().communicate(dataHandle,
                                                  Dune::InteriorBorder_All_Interface,
                                                  Dune::ForwardCommunication);
        }
        else
        {
            LoadSolutionDataHandle<SolutionVector, typename GridGeometry::VertexMapper, GridView::dimension>
                dataHandle(sol, gridGeometry.vertexMapper());
            gridGeometry.gridView().communicate(dataHandle,
                                                  Dune::InteriorBorder_All_Interface,
                                                  Dune::ForwardCommunication);
        }
    }
}
} // end namespace Dumux

#endif
