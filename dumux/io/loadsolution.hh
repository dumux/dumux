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
#include <unordered_map>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/io/vtk/vtkreader.hh>

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

    bool fixedsize(int dim, int cd) const
    { return true; }

    template<class EntityType>
    size_t size (const EntityType &e) const
    { return 1; }

    template<class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp& buff, const EntityType& e) const
    {
        const auto vIdx = mapper_.index(e);
        buff.write(container_[vIdx]);
    }

    template<class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp& buff, const EntityType& e, size_t n)
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
template <class SolutionVector, class PvNameFunc, class FVGridGeometry>
auto loadSolutionFromVtkFile(SolutionVector& sol,
                             const std::string fileName,
                             PvNameFunc&& pvNameFunc,
                             const FVGridGeometry& fvGridGeometry,
                             const VTKReader::DataType& dataType)
-> typename std::enable_if_t<!decltype(isValid(Detail::hasState())(sol[0]))::value, void>
{
    VTKReader vtu(fileName);
    using PrimaryVariables = typename SolutionVector::block_type;
    using Scalar = typename PrimaryVariables::field_type;
    constexpr auto dim = FVGridGeometry::GridView::dimension;

    for (size_t pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
    {
        const auto pvName = pvNameFunc(pvIdx);
        auto vec = vtu.readData<std::vector<Scalar>>(pvName, dataType);
        if (vec.size() != sol.size())
            DUNE_THROW(Dune::IOError, "Size mismatch between solution vector and read data (" << sol.size() << " != " << vec.size() << ")");

        if (dataType == VTKReader::DataType::cellData)
        {
            std::size_t i = 0;
            for (const auto& element : elements(fvGridGeometry.gridView(), Dune::Partitions::interior))
            {
                const auto eIdx = fvGridGeometry.elementMapper().index(element);
                sol[eIdx][pvIdx] = vec[i++];
            }
        }
        else
        {
            std::size_t i = 0;
            std::vector<bool> visited(fvGridGeometry.gridView().size(dim), false);
            for (const auto& element : elements(fvGridGeometry.gridView(), Dune::Partitions::interior))
            {
                for (int vIdxLocal = 0; vIdxLocal < element.subEntities(dim); ++vIdxLocal)
                {
                    const auto vIdxGlobal = fvGridGeometry.vertexMapper().subIndex(element, vIdxLocal, dim);
                    if (!visited[vIdxGlobal])
                    {
                        sol[vIdxGlobal][pvIdx] = vec[i++];
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
template <class SolutionVector, class PvNameFunc, class FVGridGeometry>
auto loadSolutionFromVtkFile(SolutionVector& sol,
                             const std::string fileName,
                             PvNameFunc&& pvNameFunc,
                             const FVGridGeometry& fvGridGeometry,
                             const VTKReader::DataType& dataType)
-> typename std::enable_if_t<decltype(isValid(Detail::hasState())(sol[0]))::value, void>
{
    VTKReader vtu(fileName);

    // get states at each dof location
    const auto stateAtDof = vtu.readData<std::vector<int>>("phase presence", dataType);

    // determine all states that are present
    std::unordered_set<int> states;
    for (size_t i = 0; i < stateAtDof.size(); ++i)
        states.insert(stateAtDof[i]);

    using PrimaryVariables = typename SolutionVector::block_type;
    using Scalar = typename PrimaryVariables::field_type;
    for (size_t pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
    {
        // check if the primary variable is state invariant
        bool isStateInvariant = true;
        for (const auto& state : states)
            isStateInvariant = isStateInvariant && pvNameFunc(pvIdx, state) == pvNameFunc(pvIdx, *states.begin());

        std::unordered_map<int, std::vector<Scalar>> data;
        if (isStateInvariant)
            data[0] = vtu.readData<std::vector<Scalar>>(pvNameFunc(pvIdx, *states.begin()), dataType);

        // the primary variable changes with the state
        // read the data for all occuring states from the file
        else
            for (const auto& state : states)
                data[state] = vtu.readData<std::vector<Scalar>>(pvNameFunc(pvIdx, state), dataType);

        // sanity check
        for (const auto& d : data)
            if (d.second.size() != sol.size())
                DUNE_THROW(Dune::IOError, "Size mismatch between solution vector and read data (" << sol.size() << " != " << d.second.size() << ")");

        if (dataType == VTKReader::DataType::cellData)
        {
            std::size_t i = 0;
            for (const auto& element : elements(fvGridGeometry.gridView(), Dune::Partitions::interior))
            {
                const auto eIdx = fvGridGeometry.elementMapper().index(element);
                const auto state = isStateInvariant ? 0 : stateAtDof[i];
                sol[eIdx][pvIdx] = data[state][i++];
                sol[eIdx][pvIdx].setState(state);
            }
        }
        else
        {
            std::size_t i = 0;
            constexpr int dim = FVGridGeometry::GridView::dimension;
            std::vector<bool> visited(fvGridGeometry.gridView().size(dim), false);
            for (const auto& element : elements(fvGridGeometry.gridView(), Dune::Partitions::interior))
            {
                for (int vIdxLocal = 0; vIdxLocal < element.subEntities(dim); ++vIdxLocal)
                {
                    const auto vIdxGlobal = fvGridGeometry.vertexMapper().subIndex(element, vIdxLocal, dim);
                    if (!visited[vIdxGlobal])
                    {
                        const auto state = isStateInvariant ? 0 : stateAtDof[i];
                        sol[vIdxGlobal][pvIdx] = data[state][i++];
                        sol[vIdxGlobal][pvIdx].setState(state);
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
template<class ModelTraits, class FluidSystem>
std::string pvNameWithState(int pvIdx, int state, const std::string& paramGroup = "")
{
    static auto numStates = (1 << ModelTraits::numPhases()) - 1;
    const auto paramNameWithState = "LoadSolution.PriVarNamesState" + std::to_string(state);

    if (hasParamInGroup(paramGroup, "LoadSolution.PriVarNames") && !hasParamInGroup(paramGroup, paramNameWithState))
        DUNE_THROW(Dune::NotImplemented, "please provide LoadSolution.PriVarNamesState1..." << numStates
                   << " or remove LoadSolution.PriVarNames to use the model's default primary variable names");

    else if (hasParamInGroup(paramGroup, paramNameWithState))
    {
        static const auto pvNames = getParamFromGroup<std::vector<std::string>>(paramGroup, paramNameWithState);
        return pvNames[pvIdx];
    }

    else
        return ModelTraits::template primaryVariableName<FluidSystem>(pvIdx, state);
}

/*!
 * \ingroup InputOutput
 * \brief helper function to determine the primary variable names of a model without state
 * \note use this as input for the load solution function
 */
template<class ModelTraits>
std::string pvName(int pvIdx, int state, const std::string& paramGroup = "")
{
    if (hasParamInGroup(paramGroup, "LoadSolution.PriVarNames"))
    {
        static const auto pvNames = getParamFromGroup<std::vector<std::string>>(paramGroup, "LoadSolution.PriVarNames");
        return pvNames[pvIdx];
    }
    else
        return ModelTraits::primaryVariableName(pvIdx, 0);
}

/*!
 * \ingroup InputOutput
 * \brief load a solution vector from file
 * \note Supports the following file extensions: *.vtu *.vtp *.pvtu, *.pvtp
 * \param sol the solution vector to read from file
 * \param fileName the file name of the file to read from
 * \param pvNameFunc a function with the signature std::string(int pvIdx)
 *        in case the primary variables have a state the signature is std::string(int pvIdx, int state)
 * \param fvGridGeometry the grid geometry of the discretization method used
 */
template <class SolutionVector, class PvNameFunc, class FVGridGeometry>
void loadSolution(SolutionVector& sol,
                  const std::string& fileName,
                  PvNameFunc&& pvNameFunc,
                  const FVGridGeometry& fvGridGeometry)
{
    const auto extension = fileName.substr(fileName.find_last_of(".") + 1);
    auto dataType = FVGridGeometry::discMethod == DiscretizationMethod::box ?
                    VTKReader::DataType::pointData : VTKReader::DataType::cellData;

    if (extension == "vtu" || extension == "vtp")
    {
        if (FVGridGeometry::discMethod == DiscretizationMethod::staggered && extension == "vtp")
            dataType = VTKReader::DataType::pointData;

        loadSolutionFromVtkFile(sol, fileName, pvNameFunc, fvGridGeometry, dataType);
    }
    else if (extension == "pvtu" || extension == "pvtp")
    {
        if (FVGridGeometry::discMethod == DiscretizationMethod::staggered)
            DUNE_THROW(Dune::NotImplemented, "reading staggered solution from a parallel vtk file");

        loadSolutionFromVtkFile(sol, fileName, pvNameFunc, fvGridGeometry, dataType);
    }
    else
        DUNE_THROW(Dune::NotImplemented, "loadSolution for file with extension " << extension);

    // communicate solution on ghost and overlap dofs
    if (fvGridGeometry.gridView().comm().size() > 1)
    {
        using GridView = typename FVGridGeometry::GridView;
        if (dataType == VTKReader::DataType::cellData)
        {
            LoadSolutionDataHandle<SolutionVector, typename FVGridGeometry::ElementMapper, 0>
                dataHandle(sol, fvGridGeometry.elementMapper());
            fvGridGeometry.gridView().communicate(dataHandle,
                                                  Dune::InteriorBorder_All_Interface,
                                                  Dune::ForwardCommunication);
        }
        else
        {
            LoadSolutionDataHandle<SolutionVector, typename FVGridGeometry::VertexMapper, GridView::dimension>
                dataHandle(sol, fvGridGeometry.vertexMapper());
            fvGridGeometry.gridView().communicate(dataHandle,
                                                  Dune::InteriorBorder_All_Interface,
                                                  Dune::ForwardCommunication);
        }
    }
}

} // end namespace Dumux

#endif
