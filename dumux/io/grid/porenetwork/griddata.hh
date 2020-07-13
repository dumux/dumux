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
 * \brief Class for grid data attached to dgf or gmsh grid files
 */
#ifndef DUMUX_IO_PORENETWORKGRID_DATA_HH
#define DUMUX_IO_PORENETWORKGRID_DATA_HH

#include <algorithm>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>

// FoamGrid specific includes
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#endif

#include <dumux/common/indextraits.hh>
#include <dumux/porenetworkflow/common/geometry.hh>
#include <dumux/porenetworkflow/common/throatproperties.hh>

#include "parametersforgeneratedgrid.hh"

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Class for grid data attached to dgf or gmsh grid files
 */
template <class Grid>
class PoreNetworkGridData
{
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using Intersection = typename Grid::LeafIntersection;
    using Element = typename Grid::template Codim<0>::Entity;
    using Vertex = typename Grid::template Codim<dim>::Entity;
    using GridView = typename Grid::LeafGridView;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SmallLocalIndex = typename IndexTraits<GridView>::SmallLocalIndex;
    using StringVector = std::vector<std::string>;

    using Scalar = double;
    using PersistentParameterContainer = Dune::PersistentContainer<Grid, std::vector<typename Grid::ctype>>;
    using ParametersForGeneratedGrid = Dumux::ParametersForGeneratedGrid<Grid, Scalar>;

public:

    //! constructor for dgf grid data
    PoreNetworkGridData(Dune::GridPtr<Grid> grid, const std::string& paramGroup)
    : dgfGrid_(grid)
    , isDgfData_(true)
    , paramGroup_(paramGroup)
    , numSubregions_(0)
    {
        setParameterIndices_();
    }

    //! constructor for non-dgf grid data
    PoreNetworkGridData(std::shared_ptr<Grid> grid, const std::string& paramGroup)
    : factoryGrid_(grid)
    , isDgfData_(false)
    , paramGroup_(paramGroup)
    {
        numSubregions_ = getParamFromGroup<std::size_t>(paramGroup_, "Grid.NumSubregions", 0);
        setParameterIndices_();
        parametersForGeneratedGrid_ = std::make_unique<ParametersForGeneratedGrid>(gridView_(), paramGroup);
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for vertex data
     * \note You can only pass vertices that exist on level 0!
     */
    const std::vector<double>& parameters(const Vertex& vertex) const
    {
        if (isDgfData_ && !useCopiedDgfData_)
            return dgfGrid_.parameters(vertex);
        else
        {
            assert(!(*vertexParameters_)[vertex].empty() && "No parameters available.");
            return (*vertexParameters_)[vertex];
        }
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for element data
     */
    const std::vector<double>& parameters(const Element& element) const
    {
        if (isDgfData_ && !useCopiedDgfData_)
        {
            if (element.hasFather())
            {
                auto level0Element = element;
                while(level0Element.hasFather())
                    level0Element = level0Element.father();

                return dgfGrid_.parameters(level0Element);
            }
            else
            {
                return dgfGrid_.parameters(element);
            }
        }
        else
        {
            assert(!(*elementParameters_)[element].empty() && "No parameters available.");
            return (*elementParameters_)[element];
        }
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available
     */
    template <class GridImp, class IntersectionImp>
    const Dune::DGFBoundaryParameter::type& parameters(const Dune::Intersection<GridImp, IntersectionImp>& intersection) const
    {
        if (isDgfData_)
            return dgfGrid_.parameters(intersection);
        else
            DUNE_THROW(Dune::InvalidStateException, "The parameters method is only available if the grid was constructed with a DGF file.");
    }

    /*!
     * \brief Returns the coordination numbers for all pore bodies.
     */
    std::vector<SmallLocalIndex> getCoordinationNumbers() const
    {
        std::vector<SmallLocalIndex>  coordinationNumbers(gridView_().size(dim), 0);

        for (const auto &element : elements(gridView_()))
        {
            for (SmallLocalIndex vIdxLocal = 0; vIdxLocal < 2; ++vIdxLocal)
            {
                const auto vIdxGlobal = gridView_().indexSet().subIndex(element, vIdxLocal, dim);
                coordinationNumbers[vIdxGlobal] += 1;
            }
        }

        if (std::any_of(coordinationNumbers.begin(), coordinationNumbers.end(), [](auto i){ return i == 0; }))
            DUNE_THROW(Dune::InvalidStateException, "One of the pores is not connected to another pore. SanitizeGrid will not help in this case. Check your grid file");

        return coordinationNumbers;
    }

    /*!
     * \brief Assign parameters for generically created grids
     */
    void assignParameters()
    {
        if (isDgfData_)
            DUNE_THROW(Dune::InvalidStateException, "Assigning parameter not possible for dgf gids");

        const auto numVertexParams = vertexParameterNames_.size();
        const auto numElementParams = elementParameterNames_.size();
        vertexParameters_ = makeParamContainer_(*factoryGrid_, numVertexParams, 1/*codim*/);
        elementParameters_ = makeParamContainer_(*factoryGrid_, numElementParams, 0/*codim*/);

        auto setParamHelper = [&](const auto& entity, const std::string& param, const Scalar value)
        {
            setParameter_(entity, param, value);
        };

        auto getParamHelper = [&](const auto& entity, const std::string& param)
        {
            return getParameter(entity, param);
        };

        parametersForGeneratedGrid_->assignParameters(setParamHelper, getParamHelper, numSubregions_);
    }

    void resizeParameterContainers()
    {
        // resize the parameters
        vertexParameters_->resize();
        elementParameters_->resize();
        vertexParameters_->shrinkToFit();
        elementParameters_->shrinkToFit();
    }

    void copyDgfData()
    {
        if (!isDgfData_)
            DUNE_THROW(Dune::InvalidStateException, "copying dgf data only works when a dgf grid is actually used");

        useCopiedDgfData_ = true;
        const auto someVertex = *(vertices(gridView_()).begin());
        const auto someElement = *(elements(gridView_()).begin());
        const auto numVertexParams = dgfGrid_.parameters(someVertex).size();
        const auto numElementParams = dgfGrid_.parameters(someElement).size();
        vertexParameters_ = makeParamContainer_(*dgfGrid_, numVertexParams, 1);
        elementParameters_ = makeParamContainer_(*dgfGrid_, numElementParams, 0);

        for (const auto& element : elements(gridView_()))
        {
            for (int i = 0; i < numElementParams; ++i)
                (*elementParameters_)[element][i] = dgfGrid_.parameters(element)[i];
        }

        for (const auto& vertex : vertices(gridView_()))
        {
            for (int i = 0; i < numVertexParams; ++i)
                (*vertexParameters_)[vertex][i] = dgfGrid_.parameters(vertex)[i];
        }
    }

    /*!
     * \brief Return the index for a given parameter name
     */
    int parameterIndex(const std::string& paramName) const
    {
        // make sure the string is present in the map, throw a Dumux exception otherwise (and not a std one)
        // the [] operator can't be used here due to const correctness
        if (parameterIndex_.count(paramName))
            return parameterIndex_.at(paramName);
        else
        {
            std::string msg;
            if (paramName.find("Throat") != std::string::npos)
                msg = "Make sure to include it in the vector of parameter names ElementParameters = " + paramName + " ... ...";
            else if (paramName.find("Pore") != std::string::npos)
                msg = "Make sure to include it in the vector of parameter names VertexParameters = " + paramName + " ... ...";

            DUNE_THROW(Dumux::ParameterException, paramName << " not set in the input file. \n" << msg);
        }
    }

    /*!
     * \brief Return the parameter group
     */
    const std::string& paramGroup() const
    { return paramGroup_; }

    /*!
     * \brief Return if a given element parameter is provided by the grid
     */
    bool gridHasElementParameter(const std::string& param) const
    {
        return std::any_of(elementParameterNames_.begin(), elementParameterNames_.end(), [&param]( const auto& i ){ return (i == param); });
    }

    /*!
     * \brief Return if a given vertex parameter is provided by the grid
     */
    bool gridHasVertexParameter(const std::string& param) const
    {
        return std::any_of(vertexParameterNames_.begin(), vertexParameterNames_.end(), [&param]( const auto& i ){ return (i == param); });
    }

    /*!
     * \brief Returns the value of an element parameter
     */
    Scalar getParameter(const Element& element, const std::string& param) const
    { return (*elementParameters_)[element][parameterIndex(param)]; }

    /*!
     * \brief Returns the value of an vertex parameter
     */
    Scalar getParameter(const Vertex& vertex, const std::string& param) const
    { return (*vertexParameters_)[vertex][parameterIndex(param)]; }

    /*!
     * \brief Returns the pore label at a given position for a generic grid.
     *        This is needed by the grid creator in case not all parameters are initialized yet.
     */
    int poreLabelAtPosForGenericGrid(const GlobalPosition& pos) const
    { return parametersForGeneratedGrid_->boundaryFaceMarkerAtPos(pos); }

    /*!
     * \brief Returns the names of the vertex parameters
     */
    const std::vector<std::string>& vertexParameterNames() const
    { return vertexParameterNames_; }

    /*!
     * \brief Returns the names of the element parameters
     */
    const std::vector<std::string>& elementParameterNames() const
    { return elementParameterNames_; }

private:

    void setParameter_(const Element& element, const std::string& param, const Scalar value)
    { (*elementParameters_)[element][parameterIndex(param)] = value; }

    void setParameter_(const Vertex& vertex, const std::string& param, const Scalar value)
    { (*vertexParameters_)[vertex][parameterIndex(param)] = value; }

    void setParameterIndices_()
    {
        if (!isDgfData_)
        {
            vertexParameterNames_ = StringVector{"PoreRadius", "PoreVolume", "PoreLabel"};
            elementParameterNames_ = StringVector{"ThroatRadius", "ThroatLength", "ThroatLabel"};
            if (numSubregions_ > 0)
            {
                vertexParameterNames_.push_back("PoreRegionId");
                elementParameterNames_.push_back("ThroatRegionId");
            }
        }
        else // DGF grid data
        {
            // treat vertex parameter names
            if (auto inputFileVertexParameterNames = getParamFromGroup<StringVector>(paramGroup_, "Grid.VertexParameters", StringVector{}); !inputFileVertexParameterNames.empty())
            {
                // TODO Remove at some point
                std::cout << "\n***\nWARNING: Setting Grid.VertexParameters in the input file is deprecated. Set '% Vertex parameters: Param1 Param2 ...' in the dgf file instead\n***\n" << std::endl;
                vertexParameterNames_ = std::move(inputFileVertexParameterNames);
            }
            else
            {
                if (auto dgfFileVertexParameterNames = dgfFileParameterNames_("Vertex"); dgfFileVertexParameterNames.empty())
                    DUNE_THROW(Dune::InvalidStateException, "No vertex parameter names specified in dgf file. Set '% Vertex parameters: Param1 Param2 ...'");
                else
                    vertexParameterNames_ = std::move(dgfFileVertexParameterNames);
            }

            // treat element parameter names
            if (auto inputFileElementParameterNames = getParamFromGroup<StringVector>(paramGroup_, "Grid.ElementParameters", StringVector{}); !inputFileElementParameterNames.empty())
            {
                // TODO Remove at some point
                std::cout << "\n***\nWARNING: Setting Grid.ElementParameters in the input file is deprecated. Set '% Element parameters: Param1 Param2 ...' in the dgf file instead\n***\n" << std::endl;
                elementParameterNames_ = std::move(inputFileElementParameterNames);
            }
            else
            {
                if (auto dgfFileElementParameterNames =  dgfFileParameterNames_("Element"); dgfFileElementParameterNames.empty())
                    DUNE_THROW(Dune::InvalidStateException, "No element parameter names specified in dgf file. Set '% Element parameters: Param1 Param2 ...'");
                else
                    elementParameterNames_ = std::move(dgfFileElementParameterNames);
            }

            // make sure that the number of specified parameters matches with the dgf file
            if (const auto& someElement = *(elements(gridView_()).begin()); elementParameterNames_.size() != dgfGrid_.nofParameters(someElement))
                DUNE_THROW(Dune::InvalidStateException, "Number of user-specified element parameters (" << elementParameterNames_.size()
                            << ") does not match number of element paramters in dgf file (" << dgfGrid_.nofParameters(someElement) << ")");

            if (const auto& someVertex = *(vertices(gridView_()).begin()); vertexParameterNames_.size() != dgfGrid_.nofParameters(someVertex))
                DUNE_THROW(Dune::InvalidStateException, "Number of user-specified vertex parameters (" << vertexParameterNames_.size()
                            << ") does not match number of vertex paramters in dgf file (" << dgfGrid_.nofParameters(someVertex) << ")");
        }

        for (int i = 0; i < vertexParameterNames_.size(); ++i)
        {
            std::cout << vertexParameterNames_[i] << " is vertex parameter " << i << std::endl;
            parameterIndex_[vertexParameterNames_[i]] = i;
        }

        for (int i = 0; i < elementParameterNames_.size(); ++i)
        {
            std::cout << elementParameterNames_[i] << " is element parameter " << i << std::endl;
            parameterIndex_[elementParameterNames_[i]] = i;
        }
    }

    /*!
     * \brief Initializes and returns a container for vertex (codim dim) or element (codim 0) data
     *
     * \param grid The grid
     * \param numParams The number of paramters
     * \param codim The codimension
     */
    auto makeParamContainer_(const Grid& grid, int numParams, int codim) const
    {
        auto parameters = std::make_unique<PersistentParameterContainer>(grid, codim);
        (*parameters).resize();
        for (auto&& v : (*parameters))
            v.resize(numParams);
        return parameters;
    }

    StringVector dgfFileParameterNames_(const std::string& entity) const
    {
        std::ifstream gridFile(getParamFromGroup<std::string>(paramGroup_, "Grid.File"));
        std::string line;
        while (getline(gridFile, line))
        {
            if (line.find(entity + " parameters:", 0) != std::string::npos)
            {
                std::string args = line.substr(line.find(":")+1, std::string::npos);
                StringVector paramNames;
                std::istringstream iss(args);
                std::string item;
                while (std::getline(iss, item, ' '))
                    if (!item.empty())
                        *std::back_inserter(paramNames)++ = item;

                return paramNames;
            }
        }
        return StringVector();
    }

    /*!
     * \brief Return the gridView this grid geometry object lives on
     */
    const GridView gridView_() const
    { return isDgfData_ ? dgfGrid_->leafGridView() : factoryGrid_->leafGridView(); }

    // dgf grid data
    Dune::GridPtr<Grid> dgfGrid_;

    std::shared_ptr<Grid> factoryGrid_;

    bool isDgfData_ = false;
    bool useCopiedDgfData_ = false;
    std::string paramGroup_;

    std::unique_ptr<ParametersForGeneratedGrid> parametersForGeneratedGrid_;

    std::vector<std::string> vertexParameterNames_;
    std::vector<std::string> elementParameterNames_;

    std::unique_ptr<PersistentParameterContainer> vertexParameters_;
    std::unique_ptr<PersistentParameterContainer> elementParameters_;

    std::size_t numSubregions_;

    std::unordered_map<std::string, int> parameterIndex_;
};

} // namespace Dumux

#endif
