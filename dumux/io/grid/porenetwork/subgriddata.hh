// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Class for sub grid data attached to dgf or gmsh grid files
 */
#ifndef DUMUX_IO_GRID_PORENETWORK_SUBGRID_DATA_HH
#define DUMUX_IO_GRID_PORENETWORK_SUBGRID_DATA_HH

#include <vector>
#include <memory>
#include "griddata.hh"

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModels
 * \brief wrapper for subgrid data
 */
template<class HostGrid, class SubGrid>
class SubGridData
{
    static constexpr int dim = SubGrid::dimension;
    static constexpr int dimWorld = SubGrid::dimensionworld;
    using Intersection = typename SubGrid::LeafIntersection;
    using Element = typename SubGrid::template Codim<0>::Entity;
    using Vertex = typename SubGrid::template Codim<dim>::Entity;
    using GridView = typename SubGrid::LeafGridView;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SmallLocalIndex = typename IndexTraits<GridView>::SmallLocalIndex;

public:
    SubGridData(const SubGrid& subGrid,
                std::shared_ptr<const GridData<HostGrid>> hostGridData)
    : subGrid_(subGrid)
    , hostGridData_(hostGridData)
    {}

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for vertex data
     * \note You can only pass vertices that exist on level 0!
     */
    const std::vector<double>& parameters(const Vertex& vertex) const
    { return hostGridData_->parameters(vertex.impl().hostEntity()); }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for element data
     */
    const std::vector<double>& parameters(const Element& element) const
    { return hostGridData_->parameters(element.impl().hostEntity()); }

    /*!
     * \brief Returns the value of an element parameter
     */
    auto getParameter(const Element& element, const std::string& param) const
    { return hostGridData_->getParameter(element.impl().hostEntity(), param); }

    /*!
     * \brief Returns the value of an vertex parameter
     */
    auto getParameter(const Vertex& vertex, const std::string& param) const
    { return hostGridData_->getParameter(vertex.impl().hostEntity(), param); }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available
     */
    template <class GridImp, class IntersectionImp>
    const Dune::DGFBoundaryParameter::type& parameters(const Dune::Intersection<GridImp, IntersectionImp>& intersection) const
    { return hostGridData_->parameters(intersection); }

    /*!
     * \brief Computes and returns the label of a given throat
     *
     * \param element The element (throat)
     */
    auto throatLabel(const Element& element) const
    { return hostGridData_->throatLabel(element.impl().hostEntity()); }

    /*!
     * \brief Returns the boundary face marker index at given position
     *
     * \param pos The current position
     */
    int boundaryFaceMarkerAtPos(const GlobalPosition& pos) const
    { return hostGridData_->boundaryFaceMarkerAtPos(pos); }

    /*!
     * \brief Returns the coordination numbers for all pore bodies.
     */
    std::vector<SmallLocalIndex> getCoordinationNumbers() const
    {
        const auto gridView = subGrid_.leafGridView();

        std::vector<SmallLocalIndex>  coordinationNumbers(gridView.size(dim), 0);

        for (const auto &element : elements(gridView))
        {
            for (SmallLocalIndex vIdxLocal = 0; vIdxLocal < 2; ++vIdxLocal)
            {
                const auto vIdxGlobal = gridView.indexSet().subIndex(element, vIdxLocal, dim);
                coordinationNumbers[vIdxGlobal] += 1;
            }
        }

        if (std::any_of(coordinationNumbers.begin(), coordinationNumbers.end(), [](auto i){ return i == 0; }))
            DUNE_THROW(Dune::InvalidStateException, "One of the pores is not connected to another pore. SanitizeGrid will not help in this case. Check your grid file");

        return coordinationNumbers;
    }

    /*!
     * \brief Return the index for a given parameter name
     */
    int parameterIndex(const std::string& paramName) const
    { return hostGridData_->parameterIndex(paramName); }

    /*!
     * \brief Return the parameter group
     * \todo For now we don't support two different parameter groups for the subgrids
     */
    const std::string& paramGroup() const
    { return hostGridData_->paramGroup(); }

    /*!
     * \brief Return if a given element parameter is provided by the grid
     */
    bool gridHasElementParameter(const std::string& param) const
    { return hostGridData_->gridHasElementParameter(param); }

    /*!
     * \brief Return if a given vertex parameter is provided by the grid
     */
    bool gridHasVertexParameter(const std::string& param) const
    { return hostGridData_->gridHasVertexParameter(param); }

    /*!
     * \brief Returns the names of the vertex parameters
     */
    const std::vector<std::string>& vertexParameterNames() const
    { return hostGridData_->vertexParameterNames() ; }

    /*!
     * \brief Returns the names of the element parameters
     */
    const std::vector<std::string>& elementParameterNames() const
    { return hostGridData_->elementParameterNames() ; }

private:
    const SubGrid& subGrid_;
    std::shared_ptr<const GridData<HostGrid>> hostGridData_;
};

} // end namespace Dumux::PoreNetwork

#endif
