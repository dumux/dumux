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
#ifndef DUMUX_IO_GRID_DATA_HH
#define DUMUX_IO_GRID_DATA_HH

#include <vector>
#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dumux/io/vtk/vtkreader.hh>

// UGGrid specific includes
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include "gmshgriddatahandle.hh"

namespace Dumux {

namespace Detail {

template<class Grid>
struct isUG : public std::false_type {};

#if HAVE_UG
template<int dim>
struct isUG<Dune::UGGrid<dim>> : public std::true_type {};
#endif

} // end namespace Details

/*!
 * \ingroup InputOutput
 * \brief Class for grid data attached to dgf or gmsh grid files
 */
template <class Grid>
class GridData
{
    using Intersection = typename Grid::LeafIntersection;
    using Element = typename Grid::template Codim<0>::Entity;
    using Vertex = typename Grid::template Codim<Grid::dimension>::Entity;
    using DataHandle = GmshGridDataHandle<Grid, Dune::GridFactory<Grid>, std::vector<int>>;

    enum DataSourceType { dgf, gmsh, vtk };

public:
    //! constructor for gmsh grid data
    GridData(std::shared_ptr<Grid> grid, std::shared_ptr<Dune::GridFactory<Grid>> factory,
             std::vector<int>&& elementMarkers, std::vector<int>&& boundaryMarkers, std::vector<int>&& faceMarkers = std::vector<int>{})
    : gridPtr_(grid)
    , gridFactory_(factory)
    , elementMarkers_(elementMarkers)
    , boundaryMarkers_(boundaryMarkers)
    , faceMarkers_(faceMarkers)
    , dataSourceType_(DataSourceType::gmsh)
    {}

    //! constructor for dgf grid data
    GridData(Dune::GridPtr<Grid> grid)
    : dgfGrid_(grid)
    , dataSourceType_(DataSourceType::dgf)
    {}

    //! constructor for gmsh grid data
    GridData(std::shared_ptr<Grid> grid, std::shared_ptr<Dune::GridFactory<Grid>> factory,
             VTKReader::Data&& cellData, VTKReader::Data&& pointData)
    : gridPtr_(grid)
    , gridFactory_(factory)
    , cellData_(cellData)
    , pointData_(pointData)
    , dataSourceType_(DataSourceType::vtk)
    {}


    /*!
     * \name DGF interface functions
     */
    // \{

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for vertex data
     * \note You can only pass vertices that exist on level 0!
     */
    const std::vector<double>& parameters(const Vertex& vertex) const
    {
        if (dataSourceType_ == DataSourceType::dgf)
        {
            if (vertex.level() != 0)
                DUNE_THROW(Dune::IOError, "You can only obtain parameters for level 0 vertices!");

            return dgfGrid_.parameters(vertex);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "The parameters method is only available if the grid was constructed with a DGF file.");
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for element data
     */
    const std::vector<double>& parameters(const Element& element) const
    {
        if (dataSourceType_ == DataSourceType::dgf)
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
            DUNE_THROW(Dune::InvalidStateException, "The parameters method is only available if the grid was constructed with a DGF file.");
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available
     */
    template <class GridImp, class IntersectionImp>
    const Dune::DGFBoundaryParameter::type& parameters(const Dune::Intersection<GridImp, IntersectionImp>& intersection) const
    {
        if (dataSourceType_ == DataSourceType::dgf)
            return dgfGrid_.parameters(intersection);
        else
            DUNE_THROW(Dune::InvalidStateException, "The parameters method is only available if the grid was constructed with a DGF file.");
    }

    // \}

    /*!
     * \name Gmsh interface functions
     */
    // \{

    /*!
     * \brief Return the boundary domain marker (Gmsh physical entity number) of an intersection
              Only available when using Gmsh with GridParameterGroup.DomainMarkers = 1.
     * \param boundarySegmentIndex The boundary segment index of the intersection (intersection.boundarySegmentIndex()
     */
    int getBoundaryDomainMarker(int boundarySegmentIndex) const
    {
        if (dataSourceType_ != DataSourceType::gmsh)
            DUNE_THROW(Dune::InvalidStateException, "Domain markers are only available for gmsh grids.");
        if (boundarySegmentIndex >= boundaryMarkers_.size())
            DUNE_THROW(Dune::RangeError, "Boundary segment index "<< boundarySegmentIndex << " bigger than number of boundary segments in grid.\n"
                                         "Make sure to call this function only for boundaries that were defined as physical entities in gmsh.");
        return boundaryMarkers_[boundarySegmentIndex];
    }

    /*!
     * \brief Return the boundary domain marker (Gmsh physical entity number) of an intersection
              Only available when using Gmsh with GridParameterGroup.DomainMarkers = 1.
     * \param intersection The intersection to be evaluated
     */
    int getBoundaryDomainMarker(const Intersection& intersection) const
    { return getBoundaryDomainMarker(intersection.boundarySegmentIndex()); }

    /*!
     * \brief Returns true if an intersection was inserted during grid creation
     */
    bool wasInserted(const Intersection& intersection) const
    { return gridFactory_->wasInserted(intersection); }

    /*!
     * \brief Return the element domain marker (Gmsh physical entity number) of an element.
              Only available when using Gmsh with GridParameterGroup.DomainMarkers = 1.
     * \param element The element to be evaluated
     */
    int getElementDomainMarker(const Element& element) const
    {
        if (dataSourceType_ != DataSourceType::gmsh)
            DUNE_THROW(Dune::InvalidStateException, "Domain markers are only available for gmsh grids.");

        // parameters are only given for level 0 elements
        auto level0element = element;
        while (level0element.hasFather())
            level0element = level0element.father();

        // in the parallel case the data is load balanced and then accessed with indices of the index set
        // for UGGrid element data is read on all processes since UGGrid can't communicate element data (yet)
        if (gridPtr_->comm().size() > 1 && !Detail::isUG<Grid>::value)
            return elementMarkers_[gridPtr_->levelGridView(0).indexSet().index(level0element)];
        else
            return elementMarkers_[gridFactory_->insertionIndex(level0element)];
    }

    /*!
     * \brief Create a data handle for communication of the data in parallel simulations
     * \note this data hande is the default
     */
    template<bool ug = Detail::isUG<Grid>::value, typename std::enable_if_t<!ug, int> = 0>
    DataHandle createGmshDataHandle()
    {
        return DataHandle(*gridPtr_, *gridFactory_, elementMarkers_, boundaryMarkers_, faceMarkers_);
    }

    /*!
     * \brief Create a data handle for communication of the data in parallel simulations
     * \note this data hande is the specialized for UGGrid since UGGrid can't communicate element data (yet)
     */
    template<bool ug = Detail::isUG<Grid>::value, typename std::enable_if_t<ug, int> = 0>
    DataHandle createGmshDataHandle()
    {
        return DataHandle(*gridPtr_, *gridFactory_, elementMarkers_, boundaryMarkers_);
    }

    // \}

    /*!
     * \name VTK interface functions
     */
    // \{

    /*!
     * \brief Get a element parameter
     * \param element the element
     * \param fieldName the name of the field to read from the vtk data
     */
    double getParameter(const Element& element, const std::string& fieldName) const
    {
        if (dataSourceType_ != DataSourceType::vtk)
            DUNE_THROW(Dune::InvalidStateException, "This access function is only available for data from VTK files.");

        if (cellData_.count(fieldName) == 0)
            DUNE_THROW(Dune::IOError, "No field with the name " << fieldName << " found in cell data");

        // parameters are only given for level 0 elements
        auto level0element = element;
        while (level0element.hasFather())
            level0element = level0element.father();

        return cellData_.at(fieldName)[gridFactory_->insertionIndex(level0element)];
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for vertex data
     * \param vertex the vertex
     * \param fieldName the name of the field to read from the vtk data
     * \note You can only pass vertices that exist on level 0!
     */
    double getParameter(const Vertex& vertex, const std::string& fieldName) const
    {
        if (dataSourceType_ != DataSourceType::vtk)
            DUNE_THROW(Dune::InvalidStateException, "This access function is only available for data from VTK files.");

        if (vertex.level() != 0)
            DUNE_THROW(Dune::IOError, "You can only obtain parameters for level 0 vertices!");

        if (pointData_.count(fieldName) == 0)
            DUNE_THROW(Dune::IOError, "No field with the name " << fieldName << " found in point data");

        return pointData_.at(fieldName)[gridFactory_->insertionIndex(vertex)];
    }

    // \}

private:
    // grid and grid factor for gmsh grid data / vtk grid data
    std::shared_ptr<Grid> gridPtr_;
    std::shared_ptr<Dune::GridFactory<Grid>> gridFactory_;

    /*!
     * \brief Element and boundary domain markers obtained from Gmsh physical entities
     *        They map from element indices / boundary ids to the physical entity number
     */
    std::vector<int> elementMarkers_;
    std::vector<int> boundaryMarkers_;
    std::vector<int> faceMarkers_;

    /*!
     * \brief Cell and vertex data obtained from VTK files
     */
    VTKReader::Data cellData_, pointData_;

    // dgf grid data
    Dune::GridPtr<Grid> dgfGrid_;

    // specify which type of data we have
    // TODO unfortunately all grid readers provide different data types, should be streamlined (changes in Dune)
    DataSourceType dataSourceType_;
};

} // namespace Dumux

#endif
