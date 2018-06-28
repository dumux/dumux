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
 * \brief Class for grid data attached to dgf or gmsh grid files
 */
#ifndef DUMUX_IO_GRID_DATA_HH
#define DUMUX_IO_GRID_DATA_HH

#include <vector>
#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/gridfactory.hh>

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
    using DataHandle = GmshGridDataHandle<Grid, Dune::GridFactory<Grid>, std::vector<int>>;

public:
    //! constructor for gmsh grid data
    GridData(std::shared_ptr<Grid> grid, std::shared_ptr<Dune::GridFactory<Grid>> factory,
             std::vector<int>&& elementMarkers, std::vector<int>&& boundaryMarkers, std::vector<int>&& faceMarkers = std::vector<int>{})
    : gmshGrid_(grid)
    , gridFactory_(factory)
    , elementMarkers_(elementMarkers)
    , boundaryMarkers_(boundaryMarkers)
    , faceMarkers_(faceMarkers)
    {}

    //! constructor for dgf grid data
    GridData(Dune::GridPtr<Grid> grid)
    : dgfGrid_(grid)
    , isDgfData_(true)
    {}


    /*!
     * \brief Call the parameters function of the DGF grid pointer if available
     */
    template <class Entity>
    const std::vector<double>& parameters(const Entity& entity) const
    {
        if (isDgfData_)
        {
            if (entity.hasFather())
            {
                auto level0entity = entity;
                while(level0entity.hasFather())
                    level0entity = level0entity.father();


                return dgfGrid_.parameters(level0entity);
            }
            else
            {
                return dgfGrid_.parameters(entity);
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
        if (isDgfData_)
            return dgfGrid_.parameters(intersection);
        else
            DUNE_THROW(Dune::InvalidStateException, "The parameters method is only available if the grid was constructed with a DGF file.");
    }

    /*!
     * \brief Return the boundary domain marker (Gmsh physical entity number) of an intersection
              Only available when using Gmsh with GridParameterGroup.DomainMarkers = 1.
     * \param boundarySegmentIndex The boundary segment index of the intersection (intersection.boundarySegmentIndex()
     */
    int getBoundaryDomainMarker(int boundarySegmentIndex) const
    {
        if (!gmshGrid_)
            DUNE_THROW(Dune::InvalidStateException, "Domain markers are only available for gmsh grids.");
        if (boundarySegmentIndex >= boundaryMarkers_.size())
            DUNE_THROW(Dune::RangeError, "Boundary segment index "<< boundarySegmentIndex << " bigger than number of boundary segments in grid.");
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
     * \brief Return the element domain marker (Gmsh physical entity number) of an element.
              Only available when using Gmsh with GridParameterGroup.DomainMarkers = 1.
     * \param element The element to be evaluated
     */
    int getElementDomainMarker(const Element& element) const
    {
        if (!gmshGrid_)
            DUNE_THROW(Dune::InvalidStateException, "Domain markers are only available for gmsh grids.");

        // parameters are only given for level 0 elements
        auto level0element = element;
        while (level0element.hasFather())
            level0element = level0element.father();

        // in the parallel case the data is load balanced and then accessed with indices of the index set
        // for UGGrid element data is read on all processes since UGGrid can't communicate element data (yet)
        if (gmshGrid_->comm().size() > 1 && !Detail::isUG<Grid>::value)
            return elementMarkers_[gmshGrid_->levelGridView(0).indexSet().index(level0element)];
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
        return DataHandle(*gmshGrid_, *gridFactory_, elementMarkers_, boundaryMarkers_, faceMarkers_);
    }

    /*!
     * \brief Create a data handle for communication of the data in parallel simulations
     * \note this data hande is the specialized for UGGrid since UGGrid can't communicate element data (yet)
     */
    template<bool ug = Detail::isUG<Grid>::value, typename std::enable_if_t<ug, int> = 0>
    DataHandle createGmshDataHandle()
    {
        return DataHandle(*gmshGrid_, *gridFactory_, elementMarkers_);
    }

private:
    // grid and grid factor for gmsh grid data
    std::shared_ptr<Grid> gmshGrid_;
    std::shared_ptr<Dune::GridFactory<Grid>> gridFactory_;

    /*!
    * \brief Element and boundary domain markers obtained from Gmsh physical entities
    *        They map from element indices / boundary ids to the physical entity number
    */
    std::vector<int> elementMarkers_;
    std::vector<int> boundaryMarkers_;
    std::vector<int> faceMarkers_;

    // dgf grid data
    Dune::GridPtr<Grid> dgfGrid_;
    bool isDgfData_ = false;
};

} // namespace Dumux

#endif
