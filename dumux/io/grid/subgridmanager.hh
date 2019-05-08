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
 * \brief  A grid manager for dune-subgrid.
 */
#ifndef DUMUX_SUBGRID_MANAGER_HH
#define DUMUX_SUBGRID_MANAGER_HH

#if HAVE_DUNE_SUBGRID

#include <memory>

#include <dune/subgrid/subgrid.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/dgfparser/dgfwriter.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/boundaryflag.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A grid manager for dune-subgrid.
 */
template <class HostGrid>
class SubgridManager
{
    static constexpr auto dim = HostGrid::dimension;

public:
    using Grid = Dune::SubGrid<dim, HostGrid>;

    /*!
     * \brief Make the subgrid.
     */
    template<class ElementSelector>
    static std::unique_ptr<Grid> makeGrid(HostGrid& hostgrid,
                                          const ElementSelector& selector,
                                          const std::string& modelParamGroup = "")
    {
        // A unique pointer to the subgrid.
        auto subgridPtr = std::make_unique<Grid>(hostgrid);

        // A container to store the host grid elements' ids.
        std::set<typename HostGrid::Traits::GlobalIdSet::IdType> elementsForSubgrid;
        const auto& globalIDset = subgridPtr->getHostGrid().globalIdSet();

        // Construct the subgrid.
        subgridPtr->createBegin();

        // Loop over all elements of the host grid and use the selector to
        // choose which elements to add to the subgrid.
        auto hostGridView = subgridPtr->getHostGrid().leafGridView();
        for (const auto& e : elements(hostGridView))
            if(selector(e))
                elementsForSubgrid.insert(globalIDset.template id<0>(e));

        subgridPtr->insertSetPartial(elementsForSubgrid);
        subgridPtr->createEnd();

        // If desired, write out the final subgrid as a dgf file.
        if(getParamFromGroup<bool>(modelParamGroup, "Grid.WriteSubGridToDGF", false))
        {
            const auto postfix = getParamFromGroup<std::string>(modelParamGroup, "Problem.Name", "");
            const std::string name = postfix == "" ? "subgrid" : "subgrid_" + postfix;
            Dune::DGFWriter<typename Grid::LeafGridView> writer(subgridPtr->leafGridView());
            writer.write(name + ".dgf");
        }

        // If desired, write out the hostgrid as vtk file.
        if(getParamFromGroup<bool>(modelParamGroup, "Grid.WriteSubGridToVtk", false))
        {
            const auto postfix = getParamFromGroup<std::string>(modelParamGroup, "Problem.Name", "");
            const std::string name = postfix == "" ? "subgrid" : "subgrid_" + postfix;
            Dune::VTKWriter<typename Grid::LeafGridView> vtkWriter(subgridPtr->leafGridView());
            vtkWriter.write(name);
        }

        // Return a unique pointer to the subgrid.
        return subgridPtr;
    }
};

//! dune-subgrid doesn't have this implemented
template<int dim, class HostGrid>
class BoundaryFlag<Dune::SubGrid<dim, HostGrid>>
{
public:
    BoundaryFlag() : flag_(-1) {}

    template<class Intersection>
    BoundaryFlag(const Intersection& i) : flag_(-1) {}

    using value_type = int;

    value_type get() const
    { DUNE_THROW(Dune::NotImplemented, "Sub-grid doesn't implement boundary segment indices!"); }

private:
    int flag_;
};

} // end namespace Dumux

#endif // HAVE_DUNE_SUBGRID
#endif // DUMUX_SUBGRID_MANAGER_HH
