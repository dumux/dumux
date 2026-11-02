// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief A grid geometry with a prescribed periodic degree of freedom map for assembler tests
 */
#ifndef DUMUX_TEST_ASSEMBLY_PERIODIC_DOF_MAP_HH
#define DUMUX_TEST_ASSEMBLY_PERIODIC_DOF_MAP_HH

#include <cstddef>
#include <ranges>
#include <unordered_map>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/rangegenerators.hh>

namespace Dumux::Test {

/*!
 * \ingroup Assembly
 * \brief A grid geometry that identifies prescribed sets of degrees of freedom as periodic
 * \note The identification only enters the queries of the assemblers through the grid geometry type
 *       (Jacobian pattern and periodic constraints). The discretization itself is not periodic, as the
 *       local views refer to the base grid geometry.
 */
template<class GridGeometry>
class GridGeometryWithPeriodicDofMap : public GridGeometry
{
    using PeriodicDofMap = std::unordered_map<std::size_t, std::vector<std::size_t>>;

public:
    using GridGeometry::GridGeometry;

    //! Identify the given degrees of freedom with each other
    void addPeriodicDofs(const std::vector<std::size_t>& dofs)
    {
        for (const auto dofI : dofs)
            for (const auto dofJ : dofs)
                if (dofI != dofJ)
                    periodicDofMap_[dofI].push_back(dofJ);
    }

    bool isPeriodic() const
    { return !periodicDofMap_.empty(); }

    bool dofOnPeriodicBoundary(std::size_t dofIdx) const
    { return periodicDofMap_.count(dofIdx); }

    const std::ranges::range auto periodicallyMappedDofs(std::size_t dofIdx) const
    { return std::views::all(periodicDofMap_.at(dofIdx)); }

    const PeriodicDofMap& periodicDofMap() const
    { return periodicDofMap_; }

private:
    PeriodicDofMap periodicDofMap_;
};

/*!
 * \ingroup Assembly
 * \brief Identify two pairs and a triple of degrees of freedom away from the boundary
 * \return The number of constraint rows of the identification
 */
template<class GridGeometry>
std::size_t addInteriorPeriodicDofs(GridGeometryWithPeriodicDofMap<GridGeometry>& gridGeometry)
{
    std::vector<bool> onBoundary(gridGeometry.numDofs(), false);
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
            if (scv.boundary())
                onBoundary[scv.dofIndex()] = true;
    }

    std::vector<std::size_t> interiorDofs;
    for (std::size_t dofIdx = 0; dofIdx < onBoundary.size(); ++dofIdx)
        if (!onBoundary[dofIdx])
            interiorDofs.push_back(dofIdx);

    const auto n = interiorDofs.size();
    if (n < 9)
        DUNE_THROW(Dune::InvalidStateException, "Not enough interior degrees of freedom");

    gridGeometry.addPeriodicDofs({interiorDofs[0], interiorDofs[n/2]});
    gridGeometry.addPeriodicDofs({interiorDofs[1], interiorDofs[n/2 + 1]});
    gridGeometry.addPeriodicDofs({interiorDofs[2], interiorDofs[n/3], interiorDofs[2*n/3]});
    return 4;
}

} // end namespace Dumux::Test

#endif
