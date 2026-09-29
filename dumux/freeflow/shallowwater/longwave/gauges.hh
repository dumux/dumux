// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 * \brief Discharge at interior cross-sections of a channel network
 *
 * A gauging station along a channel is not on the boundary, so its discharge is the flux
 * across an interior face. On a control-volume finite element scheme, a one-dimensional
 * element has exactly one interior sub-control-volume face, at its midpoint, so a station is
 * an element and it measures the interior flux of the model through that face.
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_GAUGES_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_GAUGES_HH

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/discretization/localview.hh>

#include "discharge.hh"

namespace Dumux::LongWave {

/*!
 * \ingroup LongWaveModel
 * \brief The stations of a gauging network, located once on the channel geometry
 *
 * Locating the stations only depends on the mesh: the face a station sits on, and which
 * direction along it is downstream, are fixed, while the discharge through that face changes
 * every time step.
 */
template<class GridGeometry>
class StreamGauges
{
    using GridView = typename GridGeometry::GridView;
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    struct Station
    {
        std::string name;
        std::size_t element;
        std::size_t scvf;
        //! +1 if the inside sub-control volume of the face is the upstream one, -1 otherwise
        int sign;
        //! distance between the face and the surveyed position
        Scalar offset;
    };

    /*!
     * \brief Snap each surveyed position onto the nearest interior face of the network
     *
     * \param gridGeometry the grid geometry of the channel network
     * \param names the station names
     * \param positions the surveyed station positions
     * \param elevations the bed elevation of every degree of freedom, which determines the
     *        downstream direction and hence the sign of the reported discharge
     * \param snapRadius the maximum distance between a station and the channel network
     */
    StreamGauges(const GridGeometry& gridGeometry,
                 const std::vector<std::string>& names,
                 const std::vector<GlobalPosition>& positions,
                 const std::vector<Scalar>& elevations,
                 const Scalar snapRadius)
    {
        if (names.size() != positions.size())
            DUNE_THROW(Dune::InvalidStateException,
                       "The gauging network names " << names.size() << " stations but gives "
                       << positions.size() << " positions");
        if (elevations.size() != gridGeometry.numDofs())
            DUNE_THROW(Dune::InvalidStateException,
                       "The channel has " << gridGeometry.numDofs() << " dofs but "
                       << elevations.size() << " bed elevations");

        const auto& mapper = gridGeometry.elementMapper();
        auto fvGeometry = localView(gridGeometry);
        for (std::size_t k = 0; k < names.size(); ++k)
        {
            Station closest{names[k], 0, 0, 1, std::numeric_limits<Scalar>::max()};
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                fvGeometry.bind(element);
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    if (scvf.boundary())
                        continue;
                    const auto distance = (scvf.ipGlobal() - positions[k]).two_norm();
                    if (distance >= closest.offset)
                        continue;

                    const auto& inside = fvGeometry.scv(scvf.insideScvIdx());
                    const auto& outside = fvGeometry.scv(scvf.outsideScvIdx());
                    const auto fall = elevations[inside.dofIndex()] - elevations[outside.dofIndex()];
                    closest = Station{names[k], mapper.index(element), scvf.index(),
                                      fall >= 0.0 ? 1 : -1, distance};
                }
            }
            if (closest.offset > snapRadius)
                DUNE_THROW(Dune::InvalidStateException,
                           "Gauge " << names[k] << " at " << positions[k] << " is "
                           << closest.offset << " m away from the nearest interior face of the "
                           "channel network, which exceeds the snap radius of " << snapRadius
                           << ". The station may be on a tributary that is not part of the "
                           "network, or its position may be given in a different coordinate system.");
            stations_.push_back(closest);
        }
    }

    const std::vector<Station>& stations() const
    { return stations_; }

    std::size_t size() const
    { return stations_.size(); }

    bool empty() const
    { return stations_.empty(); }

    //! Column headers in the order returned by `discharges`
    std::vector<std::string> columns() const
    {
        std::vector<std::string> names;
        for (const auto& station : stations_)
            names.push_back(station.name + "[m^3/s]");
        return names;
    }

    /*!
     * \brief The discharge at every station, positive in downstream direction
     *
     * The sign is independent of the orientation of the mesh, which determines which side
     * of a face is the inside sub-control volume.
     */
    template<class Problem, class GridVariables, class SolutionVector>
    std::vector<double> discharges(const Problem& problem,
                                   const GridVariables& gridVariables,
                                   const SolutionVector& sol) const
    {
        std::vector<double> values(stations_.size(), 0.0);
        if (stations_.empty())
            return values;

        const auto& gridGeometry = problem.gridGeometry();
        const auto& mapper = gridGeometry.elementMapper();
        auto fvGeometry = localView(gridGeometry);
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto index = mapper.index(element);
            bool wanted = false;
            for (const auto& station : stations_)
                wanted = wanted || station.element == index;
            if (!wanted)
                continue;

            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                    continue;
                for (std::size_t k = 0; k < stations_.size(); ++k)
                    if (stations_[k].element == index && stations_[k].scvf == scvf.index())
                        values[k] = stations_[k].sign*discharge(problem, element, fvGeometry,
                                                                elemVolVars, scvf,
                                                                elemFluxVarsCache);
            }
        }
        return values;
    }

private:
    std::vector<Station> stations_;
};

} // end namespace Dumux::LongWave

#endif
