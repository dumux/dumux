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
 * \ingroup CCDiscretization
 * \brief The grid volume variables class for cell centered models
 */
#ifndef DUMUX_DISCRETIZATION_CC_GRID_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CC_GRID_VOLUMEVARIABLES_HH

#include <vector>
#include <type_traits>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Base class for the grid volume variables
 * \note This class has a cached version and a non-cached version
 * \tparam Traits the traits class injecting the problem, volVar and elemVolVars type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class Traits, bool cachingEnabled = false>
class CCGridVolumeVariables {};

//! specialization in case of storing the volume variables
template<class Traits>
class CCGridVolumeVariables<Traits, /*cachingEnabled*/true>
{
    using ThisType = CCGridVolumeVariables<Traits, true>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the volume variables type
    using VolumeVariables = typename Traits::VolumeVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CCGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        const auto numScv = gridGeometry.numScv();
        volumeVariables_.resize(numScv);

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto elemSol = elementSolution(element, sol, gridGeometry);
                volumeVariables_[scv.dofIndex()].update(elemSol, problem(), element, scv);
            }
        }
    }

    const VolumeVariables& volVars(const std::size_t scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& volVars(const std::size_t scvIdx)
    { return volumeVariables_[scvIdx]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& volVars(const SubControlVolume& scv) const
    { return volumeVariables_[scv.dofIndex()]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& volVars(const SubControlVolume& scv)
    { return volumeVariables_[scv.dofIndex()]; }

    // required for compatibility with the box method
    const VolumeVariables& volVars(const std::size_t scvIdx, const std::size_t localIdx) const
    { return volumeVariables_[scvIdx]; }

    // required for compatibility with the box method
    VolumeVariables& volVars(const std::size_t scvIdx, const std::size_t localIdx)
    { return volumeVariables_[scvIdx]; }

    //! The problem we are solving
    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
    std::vector<VolumeVariables> volumeVariables_;
};


//! Specialization when the current volume variables are not stored globally
template<class Traits>
class CCGridVolumeVariables<Traits, /*cachingEnabled*/false>
{
    using ThisType = CCGridVolumeVariables<Traits, false>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the volume variables type
    using VolumeVariables = typename Traits::VolumeVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CCGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol) {}

    //! The problem we are solving
    const Problem& problem() const
    { return *problemPtr_;}

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
