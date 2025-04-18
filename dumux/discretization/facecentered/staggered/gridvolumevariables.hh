// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief The grid volume variables class for face-centered staggered models
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_GRID_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_GRID_VOLUMEVARIABLES_HH

#include <vector>
#include <type_traits>

#include <dumux/parallel/parallel_for.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/facecentered/staggered/elementsolution.hh>
#include <dumux/discretization/facecentered/staggered/elementvolumevariables.hh>

namespace Dumux {

template<class P, class VV>
struct FaceCenteredStaggeredDefaultGridVolumeVariablesTraits
{
    using Problem = P;
    using VolumeVariables = VV;

    template<class GridVolumeVariables, bool cachingEnabled>
    using LocalView = FaceCenteredStaggeredElementVolumeVariables<GridVolumeVariables, cachingEnabled>;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Base class for the grid volume variables
 * \note This class has a cached version and a non-cached version
 * \tparam Traits the traits class injecting the problem, volVar and elemVolVars type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class Traits, bool cachingEnabled = false>
class FaceCenteredStaggeredGridVolumeVariables;

//! specialization in case of storing the volume variables
template<class Traits>
class FaceCenteredStaggeredGridVolumeVariables<Traits, /*cachingEnabled*/true>
{
    using ThisType = FaceCenteredStaggeredGridVolumeVariables<Traits, true>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the volume variables type
    using VolumeVariables = typename Traits::VolumeVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    FaceCenteredStaggeredGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        volumeVariables_.resize(gridGeometry.numScv());
        Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
        {
            const auto element = gridGeometry.element(eIdx);
            const auto fvGeometry = localView(gridGeometry).bindElement(element);
            const auto elemSol = elementSolution(element, sol, gridGeometry);

            for (const auto& scv : scvs(fvGeometry))
                volumeVariables_[scv.index()].update(elemSol, problem, element, scv);
        });
    }

    const VolumeVariables& volVars(const std::size_t scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& volVars(const std::size_t scvIdx)
    { return volumeVariables_[scvIdx]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& volVars(const SubControlVolume scv) const
    { return volumeVariables_[scv.index()]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& volVars(const SubControlVolume scv)
    { return volumeVariables_[scv.index()]; }

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
class FaceCenteredStaggeredGridVolumeVariables<Traits, /*cachingEnabled*/false>
{
    using ThisType = FaceCenteredStaggeredGridVolumeVariables<Traits, false>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the volume variables type
    using VolumeVariables = typename Traits::VolumeVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    FaceCenteredStaggeredGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

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
