// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief The grid volume variables class for control-volume finite element methods
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_GRID_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CVFE_GRID_VOLUMEVARIABLES_HH

#include <vector>
#include <type_traits>

#include <dumux/parallel/parallel_for.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/elementvolumevariables.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>

namespace Dumux {

template<class P, class VV>
struct CVFEDefaultGridVolumeVariablesTraits
{
    using Problem = P;
    using VolumeVariables = VV;

    template<class GridVolumeVariables, bool cachingEnabled>
    using LocalView = CVFEElementVolumeVariables<GridVolumeVariables, cachingEnabled>;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Base class for the grid volume variables
 */
template<class Traits, bool enableCaching>
class CVFEGridVolumeVariables;

// specialization in case of storing the volume variables
template<class Traits>
class CVFEGridVolumeVariables<Traits, /*cachingEnabled*/true>
{
    using ThisType = CVFEGridVolumeVariables<Traits, true>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the volume variables type
    using VolumeVariables = typename Traits::VolumeVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        volumeVariables_.resize(gridGeometry.gridView().size(0));
        Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
        {
            const auto element = gridGeometry.element(eIdx);
            const auto fvGeometry = localView(gridGeometry).bindElement(element);

            // get the element solution
            auto elemSol = elementSolution(element, sol, gridGeometry);

            // update the volvars of the element
            volumeVariables_[eIdx].resize(fvGeometry.numScv());
            for (const auto& scv : scvs(fvGeometry))
                volumeVariables_[eIdx][scv.indexInElement()].update(elemSol, problem, element, scv);
        });
    }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& volVars(const SubControlVolume& scv) const
    { return volumeVariables_[scv.elementIndex()][scv.indexInElement()]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& volVars(const SubControlVolume& scv)
    { return volumeVariables_[scv.elementIndex()][scv.indexInElement()]; }

    const VolumeVariables& volVars(const std::size_t eIdx, const std::size_t scvIdx) const
    { return volumeVariables_[eIdx][scvIdx]; }

    VolumeVariables& volVars(const std::size_t eIdx, const std::size_t scvIdx)
    { return volumeVariables_[eIdx][scvIdx]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
    std::vector<std::vector<VolumeVariables>> volumeVariables_;
};


// Specialization when the current volume variables are not stored
template<class Traits>
class CVFEGridVolumeVariables<Traits, /*cachingEnabled*/false>
{
    using ThisType = CVFEGridVolumeVariables<Traits, false>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the volume variables type
    using VolumeVariables = typename Traits::VolumeVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_;}

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
