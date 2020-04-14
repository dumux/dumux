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
 * \ingroup Adaptive
 * \brief Class defining an initialization indicator for grid adaption
 */
#ifndef DUMUX_GRIDADAPTINITIALIZATIONINDICATOR_HH
#define DUMUX_GRIDADAPTINITIALIZATIONINDICATOR_HH

#include <memory>

#include <dune/geometry/type.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Adaptive
 * \brief Class defining an initialization indicator for grid adaption.
 *        Refines at sources and boundaries. Use before computing on the grid.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridAdaptInitializationIndicator
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::Traits::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    static constexpr bool isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box;

public:

    /*!
     * \brief Constructor
     * \note Reads the folloing parameters from the parameter tree
     *       - Adaptive.MinLevel The minimum refinement level
     *       - Adaptive.MaxLevel The maximum refinement level
     *       - Adaptive.RefineAtDirichletBC If to refine at Dirichlet boundaries (default: true)
     *       - Adaptive.RefineAtFluxBC If to refine at Neumann/Robin boundaries (default: true)
     *       - Adaptive.RefineAtSource If to refine where source terms are specified (default: true)
     *       - Adaptive.BCRefinementThreshold The threshold above which fluxes are treated as non-zero (default: 1e-10)
     * \param problem The problem object
     * \param gridGeometry The finite volume geometry of the grid
     * \param gridVariables The secondary variables on the grid
     */
    GridAdaptInitializationIndicator(std::shared_ptr<const Problem> problem,
                                     std::shared_ptr<const GridGeometry> gridGeometry,
                                     std::shared_ptr<const GridVariables> gridVariables)
    : problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , minLevel_(getParamFromGroup<int>(problem->paramGroup(), "Adaptive.MinLevel"))
    , maxLevel_(getParamFromGroup<int>(problem->paramGroup(), "Adaptive.MaxLevel"))
    , refineAtDirichletBC_(getParamFromGroup<bool>(problem->paramGroup(), "Adaptive.RefineAtDirichletBC", true))
    , refineAtFluxBC_(getParamFromGroup<bool>(problem->paramGroup(), "Adaptive.RefineAtFluxBC", true))
    , refineAtSource_(getParamFromGroup<bool>(problem->paramGroup(), "Adaptive.RefineAtSource", true))
    , eps_(getParamFromGroup<Scalar>(problem->paramGroup(), "Adaptive.BCRefinementThreshold", 1e-10))
    {}

    /*!
     * \brief Function to set the minimum allowed level.
     */
    void setMinLevel(std::size_t minLevel)
    {
        minLevel_ = minLevel;
    }

    /*!
     * \brief Function to set the maximum allowed level.
     */
    void setMaxLevel(std::size_t maxLevel)
    {
        maxLevel_ = maxLevel;
    }

    /*!
     * \brief Function to set the minumum/maximum allowed levels.
     */
    void setLevels(std::size_t minLevel, std::size_t maxLevel)
    {
        minLevel_ = minLevel;
        maxLevel_ = maxLevel;
    }

    /*!
     * \brief Function to set if refinement at Dirichlet boundaries is to be used.
     */
    void setRefinementAtDirichletBC(bool refine)
    {
        refineAtDirichletBC_ = refine;
    }

    /*!
     * \brief Function to set if refinement at sources is to be used.
     */
    void setRefinementAtSources(bool refine)
    {
        refineAtSource_ = refine;
    }

    /*!
     * \brief Function to set if refinement at sources is to be used.
     */
    void setRefinementAtFluxBoundaries(bool refine)
    {
        refineAtFluxBC_ = refine;
    }

    /*!
     * \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     * \param sol The solution for which indicator is to be calculated
     */
    template<class SolutionVector>
    void calculate(const SolutionVector& sol)
    {
        //! prepare an indicator for refinement
        indicatorVector_.assign(gridGeometry_->gridView().size(0), false);

        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            const auto eIdx = gridGeometry_->elementMapper().index(element);

            //! refine any element being below the minimum level
            if (element.level() < minLevel_)
            {
                indicatorVector_[eIdx] = true;
                continue; // proceed to the next element
            }

            // If refinement at sources/BCs etc is deactivated, skip the rest
            if (!refineAtSource_ && !refineAtFluxBC_ && !refineAtDirichletBC_)
                continue;

            // if the element is already on the maximum permissive level, skip rest
            if (element.level() == maxLevel_)
                continue;

            // get the fvGeometry and elementVolVars needed for the bc and source interfaces
            auto fvGeometry = localView(*gridGeometry_);
            fvGeometry.bind(element);

            auto elemVolVars = localView(gridVariables_->curGridVolVars());
            elemVolVars.bind(element, fvGeometry, sol);

            // elemFluxVarsCache for neumann interface
            auto elemFluxVarsCache = localView(gridVariables_->gridFluxVarsCache());
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            //! Check if we have to refine around a source term
            if (refineAtSource_)
            {
                for (const auto& scv : scvs(fvGeometry))
                {
                    auto source = problem_->source(element, fvGeometry, elemVolVars, scv);
                    auto pointSource = problem_->scvPointSources(element, fvGeometry, elemVolVars, scv);
                    if (source.infinity_norm() + pointSource.infinity_norm() > eps_)
                    {
                        indicatorVector_[eIdx] = true;
                        break; // element is marked, escape scv loop
                    }
                }
            }

            //! Check if we have to refine at the boundary
            if (!indicatorVector_[eIdx]                       // proceed if element is not already marked
                && element.hasBoundaryIntersections()         // proceed if element is on boundary
                && (refineAtDirichletBC_ || refineAtFluxBC_)) // proceed if boundary refinement is active
            {
                // cell-centered schemes
                if (!isBox)
                {
                    for (const auto& scvf : scvfs(fvGeometry))
                    {
                        // skip non-boundary scvfs
                        if (!scvf.boundary())
                            continue;

                        const auto bcTypes = problem_->boundaryTypes(element, scvf);
                        // We assume pure BCs, mixed boundary types are not allowed anyway!
                        if(bcTypes.hasOnlyDirichlet() && refineAtDirichletBC_)
                        {
                            indicatorVector_[eIdx] = true;
                            break; // element is marked, escape scvf loop
                        }

                        // we are on a pure Neumann boundary
                        else if(refineAtFluxBC_)
                        {
                            const auto fluxes = problem_->neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                            if (fluxes.infinity_norm() > eps_)
                            {
                                indicatorVector_[eIdx] = true;
                                break; // element is marked, escape scvf loop
                            }
                        }
                    }
                }
                // box-scheme
                else
                {
                    // container to store bcTypes
                    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
                    std::vector<BoundaryTypes> bcTypes(fvGeometry.numScv());

                    // Get bcTypes and maybe mark for refinement on Dirichlet boundaries
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        bcTypes[scv.localDofIndex()] = problem_->boundaryTypes(element, scv);
                        if (refineAtDirichletBC_ && bcTypes[scv.localDofIndex()].hasDirichlet())
                        {
                            indicatorVector_[eIdx] = true;
                            break; // element is marked, escape scv loop
                        }
                    }

                    // If element hasn't been marked because of Dirichlet BCS, check Neumann BCs
                    if (!indicatorVector_[eIdx] && refineAtFluxBC_)
                    {
                        for (const auto& scvf : scvfs(fvGeometry))
                        {
                            //! check if scvf is on Neumann boundary
                            if (scvf.boundary() && bcTypes[scvf.insideScvIdx()].hasNeumann())
                            {
                                const auto fluxes = problem_->neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                                if (fluxes.infinity_norm() > eps_)
                                {
                                    indicatorVector_[eIdx] = true;
                                    break; // element is marked, escape scvf loop
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*! \brief function call operator to return mark
     *
     *  \return  1 if an element should be refined
     *          -1 if an element should be coarsened
     *           0 otherwise
     *
     *  \note In this initialization indicator implementation
     *        element coarsening is not considered.
     *
     *  \param element A grid element
     */
    int operator() (const Element& element) const
    {
        if (indicatorVector_[gridGeometry_->elementMapper().index(element)])
            return 1;
        return 0;
    }

private:
    std::shared_ptr<const Problem> problem_;               //!< The problem to be solved
    std::shared_ptr<const GridGeometry> gridGeometry_; //!< The finite volume grid geometry
    std::shared_ptr<const GridVariables> gridVariables_;   //!< The secondary variables on the grid
    std::vector<bool> indicatorVector_;                    //!< Indicator for BCs/sources

    int minLevel_;             //!< The minimum allowed level
    int maxLevel_;             //!< The maximum allowed level
    bool refineAtDirichletBC_; //!< Specifies if it should be refined at Dirichlet BCs
    bool refineAtFluxBC_;      //!< Specifies if it should be refined at non-zero Neumann BCs
    bool refineAtSource_;      //!< Specifies if it should be refined at sources
    Scalar eps_;               //!< Threshold for refinement at sources/BCS
};

} // end namespace Dumux

#endif
