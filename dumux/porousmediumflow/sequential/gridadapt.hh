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
 * \brief Base class for h-adaptive sequential models.
 */
#ifndef DUMUX_GRIDADAPT_HH
#define DUMUX_GRIDADAPT_HH

#include <unordered_map>
#include <dune/grid/common/partitionset.hh>

#include "properties.hh"
#include "gridadaptproperties.hh"

namespace Dumux
{

/*!\ingroup IMPET
 * @brief Standard Module for h-adaptive simulations
 *
 * This class is created by the problem class with the template
 * parameters <TypeTag, true> and provides basic functionality
 * for adaptive methods:
 *
 * A standard implementation adaptGrid() will prepare everything
 * to calculate the next pressure field on the new grid.
 */
template<class TypeTag, bool adaptive>
class GridAdapt
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    using Grid = typename GridView::Grid;
    using LeafGridView = typename Grid::LeafGridView;
    using Element = typename Grid::template Codim<0>::Entity;

    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using AdaptionIndicator = GetPropType<TypeTag, Properties::AdaptionIndicator>;
    using AdaptionInitializationIndicator = GetPropType<TypeTag, Properties::AdaptionInitializationIndicator>;

    using IdSet = typename Grid::Traits::LocalIdSet;
    using IdType = typename IdSet::IdType;

public:
    /*!
     * Constructor for h-adaptive simulations (adaptive grids)
     * @param problem The problem
     */
    GridAdapt (Problem& problem)
        : problem_(problem), adaptionIndicator_(problem), marked_(0), coarsened_(0)
    {
        levelMin_ = getParam<int>("GridAdapt.MinLevel");
        levelMax_ = getParam<int>("GridAdapt.MaxLevel");
        adaptationInterval_ = getParam<int>("GridAdapt.AdaptionInterval", 1);

        if (levelMin_ < 0)
            Dune::dgrave <<  __FILE__<< ":" <<__LINE__
                         << " :  Dune cannot coarsen to gridlevels smaller 0! "<< std::endl;
    }

    /*!
     * @brief Initalization method of the h-adaptive module
     *
     * Prepares the grid for simulation after the initialization of the
     * problem. The applied indicator is selectable via the property
     * AdaptionInitializationIndicator
     */
    void init()
    {
        adaptionIndicator_.init();

        if (!getParam<bool>("GridAdapt.EnableInitializationIndicator"))
            return;

        AdaptionInitializationIndicator adaptionInitIndicator(problem_, adaptionIndicator_);

        int maxIter = 2*levelMax_;
        int iter = 0;
        while (iter <= maxIter)
        {
            adaptGrid(adaptionInitIndicator);

            if (!wasAdapted())
            {
                break;
            }

            int shouldInitialize = adaptionInitIndicator.initializeModel();
            if (problem_.grid().comm().max(shouldInitialize))
                problem_.model().initialize();

            iter++;
        }
    }

    /*!
     * @brief Standard method to adapt the grid
     *
     * This method is called from IMPETProblem::preTimeStep() if
     * adaptive grids are used in the simulation. It uses the standard
     * indicator (selected by the property AdaptionIndicator) and forwards to
     * with it to the ultimate method adaptGrid(indicator), which
     * uses a standard procedure for adaptivity:
     * 1) Determine the refinement indicator
     * 2) Mark the elements
     * 3) Store primary variables in a map
     * 4) Adapt the grid, adapt variables sizes, update mappers
     * 5) Reconstruct primary variables, regain secondary variables
     */
    void adaptGrid()
    {
        adaptGrid(adaptionIndicator_) ;
    }

    /*!
     * @brief Method to adapt the grid with individual indicator vector
     *
     * @param indicator The refinement indicator that is applied
     *
     * This method is called by an user-defined preTimeStep() of
     * the applied problem and takes a given vector with indicator
     * values.
     *
     * It uses a standard procedure for adaptivity:
     * 1) Determine the refinement indicator
     * 2) Mark the elements
     * 3) Store primary variables in a map
     * 4) Adapt the grid, adapt variables sizes, update mappers
     * 5) Reconstruct primary variables, regain secondary variables
     */
    template<class Indicator>
    void adaptGrid(Indicator& indicator)
    {
        // reset internal counter for marked elements
        marked_ = coarsened_ = 0;

        // check for adaption interval: Adapt only at certain time step indices
        if (problem_.timeManager().timeStepIndex() % adaptationInterval_ != 0)
            return;

        /**** 1) determine refining parameter if standard is used ***/
        // if not, the indicatorVector and refinement Bounds have to
        // specified by the problem through setIndicator()
        indicator.calculateIndicator();

        /**** 2) mark elements according to indicator     *********/
        markElements(indicator);

        // abort if nothing in grid is marked
        int sumMarked = problem_.grid().comm().sum(marked_);
        int sumCoarsened = problem_.grid().comm().sum(coarsened_);
        if (sumMarked == 0 && sumCoarsened == 0)
            return;
        else
            Dune::dinfo << marked_ << " cells have been marked_ to be refined, "
                        << coarsened_ << " to be coarsened." << std::endl;

        /****  2b) Do pre-adaption step    *****/
        problem_.grid().preAdapt();
        problem_.preAdapt();

        /****  3) Put primary variables in a map         *********/
        problem_.variables().storePrimVars(problem_);

        /****  4) Adapt Grid and size of variable vectors    *****/
        problem_.grid().adapt();

        //        forceRefineRatio(1);

        // update mapper to new cell indices
        problem_.variables().elementMapper().update();

        // adapt size of vectors
        problem_.variables().adaptVariableSize(problem_.variables().elementMapper().size());

        /****  5) (Re-)construct primary variables to new grid **/
        problem_.variables().reconstructPrimVars(problem_);

        // delete markers in grid
        problem_.grid().postAdapt();

        // call method in problem for potential output etc.
        problem_.postAdapt();

        return;
    }

    /*!
     * Mark Elements for grid refinement according to applied Indicator
     * @return Total ammount of marked cells
     */
    template<class Indicator>
    void markElements(Indicator& indicator)
    {
        using CoarsenMarkerType = std::unordered_map<IdType, IdType>;
        CoarsenMarkerType coarsenMarker;
        const IdSet& idSet(problem_.grid().localIdSet());

        for (const auto& element : elements(problem_.gridView()))
        {
            // only mark non-ghost elements
            if (element.partitionType() == Dune::GhostEntity)
                continue;

            // refine?
            if (indicator.refine(element) && element.level() < levelMax_)
            {
                problem_.grid().mark( 1,  element);
                ++marked_;

                // this also refines the neighbor elements
                checkNeighborsRefine_(element);
            }
            if (indicator.coarsen(element) && element.hasFather())
            {
                IdType idx = idSet.id(element.father());
                auto it = coarsenMarker.find(idx);
                if (it != coarsenMarker.end())
                {
                    ++it->second;
                }
                else
                {
                    coarsenMarker[idx] = 1;
                }
            }
        }
        // coarsen
        for (const auto& element : elements(problem_.gridView()))
        {
            // only mark non-ghost elements
            if (element.partitionType() == Dune::GhostEntity)
                continue;

            if (indicator.coarsen(element) && element.level() > levelMin_)
            {
                IdType idx = idSet.id(element.father());
                auto it = coarsenMarker.find(idx);
                if (it != coarsenMarker.end())
                {
                    if (problem_.grid().getMark(element) == 0
                        && it->second == element.geometry().corners())
                    {
                        // check if coarsening is possible
                        bool coarsenPossible = true;
                        for(const auto& intersection : intersections(problem_.gridView(), element))
                        {
                            if(intersection.neighbor())
                            {
                                auto outside = intersection.outside();
                                if ((problem_.grid().getMark(outside) > 0)
                                    || outside.level() > element.level())
                                {
                                    coarsenPossible = false;
                                }
                            }
                        }

                        if(coarsenPossible)
                        {
                            problem_.grid().mark( -1, element );
                            ++coarsened_;
                        }
                    }
                }
            }
        }
    }

    /*!
     * @brief Returns true if grid cells have been marked for adaptation
     */
    bool wasAdapted()
    {
        int sumMarked = problem_.grid().comm().sum(marked_);
        int sumCoarsened = problem_.grid().comm().sum(coarsened_);

        return (sumMarked != 0 || sumCoarsened != 0);
    }

    /*!
     * Sets minimum and maximum refinement levels
     *
     * @param levMin minimum level for coarsening
     * @param levMax maximum level for refinement
     */
    void setLevels(int levMin, int levMax)
    {
        if (levMin < 0)
            Dune::dgrave <<  __FILE__<< ":" <<__LINE__
                         << " :  Dune cannot coarsen to gridlevels smaller 0! "<< std::endl;
        levelMin_ = levMin;
        levelMax_ = levMax;
    }

    /*!
     * @brief Returns maximum refinement level
     *
     * The value is the assign maximum possible level,
     * not the actual maximum level of the grid.
     * @return levelMax_ maximum level for refinement
     */
    int getMaxLevel() const
    {
        return levelMax_;
    }
    /*!
     * @brief Returns minimum refinement level
     *
     * The value is the assign minimum possible level,
     * not the actual minimum level of the grid.
     * @return levelMin_ minimum level for coarsening
     */
    int getMinLevel() const
    {
        return levelMin_;
    }

    AdaptionIndicator& adaptionIndicator()
    {
        return adaptionIndicator_;
    }

    AdaptionIndicator& adaptionIndicator() const
    {
        return adaptionIndicator_;
    }

private:
    /*!
     * @brief Method ensuring the refinement ratio of 2:1
     *
     * For any given entity, a loop over the neighbors checks weather the
     * entities refinement would require that any of the neighbors has
     * to be refined, too.
     * This is done recursively over all levels of the grid.
     *
     * @param entity Element of interest that is to be refined
     * @param level level of the refined entity: it is at least 1
     * @return true if everything was successful
     */
    bool checkNeighborsRefine_(const Element &entity, int level = 1)
    {
        // this also refines the neighbor elements
        for(const auto& intersection : intersections(problem_.gridView(), entity))
        {
            if(!intersection.neighbor())
                continue;

            auto outside = intersection.outside();

            // only mark non-ghost elements
            if (outside.partitionType() == Dune::GhostEntity)
                continue;

            if ((outside.level() < levelMax_)
                && (outside.level() < entity.level()))
            {
                problem_.grid().mark(1, outside);
                ++marked_;

                if(level != levelMax_)
                    checkNeighborsRefine_(outside, ++level);
            }
        }
        return true;
    }


    /*!
     * \brief Enforces a given refine ratio after grid was adapted
     *
     * If the refine ratio is not taken into consideration during
     * marking, then this method ensures a certain ratio.
     *
     * @param maxLevelDelta The maximum level difference (refine ratio)
     *             between neighbors.
     */
    void forceRefineRatio(int maxLevelDelta = 1)
    {
        LeafGridView leafGridView = problem_.gridView();
        // delete all existing marks
        problem_.grid().postAdapt();
        bool done;
        do
        {
            // run through all cells
            done=true;
            for (const auto& element : elements(problem_.gridView()))
            {
                // only mark non-ghost elements
                if (element.partitionType() == Dune::GhostEntity)
                    continue;

                // run through all neighbor-cells (intersections)
                for (const auto& intersection : intersections(leafGridView, element))
                {
                    if(!intersection.neighbor())
                        continue;

                    if (element.level() + maxLevelDelta < intersection.outside().level())
                    {
                        problem_.grid().mark( 1, element );
                        done=false;
                    }
                }
            }
            if (done==false)
            {
                // adapt the grid
                problem_.grid().adapt();
                // delete marks
                problem_.grid().postAdapt();
            }
        }
        while (done!=true);
    }

    // private Variables
    Problem& problem_;
    AdaptionIndicator adaptionIndicator_;

    int marked_;
    int coarsened_;

    int levelMin_; //!< minimum allowed level
    int levelMax_; //!< maximum allowed level

    int adaptationInterval_; //!< time step interval for adaption
};

/*!
 * @brief Class for NON-adaptive simulations
 *
 * This class provides empty methods for non-adaptive simulations
 * for compilation reasons. If adaptivity is desired, create the
 * class with template arguments <TypeTag, true> instead.
 */
template<class TypeTag>
class GridAdapt<TypeTag, false>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ScalarSolutionType = typename SolutionTypes::ScalarSolution;

public:
    void init()
    {}
    void adaptGrid()
    {}
    bool wasAdapted()
    {
        return false;
    }
    void setLevels(int, int)
    {}
    void setTolerance(int, int)
    {}
    void setIndicator(const ScalarSolutionType&,
                            const Scalar&, const Scalar&)
    {}
    GridAdapt (Problem& problem)
    {}
};

}
#endif
