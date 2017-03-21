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
 * \brief Base class for h-adaptive implicit models.
 */
#ifndef DUMUX_IMPLICIT_GRIDADAPT_HH
#define DUMUX_IMPLICIT_GRIDADAPT_HH

#include "gridadaptproperties.hh"
#include <unordered_map>

#include <dune/common/exceptions.hh>

#include <dumux/common/propertysystem.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(SolutionVector);
}

/*!\ingroup ImplicitGridAdapt
 * @brief Standard Module for h-adaptive simulations
 *
 * This class is created by the problem class with the template
 * parameters <TypeTag, true> and provides basic functionality
 * for adaptive methods:
 *
 * A standard implementation adaptGrid() will prepare everything
 * to calculate the next primary variables vector on the new grid.
 */
template<class TypeTag, bool adaptive>
class ImplicitGridAdapt
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem)  Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GridView::Grid Grid;
    typedef typename Grid::LeafGridView LeafGridView;
    typedef typename Grid::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, AdaptionIndicator) AdaptionIndicator;
    typedef typename GET_PROP_TYPE(TypeTag, AdaptionInitializationIndicator) AdaptionInitializationIndicator;
    typedef typename GET_PROP_TYPE(TypeTag, AdaptionHelper) AdaptionHelper;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:
    /*!
     * Constructor for h-adaptive simulations (adaptive grids)
     * @param problem The problem
     */
    ImplicitGridAdapt (Problem& problem)
        : problem_(problem),
          adaptionHelper_(problem),
          adaptionIndicator_(problem),
          sumRefine_(0),
          sumCoarsen_(0)
    {
            levelMin_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MinLevel);
            levelMax_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
            adaptionInterval_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, AdaptionInterval);

            if (levelMin_ < 0)
            {
                DUNE_THROW(Dune::InvalidStateException, "Coarsening the level 0 entities is not possible! Choose MinLevel >= 0");
            }
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

        if (!GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, EnableInitializationIndicator))
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
            {
                problem_.model().init(problem_);
            }

            ++iter;
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
        if(levelMax_ > levelMin_)
            adaptGrid(adaptionIndicator_) ;
    }

    /*!
     * @brief Returns true if grid cells have been marked for adaption
     */
    bool wasAdapted() const
    {
        return (sumRefine_ != 0 || sumCoarsen_ != 0);
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
            DUNE_THROW(Dune::InvalidStateException, "Coarsening the level 0 entities is not possible!");
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
        // check for adaption interval: Adapt only at certain time step indices
        if (problem_.timeManager().timeStepIndex() % adaptionInterval_ != 0)
            return;

        /**** 1) determine refining parameter if standard is used ***/
        // if not, the indicatorVector and refinement Bounds have to
        // specified by the problem through setIndicator()
        indicator.calculateIndicator();

        /**** 2) mark elements according to indicator     *********/
        int refine, coarsen;
        std::tie(refine, coarsen) = markElements(indicator);

        // abort if nothing in grid is marked
        sumRefine_ = problem_.grid().comm().sum(refine);
        sumCoarsen_ = problem_.grid().comm().sum(coarsen);
        if (sumRefine_ == 0 && sumCoarsen_ == 0)
            return;
        else if (problem_.grid().comm().rank() == 0)
            Dune::dinfo << sumRefine_ << " cells have been marked to be refined, "
                        << sumCoarsen_ << " to be coarsened." << std::endl;

        /****  2b) Do pre-adaption step    *****/
        problem_.grid().preAdapt();
        problem_.preAdapt();

        /****  3) Put primary variables in a map         *********/
        adaptionHelper_.storePrimVars(problem_);

        /****  4) Adapt Grid and size of variable vectors    *****/
        problem_.grid().adapt();

        // update mapper to new cell indices
        problem_.elementMapper().update();
        problem_.vertexMapper().update();

        // adapt size of vectors
        problem_.model().adaptVariableSize();

        /****  5) (Re-)construct primary variables to new grid **/
        adaptionHelper_.reconstructPrimVars(problem_);

        // delete markers in grid
        problem_.grid().postAdapt();

        // call method in problem for potential output etc.
        problem_.postAdapt();

        return;
    }

    /*!
     * Mark Elements for grid refinement according to applied Indicator
     * @param indicator The refinement indicator that is applied
     */
    template<class Indicator>
    std::pair<int, int> markElements(Indicator& indicator)
    {
        int refine = 0;
        int coarsen = 0;
        for (const auto& element : elements(problem_.gridView()))
        {
            // only mark non-ghost elements
            if (element.partitionType() != Dune::GhostEntity)
            {
                // refine?
                if (indicator.refine(element) && element.level() < levelMax_)
                {
                    problem_.grid().mark( 1,  element);
                    ++refine;

                    // this also refines the neighbor elements
                    refine += checkNeighborsRefine_(element);
                }
                if (indicator.coarsen(element) && element.hasFather())
                {
                    problem_.grid().mark( -1, element );
                    ++coarsen;
                }
            }
        }
        return std::make_pair(refine, coarsen);
    }

    /*!
     * @brief Method ensuring the refinement ratio of 2:1
     *
     * For any given element, a loop over the neighbors checks weather the
     * entities refinement would require that any of the neighbors has
     * to be refined, too.
     * This is done recursively over all levels of the grid.
     *
     * @param element Element of interest that is to be refined
     * @param level level of the refined element: it is at least 1
     * @return true if everything was successful
     */
    int checkNeighborsRefine_(const Element &element, int level = 1)
    {
        // this also refines the neighbor elements
        int refine = 0;
        for(const auto& intersection : intersections(problem_.gridView(), element))
        {
            if(!intersection.neighbor())
                continue;

            auto outside = intersection.outside();

            // only mark non-ghost elements
            if (outside.partitionType() == Dune::GhostEntity)
                continue;

            if ((outside.level() < levelMax_)
                && (outside.level() < element.level()))
            {
                problem_.grid().mark(1, outside);
                ++refine;

                if(level != levelMax_)
                    refine += checkNeighborsRefine_(outside, ++level);
            }
        }
        return refine;
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
        auto leafGridView = problem_.gridView();
        // delete all existing marks
        problem_.grid().postAdapt();
        bool done;
        do
        {
            // run through all cells
            done=true;
            for (const auto& element : elements(leafGridView))
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
                        problem_.grid().mark(1, element);
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
        while (!done);
    }

    // private Variables
    Problem& problem_;
    AdaptionHelper adaptionHelper_;
    AdaptionIndicator adaptionIndicator_;

    int sumRefine_;
    int sumCoarsen_;

    int levelMin_;
    int levelMax_;

    int adaptionInterval_;
};

/*!
 * @brief Class for non-adaptive simulations
 *
 * This class provides empty methods for non-adaptive simulations
 * for compilation reasons. If adaptivity is desired, create the
 * class with template arguments <TypeTag, true> instead.
 */
template<class TypeTag>
class ImplicitGridAdapt<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem)     Problem;
    typedef typename GET_PROP(TypeTag, SolutionVector) SolutionVector;

public:
    void init()
    {}
    void adaptGrid()
    {}
    bool wasAdapted() const
    { return false; }
    void setLevels(int, int)
    {}
    void setTolerance(int, int)
    {}
    void setIndicator(const SolutionVector&,
                            const Scalar&, const Scalar&)
    {}
    ImplicitGridAdapt (Problem& problem)
    {}
};

}
#endif /* DUMUX_IMPLICIT_GRIDADAPT_HH */
