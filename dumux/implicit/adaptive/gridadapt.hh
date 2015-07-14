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
#include "adaptationhelper.hh"
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
    typedef typename LeafGridView::template Codim<0>::Iterator LeafIterator;
    typedef typename GridView::IntersectionIterator LeafIntersectionIterator;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    typedef typename GET_PROP_TYPE(TypeTag, AdaptationIndicator) AdaptationIndicator;
    typedef typename GET_PROP_TYPE(TypeTag, AdaptationInitializationIndicator) AdaptationInitializationIndicator;

public:
    /*!
     * Constructor for h-adaptive simulations (adaptive grids)
     * @param problem The problem
     */
    ImplicitGridAdapt (Problem& problem)
        : adaptationHelper_(problem.gridView()),
          problem_(problem),
          adaptationIndicator_(problem),
          marked_(0),
          coarsened_(0)
    {
        levelMin_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MinLevel);
        levelMax_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
        adaptationInterval_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, AdaptationInterval);

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
     * AdaptationInitializationIndicator
     */
    void init()
    {
        adaptationIndicator_.init();

        if (!GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, EnableInitializationIndicator))
            return;

        AdaptationInitializationIndicator adaptationInitIndicator(problem_, adaptationIndicator_);

        int maxIter = 2*levelMax_;
        int iter = 0;
        while (iter <= maxIter)
        {
            adaptGrid(adaptationInitIndicator);

            if (!wasAdapted())
            {
                break;
            }

            int shouldInitialize = adaptationInitIndicator.initializeModel();
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
     * indicator (selected by the property AdaptationIndicator) and forwards to
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
        adaptGrid(adaptationIndicator_) ;
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

        // check for adaptation interval: Adapt only at certain time step indices
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

        /****  2b) Do pre-adaptation step    *****/
        problem_.grid().preAdapt();
        problem_.preAdapt();

        /****  3) Put primary variables in a map         *********/
        adaptationHelper_.storePrimVars(problem_);

        /****  4) Adapt Grid and size of variable vectors    *****/
        problem_.grid().adapt();

        // update mapper to new cell indices
        problem_.elementMapper().update();
        problem_.vertexMapper().update();

        // adapt size of vectors
        problem_.model().adaptVariableSize();

        /****  5) (Re-)construct primary variables to new grid **/
        adaptationHelper_.reconstructPrimVars(problem_);

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
    void markElements(Indicator& indicator)
    {
        typedef std::unordered_map<int, int> CoarsenMarkerType;
        CoarsenMarkerType coarsenMarker;

        for (LeafIterator eIt = problem_.gridView().template begin<0>();
             eIt!=problem_.gridView().template end<0>(); ++eIt)
        {
            // only mark non-ghost elements
            if (eIt->partitionType() != Dune::GhostEntity)
            {
                // refine?
                if (indicator.refine(*eIt) && eIt->level() < levelMax_)
                {
                    problem_.grid().mark( 1,  *eIt);
                    ++marked_;

                    // this also refines the neighbor elements
                    checkNeighborsRefine_(*eIt);
                }
                if (indicator.coarsen(*eIt) && eIt->hasFather())
                {
                    problem_.grid().mark( -1, *eIt );
                    ++coarsened_;
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

    const bool wasAdapted() const
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
    const int getMaxLevel() const
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
    const int getMinLevel() const
    {
        return levelMin_;
    }

    AdaptationIndicator& adaptationIndicator()
    {
        return adaptationIndicator_;
    }

    AdaptationIndicator& adaptationIndicator() const
    {
        return adaptationIndicator_;
    }

private:
    AdaptationHelper<TypeTag> adaptationHelper_;


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
        LeafIntersectionIterator isend = problem_.gridView().iend(entity);
        for(LeafIntersectionIterator is = problem_.gridView().ibegin(entity); is != isend; ++is)
        {
            if(!is->neighbor())
                continue;

            ElementPointer outside = is->outside();

            // only mark non-ghost elements
            if (outside->partitionType() == Dune::GhostEntity)
                continue;

            if ((outside->level() < levelMax_)
                && (outside->level() < entity.level()))
            {
                problem_.grid().mark(1, *outside);
                ++marked_;

                if(level != levelMax_)
                    checkNeighborsRefine_(*outside, ++level);
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
            for (LeafIterator eIt = leafGridView.template begin<0>();
                 eIt!=leafGridView.template end<0>(); ++eIt)
            {
                // only mark non-ghost elements
                if (eIt->partitionType() == Dune::GhostEntity)
                    continue;

                // run through all neighbor-cells (intersections)
                LeafIntersectionIterator isItend = leafGridView.iend(*eIt);
                for (LeafIntersectionIterator isIt = leafGridView.ibegin(*eIt); isIt!= isItend; ++isIt)
                {
                    const typename LeafIntersectionIterator::Intersection intersection = *isIt;
                    if(!intersection.neighbor())
                        continue;

                    ElementPointer outside =intersection.outside();
                    if (eIt.level()+maxLevelDelta<outside.level())
                    {
                        ElementPointer entity =eIt;
                        problem_.grid().mark( 1, *entity );
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
    AdaptationIndicator adaptationIndicator_;

    int marked_;
    int coarsened_;

    int levelMin_;
    int levelMax_;

    int adaptationInterval_;
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
    bool wasAdapted()
    {
        return false;
    }
    void setLevels(int, int)
    {}
    void setTolerance(int, int)
    {}
    const void setIndicator(const SolutionVector&,
                            const Scalar&, const Scalar&)
    {}
    ImplicitGridAdapt (Problem& problem)
    {}
};

}
#endif /* DUMUX_IMPLICIT_GRIDADAPT_HH */
