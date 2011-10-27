/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle                                     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Base class for h-adaptive sequential models.
 */
#ifndef DUMUX_GIRDADAPT_HH
#define DUMUX_GIRDADAPT_HH


namespace Dumux
{

/*!
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))     Problem;


    //*******************************
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid                         Grid;
    typedef typename Grid::LeafGridView                     LeafGridView;
    typedef typename LeafGridView::template Codim<0>::Iterator LeafIterator;
    typedef typename GridView::IntersectionIterator         LeafIntersectionIterator;
    typedef typename Grid::template Codim<0>::Entity         Entity;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

public:
    /*!
     * Constructor for h-adaptive simulations (adaptive grids)
     * @param problem The problem
     * @param levelMin minimum refinement level
     * @param levelMax maximum refinement level
     */
    GridAdapt (Problem& problem, const int levelMin = 0,const int levelMax=1)
    : levelMin_(levelMin), levelMax_(levelMax), problem_(problem)
    {
        if (levelMin_ < 0)
            Dune::dgrave <<  __FILE__<< ":" <<__LINE__
            << " :  Dune cannot coarsen to gridlevels smaller 0! "<< std::endl;

        standardIndicator_ = true;
        refinetol_ = 0.05;
        coarsentol_ = 0.001;
    }

    /*!
     * @brief Standard method to adapt the grid
     *
     * This method is called in preTimeStep() of the problem if
     * adaptive grids are used in the simulation.
     *
     * It uses a standard procedure for adaptivity:
     * 1) Determine the refinement indicator
     * 2) Mark the elements
     * 3) Store primary variables in a map
     * 4) Adapt the grid, adapt variables sizes, update mappers
     * 5) Reconstruct primary variables, regain secondary variables
     */
    void adaptGrid()
    {
        // reset internal counter for marked elements
        marked_ = coarsened_ = 0;

        // prepare an indicator for refinement
        if(indicatorVector_.size()!=problem_.variables().gridSize())
        {
            indicatorVector_.resize(problem_.variables().gridSize());
            indicatorVector_ = -1e00;
        }

        /**** 1) determine refining parameter if standard is used ***/
        // if not, the indicatorVector and refinement Bounds have to
        // specified by the problem through setIndicator()
        if(standardIndicator_)
            indicator();

        /**** 2) mark elements according to indicator     *********/
        markElements(indicatorVector_, refineBound_, coarsenBound_);

        // abort if nothing in grid is marked
        if (marked_==0 && coarsened_ == 0)
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
//        problem_.grid().preAdapt();
        problem_.grid().adapt();

//        forceRefineRatio(1);

        // update mapper to new cell idice
        problem_.variables().elementMapper().update();

        // adapt size of vectors (
        problem_.variables().adaptVariableSize(problem_.variables().elementMapper().size());

        /****  5) (Re-)construct primary variables to new grid **/
        problem_.variables().reconstructPrimVars(problem_);
        // delete markers in grid
        problem_.grid().postAdapt();
        // adapt secondary variables
        problem_.pressureModel().updateMaterialLaws();

        // write out new grid
//        Dune::VTKWriter<LeafGridView> vtkwriter(leafView);
//        vtkwriter.write("latestgrid",Dune::VTKOptions::binaryappended);
        return;
    };

    /*!
     * Mark Elements for grid refinement according to applied Indicator
     * @param indicator Vector where the refinement indicator is stored
     * @param refineThreshold lower threshold where to refine
     * @param coarsenThreshold upper threshold where to coarsen
     * @return
     */
    int markElements(ScalarSolutionType &indicator, const double refineThreshold, const double coarsenThreshold)
    {
        for (LeafIterator eIt = problem_.gridView().template begin<0>();
                    eIt!=problem_.gridView().template end<0>(); ++eIt)
        {
            // Verfeinern?
            if (indicator[problem_.variables().elementMapper().map(*eIt)] > refineThreshold
                    && eIt.level()<levelMax_)
            {
                const Entity &entity =*eIt;
                problem_.grid().mark( 1, entity );
                ++marked_;

                // this also refines the neighbor elements
                checkNeighborsRefine_(entity);
            }
        }
        // coarsen
        for (LeafIterator eIt = problem_.gridView().template begin<0>();
                    eIt!=problem_.gridView().template end<0>(); ++eIt)
        {
            if (indicator[problem_.variables().elementMapper().map(*eIt)] < coarsenThreshold
                    && eIt.level()>levelMin_ && problem_.grid().getMark(*eIt) == 0)
            {
                // check if coarsening is possible
                bool coarsenPossible = true;
                LeafIntersectionIterator isend = problem_.gridView().iend(*eIt);
                for(LeafIntersectionIterator is = problem_.gridView().ibegin(*eIt); is != isend; ++is)
                {
                    if(!is->neighbor())
                        continue;

                    const Entity &outside = *is->outside();
                    if ((problem_.grid().getMark(outside) > 0)
                            || (outside.level()>eIt.level()))
                        coarsenPossible = false;
                }

                if(coarsenPossible)
                    {
                        problem_.grid().mark( -1, *eIt );
                        ++coarsened_;
                    }
            }
        }

        return marked_;
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
     * Sets minimum and maximum refinement tolerances
     *
     * @param coarsentol minimum tolerance when to coarsen
     * @param refinetol maximum tolerance when to refine
     */
    void setTolerance(Scalar coarsentol, Scalar refinetol)
    {
        if (coarsentol < 0. or refinetol < 0. )
            Dune::dgrave <<  __FILE__<< ":" <<__LINE__
            << " :  Tolerance levels out of meaningful bounds! "<< std::endl;
        if (coarsentol > refinetol )
            Dune::dgrave <<  __FILE__<< ":" <<__LINE__
            << " :  Check tolerance levels: coarsentol > refinetol! "<< std::endl;

        coarsentol_ = coarsentol;
        refinetol_ = refinetol;
    }
    /*!
     * Gets maximum refinement level
     *
     * @return levelMax_ maximum level for refinement
     */
    const int getMaxLevel() const
    {
        return levelMax_;
    }
    /*!
     * Gets minimum refinement level
     *
     * @return levelMin_ minimum level for coarsening
     */
    const int getMinLevel() const
    {
        return levelMin_;
    }
    /*!
     * @brief Adapter for external Refinement indicators.
     *
     * External indicators to refine/coarsen the grid can be set
     * from outside this container. Both indicator itsself as well as
     * the coarsening and refinement bounds have to be specified.
     * @param indicatorVector Vector holding indicator values
     * @param coarsenLowerBound bounding value where to be coarsened
     * @param refineUpperBound bounding value where to be refined
     */
    const void setIndicator(const ScalarSolutionType& indicatorVector,
            const Scalar& coarsenLowerBound, const Scalar& refineUpperBound)
    {
        // switch off usage of standard (saturation) indicator: by settting this,
        // the standard function indicator() called by adaptGrid() is disabled.
        standardIndicator_=false;

        indicatorVector_ = indicatorVector;
        refineBound_ = refineUpperBound;
        coarsenBound_ = coarsenLowerBound;
    }
private:
    /*!
     * @brief Simple standard indicator.
     *
     * Mehod computes the refinement and coarsening bounds through a
     * standard refinement criteria.
     */
    void indicator()
    {
        Scalar globalMax_(0.), globalMin_(0.);

        /**** determine refining parameter          *************/
        problem_.transportModel().indicatorSaturation(indicatorVector_, globalMin_=1e100, globalMax_=-1e100);
        Scalar globaldelta = globalMax_- globalMin_;
//            globaldelta = std::max(globaldelta,0.1);

        refineBound_ = refinetol_*globaldelta;
        coarsenBound_ = coarsentol_*globaldelta;

        return;
    }
    /*!
     * @brief Method ensuring the refinement ratio of 2:1
     *
     * For any given entity, a loop over the neighbors checks weather the
     * entities refinement would require that any of the neighbors has
     * to be refined, too.
     * This is done recursively over all levels of the grid.
     *
     * @param entity Entity of interest that is to be refined
     * @param level level of the refined entity: it is at least 1
     * @return true if everything was successful
     */
    bool checkNeighborsRefine_(const Entity &entity, int level = 1)
    {
        // this also refines the neighbor elements
        LeafIntersectionIterator isend = problem_.gridView().iend(entity);
        for(LeafIntersectionIterator is = problem_.gridView().ibegin(entity); is != isend; ++is)
        {
            if(!is->neighbor())
                continue;

            const Entity &outside =*is->outside();
            if ((outside.level()<levelMax_)
                    && (outside.level()<entity.level()))
                {
                    problem_.grid().mark(1, outside);
                    ++marked_;

                    if(level != levelMax_)
                        checkNeighborsRefine_(outside, ++level);
                }
        }
        return true;
    };


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
        LeafGridView leafView = problem_.gridView();
        // delete all existing marks
        problem_.grid().postAdapt();
        bool done;
        do
        {
            // run through all cells
            done=true;
            for (LeafIterator it = leafView.template begin<0>();
                        it!=leafView.template end<0>(); ++it)
            {
                // run through all neighbor-cells (intersections)
                LeafIntersectionIterator isend = leafView.iend(*it);
                for (LeafIntersectionIterator is = leafView.ibegin(*it); is!= isend; ++is)
                {
                    const typename LeafIntersectionIterator::Intersection intersection = *is;
                    if(!intersection.neighbor())
                        continue;

                    const Entity &outside =intersection.outside();
                    if (it.level()+maxLevelDelta<outside.level())
                            {
                            const Entity &entity =*it;
                            problem_.grid().mark( 1, entity );
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
    ScalarSolutionType indicatorVector_;
    bool standardIndicator_;
    Scalar refineBound_, coarsenBound_;
    int marked_, coarsened_;
    int levelMin_, levelMax_;
    Scalar refinetol_;
    Scalar coarsentol_;

    Problem& problem_;
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))     Problem;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

public:
    void adaptGrid()
    {};
    void setLevels(int, int)
    {};
    void setTolerance(int, int)
    {};
    const void setIndicator(const ScalarSolutionType&,
            const Scalar&, const Scalar&)
    {};
    GridAdapt (Problem& problem)
    {}
};

}
#endif /* DUMUX_GIRDADAPT_HH */
