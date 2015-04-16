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
#ifndef DUMUX_IMPLICIT_GRIDADAPTINITIALIZATIONINDICATOR_HH
#define DUMUX_IMPLICIT_GRIDADAPTINITIALIZATIONINDICATOR_HH

#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include "gridadaptproperties.hh"

/**
 * @file
 * @brief  Class defining an initialization indicator for grid adaption
 */
namespace Dumux
{
/*!\ingroup ImplicitGridAdaptInitializationIndicator
 * @brief  Class defining an initialization indicator for grid adaption
 *
 *  Uses the defined grid adaptation indicator and further accounts for sources and boundaries.
 *  Only for grid initialization!
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class ImplicitGridAdaptInitializationIndicator
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, AdaptionIndicator) AdaptionIndicator;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
    };

    enum
    {
        refineCell = 1,
        coarsenCell = -1
    };

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dim-1> LocalPositionFace;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    /*! \brief Hierarchical search for a source term
     *
     *  For every element we do virtual refinement until maxAllowedLevel
     *  and check if we have a source term at the element center. This
     *  is necessary as an element can also be only partly in a source zone.
     *
     *  \param source A primary variables vector that will return the source values
     *  \param element A grid element
     */
    void virtualHierarchicSourceSearch_(PrimaryVariables &source, const Element& element)
    {
        int level = element.level();
        const auto geometry = element.geometry();
        
        if (level == maxAllowedLevel_)
        {
            GlobalPosition globalPos = geometry.center();
            problem_.sourceAtPos(source, globalPos);
            return;
        }

        // get the number of check points in each dimension
        unsigned int numRefine = maxAllowedLevel_ - level;
        int numCheckCoords = 1 << numRefine;

        // the local position of the check point as we do this on the reference element
        LocalPosition localPos(0.0);
        GlobalPosition globalPosCheck(0.0);

        // we check for a source in the midpoint
        Scalar halfInterval = (1.0/double(numCheckCoords))/2.0;

        // TODO this only works correctly for cubes!
        
        PrimaryVariables sourceCheck(0.0); 
        // use a switch statement to let the compiler do easy optimization
        switch(dim)
        {
        case 3:
            for (int i = 1; i <= numCheckCoords; i++)
                for (int j = 1; i <= numCheckCoords; i++)
                    for (int k = 1; i <= numCheckCoords; i++)
                    {
                        localPos[0] = double(i)/double(numCheckCoords) - halfInterval;
                        localPos[1] = double(j)/double(numCheckCoords) - halfInterval;
                        localPos[2] = double(k)/double(numCheckCoords) - halfInterval;
                        globalPosCheck = geometry.global(localPos);
                        problem_.sourceAtPos(sourceCheck, globalPosCheck);

                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            if (std::abs(sourceCheck[eqIdx]) > std::abs(source[eqIdx]))
                                source[eqIdx] = sourceCheck[eqIdx];
                        }
                    }
            break;
        case 2:
            for (int i = 1; i <= numCheckCoords; i++)
                for (int j = 1; i <= numCheckCoords; i++)
                {
                    localPos[0] = double(i)/double(numCheckCoords) - halfInterval;
                    localPos[1] = double(j)/double(numCheckCoords) - halfInterval;
                    globalPosCheck = geometry.global(localPos);
                    problem_.sourceAtPos(sourceCheck, globalPosCheck);

                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        if (std::abs(sourceCheck[eqIdx]) > std::abs(source[eqIdx]))
                            source[eqIdx] = sourceCheck[eqIdx];
                    }
                }
            break;
        case 1:
            for (int i = 1; i <= numCheckCoords; i++)
            {
                localPos[0] = double(i)/double(numCheckCoords) - halfInterval;
                globalPosCheck = geometry.global(localPos);
                problem_.sourceAtPos(sourceCheck, globalPosCheck);

                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    if (std::abs(sourceCheck[eqIdx]) > std::abs(source[eqIdx]))
                        source[eqIdx] = sourceCheck[eqIdx];
                }
            }
            break;
        }
    }

    /*! \brief Hierarchical search for the value of Neumann boundary condition
     *
     *  For every intersection we do virtual refinement until maxAllowedLevel
     *  and check which boundary condition is defined on the intersection center. This
     *  is necessary as an element can partly have Neumann boundary conditions.
     *
     *  \param bcTypes The boundary condition types 
     *  \param values The value of the boundary condition. Returns the Neumann flux values
     *  \param element A grid element
     *  \param intersection The boundary intersection
     */
    void virtualHierarchicBCSearch_(BoundaryTypes &bcTypes, PrimaryVariables &values, const Element& element, const Intersection& intersection)
    {
        int level = element.level();
        const auto isGeometry = intersection.geometry();

        if (level == maxAllowedLevel_ || dim==1)
        {
            GlobalPosition globalPos = isGeometry.center();
            problem_.boundaryTypesAtPos(bcTypes, globalPos);

            if (refineAtFluxBC_)
            {
                for (int i = 0; i < numEq; i++)
                {
                    if (bcTypes.isNeumann(i))
                    {
                        PrimaryVariables fluxes;
                        problem_.neumannAtPos(fluxes, globalPos);
                        values[i] += std::abs(fluxes[i]);
                    }
                }
            }
            return;
        }

        // get the number of check points in each dimension
        unsigned int numRefine = maxAllowedLevel_ - level;
        int numCheckCoords = 1 << numRefine;

        // the local position of the check point 
        // as we do this on the reference codim 1 element
        LocalPositionFace localPos(0.0);
        GlobalPosition globalPosCheck(0.0);
        
        // we check for the boundary condition type in the midpoint
        Scalar halfInterval = (1.0/double(numCheckCoords))/2.0;

        // TODO this only works correctly for cubes!

        PrimaryVariables fluxCheck(0.0);
        // use a switch statement to let the compiler do easy optimization
        switch(dim-1)
        {
        case 2:
            for (int i = 1; i <= numCheckCoords; i++)
                for (int j = 1; i <= numCheckCoords; i++)
                {
                    localPos[0] = double(i)/double(numCheckCoords) - halfInterval;
                    localPos[1] = double(j)/double(numCheckCoords) - halfInterval;
                    globalPosCheck = isGeometry.global(localPos);
                    problem_.boundaryTypesAtPos(bcTypes, globalPosCheck);

                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        if (refineAtDirichletBC_ && bcTypes.isDirichlet(eqIdx))
                        {
                            return;
                        }
                        else if (refineAtFluxBC_ && bcTypes.isNeumann(eqIdx))
                        {
                            problem_.neumannAtPos(fluxCheck, globalPosCheck);
                            values[eqIdx] += std::abs(fluxCheck[eqIdx]);
                        }
                    }
                }
            break;
        case 1:
            for (int i = 1; i <= numCheckCoords; i++)
            {
                localPos[0] = double(i)/double(numCheckCoords) - halfInterval;
                globalPosCheck = isGeometry.global(localPos);
                    problem_.boundaryTypesAtPos(bcTypes, globalPosCheck);

                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    if (refineAtDirichletBC_ && bcTypes.isDirichlet(eqIdx))
                    {
                        return;
                    }
                    else if (refineAtFluxBC_ && bcTypes.isNeumann(eqIdx))
                    {
                        problem_.neumannAtPos(fluxCheck, globalPosCheck);
                        values[eqIdx] += std::abs(fluxCheck[eqIdx]);
                    }
                }
            }
            break;
        }
    }


public:
    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     */
    void calculateIndicator()
    {
        if (!enableInitializationIndicator_)
            return;

        //First adapt for boundary conditions and sources to get a good initial solution
        if (nextMaxLevel_ == maxAllowedLevel_)
            adaptionIndicator_.calculateIndicator();

        // prepare an indicator for refinement
        indicatorVector_.resize(problem_.model().numDofs());

        // set the default to coarsen
        indicatorVector_ = coarsenCell;

        // 1) calculate Indicator -> min, maxvalues
        // loop over all leaf elements
        const ElementIterator eEndIt = problem_.gridView().template end<0>();
        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eEndIt; ++eIt)
        {

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        	int globalIdxI = problem_.elementMapper().index(*eIt);
#else
        	int globalIdxI = problem_.elementMapper().map(*eIt);
#endif
            int level = eIt->level();
            maxLevel_ = std::max(level, maxLevel_);

            if (level < minAllowedLevel_)
            {
                nextMaxLevel_ = std::min(std::max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                indicatorVector_[globalIdxI] = refineCell;
                continue;
            }

            // Check if we have to refine around a source term
            if (refineAtSource_)
            {
                PrimaryVariables source(0.0);
                virtualHierarchicSourceSearch_(source, *eIt);
                for (int i = 0; i < numEq; i++)
                {
                    // If source term was found mark cell for refinement
                    if (std::abs(source[i]) > 1e-30)
                    {
                        nextMaxLevel_ = std::min(std::max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                        indicatorVector_[globalIdxI] = refineCell;
                        break;
                    }
                }
            }

            // Check if we have to refine at the boundary
            if (indicatorVector_[globalIdxI] != refineCell && (refineAtDirichletBC_ || refineAtFluxBC_))
            {
                // Calculate the boundary indicator for all boundary intersections
                const IntersectionIterator isItend = problem_.gridView().iend(*eIt);
                for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItend; ++isIt)
                {
                    if (isIt->boundary() && indicatorVector_[globalIdxI] != refineCell)
                    {
                        BoundaryTypes bcTypes;
                        PrimaryVariables values(0.0);
                        virtualHierarchicBCSearch_(bcTypes, values, *eIt, *isIt);

                        for (int i = 0; i < numEq; i++)
                        {
                            if (bcTypes.isDirichlet(i) && refineAtDirichletBC_)
                            {
                                nextMaxLevel_ = std::min(std::max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                                indicatorVector_[globalIdxI] = refineCell;
                                break;
                            }
                            if (std::abs(values[i]) > 1e-30)
                            {
                                nextMaxLevel_ = std::min(std::max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                                indicatorVector_[globalIdxI] = refineCell;
                                break;
                            } 
                        }
                    }
                }
            }
        }
    }

    /*! \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool refine(const Element& element)
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int idx = problem_.elementMapper().index(element);
#else
        int idx = problem_.elementMapper().map(element);
#endif
        if (indicatorVector_[idx] == refineCell)
            return true;
        else if (maxLevel_ == maxAllowedLevel_)
            return adaptionIndicator_.refine(element);
        else
            return false;
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool coarsen(const Element& element)
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int idx = problem_.elementMapper().index(element);
#else
        int idx = problem_.elementMapper().map(element);
#endif
        if (indicatorVector_[idx] == coarsenCell && maxLevel_ < maxAllowedLevel_)
            return true;
        else if (indicatorVector_[idx] == coarsenCell && !adaptionIndicator_.refine(element))
            return true;
        else
            return false;
    }

    int maxLevel()
    {
        return maxLevel_;
    }

    /*! \brief Initializes the adaption indicator class */
    void init()
    {};

    bool initializeModel()
    {
        return nextMaxLevel_ == maxAllowedLevel_;
    }

    /*! \brief Constructs a GridAdaptionIndicator instance
     *
     * This standard indicator is based on the saturation gradient. It checks the local gradient
     * compared to the maximum global gradient. The indicator is compared locally to a
     * refinement/coarsening threshold to decide whether a cell should be marked for refinement
     * or coarsening or should not be adapted.
     *
     * \param problem The problem object
     * \param adaptionIndicator Indicator whether a be adapted
     */
    ImplicitGridAdaptInitializationIndicator(Problem& problem, AdaptionIndicator& adaptionIndicator):
        problem_(problem), adaptionIndicator_(adaptionIndicator), maxLevel_(0), nextMaxLevel_(0)
    {
        minAllowedLevel_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MinLevel);
        maxAllowedLevel_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
        enableInitializationIndicator_ = GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, EnableInitializationIndicator);
        refineAtDirichletBC_ = GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, RefineAtDirichletBC);
        refineAtFluxBC_ = GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, RefineAtFluxBC);
        refineAtSource_ = GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, RefineAtSource);

        if (!refineAtDirichletBC_ && !refineAtFluxBC_ && !refineAtSource_)
        {
            nextMaxLevel_ = maxAllowedLevel_;
            maxLevel_ = maxAllowedLevel_;
        }
    }

private:
    Problem& problem_;
    AdaptionIndicator& adaptionIndicator_;
    Dune::DynamicVector<int> indicatorVector_;
    int maxLevel_;
    int nextMaxLevel_;
    int minAllowedLevel_;
    int maxAllowedLevel_;
    bool enableInitializationIndicator_;
    bool refineAtDirichletBC_;
    bool refineAtFluxBC_;
    bool refineAtSource_;
};


/*!\ingroup IMPES
 * @brief  Class defining a start indicator for grid adaption
 *
 *Default implementation
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class ImplicitGridAdaptInitializationIndicatorDefault
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, AdaptionIndicator) AdaptionIndicator;

public:
    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     */
    void calculateIndicator()
    {}

    /*! \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool refine(const Element& element)
    {
        return false;
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool coarsen(const Element& element)
    {
        return false;
    }

    bool initializeModel()
    {
        return false;
    }

    /*! \brief Initializes the adaption indicator class*/
    void init()
    {};

    /*! \brief Constructs a GridAdaptionIndicator for initialization of an adaptive grid
     *
     * Default implementation
     *
     * \param problem The problem object
     * \param adaptionIndicator Indicator whether a be adapted
     */
    ImplicitGridAdaptInitializationIndicatorDefault(Problem& problem, AdaptionIndicator& adaptionIndicator)
    {}
};

}
#endif
