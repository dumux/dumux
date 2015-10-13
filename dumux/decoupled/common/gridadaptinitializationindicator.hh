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
#ifndef DUMUX_GRIDADAPTINITIALIZATIONINDICATOR_HH
#define DUMUX_GRIDADAPTINITIALIZATIONINDICATOR_HH

#include "decoupledproperties.hh"

#include <dune/common/dynvector.hh>
/**
 * @file
 * @brief  Class defining a start indicator for grid adaption
 */
namespace Dumux
{
/*!\ingroup IMPES
 * @brief  Class defining a start indicator for grid adaption
 *
 *  Uses the defined grid adaptation indicator and further accounts for sources and boundaries.
 *  Only for grid initialization!
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridAdaptInitializationIndicator
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
            numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
        };

    enum
        {
            refineCell = 1,
            coarsenCell = -1
        };

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dim-1> LocalPositionFace;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    void virtualHierarchicSourceSearch_(PrimaryVariables &source, const Element& element)
    {
        int level = element.level();

        if (level == maxAllowedLevel_)
        {
            GlobalPosition globalPos = element.geometry().center();
            problem_.sourceAtPos(source, globalPos);

            return;
        }

        unsigned int numRefine = maxAllowedLevel_ - level;
        int numCheckCoords = pow(2, numRefine);

        LocalPosition localPos(0.0);
        GlobalPosition globalPosCheck(0.0);
        Scalar halfInterval = (1.0/double(numCheckCoords))/2.;

        PrimaryVariables sourceCheck(0.0);

        for (int i = 1; i <= numCheckCoords; i++)
        {
            for (int j = 1; j <= numCheckCoords; j++)
            {
                localPos[0] = double(i)/double(numCheckCoords) - halfInterval;
                localPos[1] = double(j)/double(numCheckCoords) - halfInterval;
                if (dim == 2)
                {
                    globalPosCheck = element.geometry().global(localPos);
                    problem_.sourceAtPos(sourceCheck, globalPosCheck);

                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        if (std::abs(sourceCheck[eqIdx]) > std::abs(source[eqIdx]))
                        {
                            source[eqIdx] = sourceCheck[eqIdx];
                        }
                    }
                }
                else if (dim == 3)
                {
                    for (int k = 1; k <= numCheckCoords; k++)
                    {
                        localPos[2] = double(k)/double(numCheckCoords) - halfInterval;
                        globalPosCheck = element.geometry().global(localPos);
                        problem_.sourceAtPos(sourceCheck, globalPosCheck);

                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            if (std::abs(sourceCheck[eqIdx]) > std::abs(source[eqIdx]))
                            {
                                source[eqIdx] = sourceCheck[eqIdx];
                            }
                        }
                    }
                }
            }
        }
    }

    void virtualHierarchicBCSearch_(BoundaryTypes &bcTypes, PrimaryVariables &values, const Element& element, const Intersection& intersection)
    {
        int level = element.level();

        if (level == maxAllowedLevel_)
        {
            GlobalPosition globalPos = intersection.geometry().center();
            problem_.boundaryTypesAtPos(bcTypes, globalPos);

            if (refineAtFluxBC_)
            {
                for (int i = 0; i < numEq; i++)
                {
                    if (bcTypes.isNeumann(i))
                    {
                        PrimaryVariables fluxes;
                        problem_.neumannAtPos(fluxes, globalPos);

                        values += fluxes;
                    }
                }
            }
            return;
        }

        unsigned int numRefine = maxAllowedLevel_ - level;
        int numCheckCoords = pow(2, numRefine);

        LocalPositionFace localPos(0.0);
        GlobalPosition globalPosCheck(0.0);
        Scalar halfInterval = (1.0/double(numCheckCoords))/2.;

        PrimaryVariables fluxCheck(0.0);

        for (int i = 1; i <= numCheckCoords; i++)
        {
            localPos[0] = double(i)/double(numCheckCoords) - halfInterval;
            if (dim == 2)
            {
                globalPosCheck = intersection.geometry().global(localPos);
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

                        values += fluxCheck;
                    }
                }
            }
            else if (dim == 3)
            {
                for (int k = 1; k <= numCheckCoords; k++)
                {
                    localPos[1] = double(k)/double(numCheckCoords) - halfInterval;
                    globalPosCheck = intersection.geometry().global(localPos);

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

                            values += fluxCheck;
                        }

                    }
                }
            }
        }
    }


public:
    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     */
    void calculateIndicator()
    {
        //First Adapt for boundary conditions and sources to get a correct pressure solution
        if (nextMaxLevel_ == maxAllowedLevel_)
            adaptionIndicator_.calculateIndicator();

        // prepare an indicator for refinement
        indicatorVector_.resize(problem_.variables().cellDataGlobal().size());

        indicatorVector_ = coarsenCell;

        if (!enableInitializationIndicator_)
            return;

        ElementIterator eEndIt = problem_.gridView().template end<0>();
        // 1) calculate Indicator -> min, maxvalues
        // Schleife über alle Leaf-Elemente
        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eEndIt;
             ++eIt)
        {
            int globalIdxI = problem_.variables().index(*eIt);

            int level = eIt->level();
            maxLevel_ = std::max(level, maxLevel_);

            if (level < minAllowedLevel_)
            {
                nextMaxLevel_ = std::min(std::max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                indicatorVector_[globalIdxI] = refineCell;
                continue;
            }

            if (refineAtSource_)
            {
                PrimaryVariables source(0.0);
                virtualHierarchicSourceSearch_(source, *eIt);
                for (int i = 0; i < numEq; i++)
                {
                    if (std::abs(source[i]) > 1e-10)
                    {
                        nextMaxLevel_ = std::min(std::max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                        indicatorVector_[globalIdxI] = refineCell;
                        break;
                    }
                }
            }

            if (indicatorVector_[globalIdxI] != refineCell && (refineAtDirichletBC_ || refineAtFluxBC_))
            {
                // Berechne Verfeinerungsindikator an allen Zellen
                IntersectionIterator isItend = problem_.gridView().iend(*eIt);
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
                        }
                        for (int j = 0; j < numPhases; j++)
                        {
                            if (std::abs(values[j]) > 1e-10)
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
        int idx = problem_.elementMapper().index(element);

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
        int idx = problem_.elementMapper().index(element);

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
    {}

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
    GridAdaptInitializationIndicator(Problem& problem, AdaptionIndicator& adaptionIndicator):
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
}
#endif
