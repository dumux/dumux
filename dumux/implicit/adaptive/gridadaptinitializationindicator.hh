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
 * @brief  Class defining an initialization indicator for grid adaptation
 */
namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(BoundaryTypes);
NEW_PROP_TAG(PrimaryVariables);
NEW_PROP_TAG(NumEq);
}

/*!\ingroup ImplicitGridAdaptInitializationIndicator
 * @brief  Class defining an initialization indicator for grid adaptation
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

    typedef typename GET_PROP_TYPE(TypeTag, AdaptationIndicator) AdaptationIndicator;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
    };

    enum { refineCell = 1 };

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dim-1> LocalPositionFace;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    /*! \brief Search for a source term
     *
     *  For every element we check if the element center or the element corners
     *  are inside a source zone with source value > 0.
     *
     *  \param element A grid element
     */
    bool hasSource_(const Element& element)
    {
        const auto geometry = element.geometry();
        PrimaryVariables source(0.0);
        const GlobalPosition &globalPos = geometry.center();

        // Check if the midpoint is in a source zone
        problem_.sourceAtPos(source, globalPos);
        if (source.infinity_norm() > eps_)
            return true;

        // Check if a corner is on a source zone
        for (int vIdx = 0; vIdx < geometry.corners(); ++vIdx)
        {
            source = 0.0;
            problem_.sourceAtPos(source, geometry.corner(vIdx));
            if (source.infinity_norm() > eps_)
                return true;
        }
        return false;
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
    bool hasRefineBC_(BoundaryTypes &bcTypes, const Element& element, const Intersection& intersection)
    {
        const auto isGeometry = intersection.geometry();
        const GlobalPosition &globalPos = isGeometry.center();

        // Check if the midpoint has matching boundary condition
        problem_.boundaryTypesAtPos(bcTypes, globalPos);
        for (int i = 0; i < numEq; i++)
        {
            if(bcTypes.isDirichlet(i) && refineAtDirichletBC_)
                return true;
            if(bcTypes.isNeumann(i) && refineAtFluxBC_)
            {
                PrimaryVariables fluxes(0.0);
                problem_.neumannAtPos(fluxes, globalPos);
                if (fluxes.infinity_norm() > eps_)
                    return true;
            }
        }

        // Check if a corner has a matching boundary condition
        for (int vIdx = 0; vIdx < isGeometry.corners(); ++vIdx)
        {
            problem_.boundaryTypesAtPos(bcTypes, isGeometry.corner(vIdx));
            for (int i = 0; i < numEq; i++)
            {
                if(bcTypes.isDirichlet(i) && refineAtDirichletBC_)
                    return true;
                if(bcTypes.isNeumann(i) && refineAtFluxBC_)
                {
                    PrimaryVariables fluxes(0.0);
                    problem_.neumannAtPos(fluxes, isGeometry.corner(vIdx));
                    if (fluxes.infinity_norm() > eps_)
                        return true;
                }
            }
        }
        return false;
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
            adaptationIndicator_.calculateIndicator();

        // prepare an indicator for refinement
        indicatorVector_.resize(problem_.gridView().size(0));
        indicatorVector_ = 0;

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
            if (indicatorVector_[globalIdxI] != refineCell && refineAtSource_)
            {
                if(hasSource_(*eIt))
                {
                    nextMaxLevel_ = std::min(std::max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                    indicatorVector_[globalIdxI] = refineCell;
                    continue;
                }
            }

            // Check if we have to refine at the boundary
            if (indicatorVector_[globalIdxI] != refineCell && (refineAtDirichletBC_ || refineAtFluxBC_))
            {
                // Calculate the boundary indicator for all boundary intersections
                const IntersectionIterator isItend = problem_.gridView().iend(*eIt);
                for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItend; ++isIt)
                {
                    if (isIt->boundary())
                    {
                        BoundaryTypes bcTypes;
                        if(hasRefineBC_(bcTypes, *eIt, *isIt))
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
            return adaptationIndicator_.refine(element);
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
        return false;
    }

    int maxLevel()
    {
        return maxLevel_;
    }

    /*! \brief Initializes the adaptation indicator class */
    void init()
    {};

    bool initializeModel()
    {
        return nextMaxLevel_ == maxAllowedLevel_;
    }

    /*! \brief Constructs a GridAdaptationIndicator instance
     *
     * This standard indicator is based on the saturation gradient. It checks the local gradient
     * compared to the maximum global gradient. The indicator is compared locally to a
     * refinement/coarsening threshold to decide whether a cell should be marked for refinement
     * or coarsening or should not be adapted.
     *
     * \param problem The problem object
     * \param adaptationIndicator Indicator whether a be adapted
     */
    ImplicitGridAdaptInitializationIndicator(Problem& problem, AdaptationIndicator& adaptationIndicator):
        problem_(problem), adaptationIndicator_(adaptationIndicator), maxLevel_(0), nextMaxLevel_(0), eps_(1e-30)
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
    AdaptationIndicator& adaptationIndicator_;
    Dune::DynamicVector<int> indicatorVector_;
    int maxLevel_;
    int nextMaxLevel_;
    int minAllowedLevel_;
    int maxAllowedLevel_;
    bool enableInitializationIndicator_;
    bool refineAtDirichletBC_;
    bool refineAtFluxBC_;
    bool refineAtSource_;
    Scalar eps_;
};


/*!\ingroup IMPES
 * @brief  Class defining a start indicator for grid adaptation
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
    typedef typename GET_PROP_TYPE(TypeTag, AdaptationIndicator) AdaptationIndicator;

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

    /*! \brief Initializes the adaptation indicator class*/
    void init()
    {};

    /*! \brief Constructs a GridAdaptationIndicator for initialization of an adaptive grid
     *
     * Default implementation
     *
     * \param problem The problem object
     * \param adaptationIndicator Indicator whether a be adapted
     */
    ImplicitGridAdaptInitializationIndicatorDefault(Problem& problem, AdaptationIndicator& adaptationIndicator)
    {}
};

}
#endif
