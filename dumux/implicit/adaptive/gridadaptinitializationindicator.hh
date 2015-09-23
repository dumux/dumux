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

#include <dune/geometry/type.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include "gridadaptproperties.hh"

/**
 * @file
 * @brief  Class defining an initialization indicator for grid adaption
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
 * @brief  Class defining an initialization indicator for grid adaption
 *
 *  Uses the defined grid adaption indicator and further accounts for sources and boundaries.
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

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
    };

    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Traits::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid::ctype CoordScalar;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;
    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, AdaptionIndicator) AdaptionIndicator;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    enum { refineCell = 1 };




    /*! \brief Search for a source term
     *
     *  For every element we check if the element center or the element corners
     *  are inside a source zone with source value > 0.
     *
     *  \param element A grid element
     */
    bool hasSource_(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVolumeVariables& elemVolVars)
    {
        PrimaryVariables source(0.0);
        for (int scvIdx = 0; scvIdx < fvGeometry.numScv; scvIdx++)
        {
            problem_.solDependentSource(source, element, fvGeometry, scvIdx, elemVolVars);
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
    bool hasRefineBC_(ElementBoundaryTypes &elemBcTypes,
                      const Element& element,
                      const Intersection& intersection,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars)
    {
        // Check for boundary conditions
        if(isBox)
        {
            Dune::GeometryType geoType = element.geometry().type();
            const ReferenceElement &refElement = ReferenceElements::general(geoType);
            int fIdx = intersection.indexInInside();
            int numFaceVerts = refElement.size(fIdx, 1, dim);
            for (int faceVertexIdx = 0; faceVertexIdx < numFaceVerts; ++faceVertexIdx)
            {
                int scvIdx = refElement.subEntity(fIdx, 1, faceVertexIdx, dim);
                BoundaryTypes bcTypes = elemBcTypes[scvIdx];
                const VertexPointer v = element.template subEntity<dim>(scvIdx);
                problemBoundaryTypes_(bcTypes, *v);
                int bfIdx = fvGeometry.boundaryFaceIndex(fIdx, faceVertexIdx);
                for (int i = 0; i < numEq; i++)
                {
                    if(bcTypes.isDirichlet(i) && refineAtDirichletBC_)
                        return true;
                    if(bcTypes.isNeumann(i) && refineAtFluxBC_)
                    {
                        PrimaryVariables fluxes(0.0);
                        problem_.solDependentNeumann(fluxes, element, fvGeometry,
                                                 intersection, scvIdx, bfIdx,
                                                 elemVolVars);
                        if (fluxes.infinity_norm() > eps_)
                            return true;
                    }
                }
            }

        }
        else
        {
            BoundaryTypes bcTypes = elemBcTypes[0];
            problem_.boundaryTypes(bcTypes, intersection);
            int bfIdx = intersection.indexInInside();
            for (int i = 0; i < numEq; i++)
            {
                if(bcTypes.isDirichlet(i) && refineAtDirichletBC_)
                    return true;
                if(bcTypes.isNeumann(i) && refineAtFluxBC_)
                {
                    PrimaryVariables fluxes(0.0);
                    problem_.solDependentNeumann(fluxes, element, fvGeometry,
                                                 intersection, 0, bfIdx,
                                                 elemVolVars);
                    if (fluxes.infinity_norm() > eps_)
                        return true;
                }
            }
        }
        return false;
    }

    // only actually call the problem method when it exists
    template <class T = TypeTag>
    void problemBoundaryTypes_(BoundaryTypes& bcTypes,
                              const typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox), Vertex>::type& v) const
    {
        problem_.boundaryTypes(bcTypes, v);
    }
    template <class T = TypeTag>
    void problemBoundaryTypes_(BoundaryTypes& bcTypes,
                              const typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox), Vertex>::type& v) const
    {}


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

            if (!(refineAtSource_ || refineAtFluxBC_ || refineAtDirichletBC_))
            	continue;

            // get the fvGeometry and elementVolVars needed for the bc and source interfaces
            FVElementGeometry fvGeometry;
            ElementVolumeVariables elemVolVars;
            fvGeometry.update(problem_.gridView(), *eIt);
            elemVolVars.update(problem_, *eIt, fvGeometry, false /* oldSol? */);

            // Check if we have to refine around a source term
            if (indicatorVector_[globalIdxI] != refineCell && refineAtSource_)
            {
                if(hasSource_(*eIt, fvGeometry, elemVolVars))
                {
                    nextMaxLevel_ = std::min(std::max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                    indicatorVector_[globalIdxI] = refineCell;
                    continue;
                }
            }

            // get the element boundary types
            ElementBoundaryTypes bcTypes;
            bcTypes.update(problem_, *eIt, fvGeometry);

            // Check if we have to refine at the boundary
            if (indicatorVector_[globalIdxI] != refineCell && (refineAtDirichletBC_ || refineAtFluxBC_))
            {
                // Calculate the boundary indicator for all boundary intersections
                const IntersectionIterator isItend = problem_.gridView().iend(*eIt);
                for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItend; ++isIt)
                {
                    if (isIt->boundary())
                    {
                        if(hasRefineBC_(bcTypes, *eIt, *isIt, fvGeometry, elemVolVars))
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
        return false;
    }

    int maxLevel()
    {
        return maxLevel_;
    }

    /*! \brief Initializes the adaption indicator class */
    void init()
    {};

    /*! \brief If the model needs to be initialized after adaption.
     *         We always need to initialize since the hasRefineBC_ method needs information from
               an up-to-date model.
     */
    bool initializeModel()
    {
        return true;
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
        problem_(problem), adaptionIndicator_(adaptionIndicator), maxLevel_(0), nextMaxLevel_(0), eps_(1e-30)
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
    Scalar eps_;
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
