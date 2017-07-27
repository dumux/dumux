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
NEW_PROP_TAG(ElementBoundaryTypes);
NEW_PROP_TAG(ElementVolumeVariables);
NEW_PROP_TAG(FVElementGeometry);
NEW_PROP_TAG(ImplicitIsBox);
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
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
    };

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Vertex = typename GridView::Traits::template Codim<dim>::Entity;
    using CoordScalar = typename GridView::Grid::ctype;
    using ReferenceElement = typename Dune::ReferenceElement<CoordScalar, dim>;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;
    using AdaptionIndicator = typename GET_PROP_TYPE(TypeTag, AdaptionIndicator);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

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
        for (auto&& scv : scvs(fvGeometry))
        {
            auto source = problem_.source(element, fvGeometry, elemVolVars, scv);
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
    bool hasRefineBC_(const Element& element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf,
                      const ElementBoundaryTypes &elemBcTypes)
    {
        // Check for boundary conditions
        if(isBox)
        {
            // Dune::GeometryType geoType = element.geometry().type();
            // const ReferenceElement &refElement = ReferenceElements::general(geoType);
            // int fIdx = intersection.indexInInside();
            // int numFaceVerts = refElement.size(fIdx, 1, dim);
            // for (int faceVertexIdx = 0; faceVertexIdx < numFaceVerts; ++faceVertexIdx)
            // {
            //     int scvIdx = refElement.subEntity(fIdx, 1, faceVertexIdx, dim);
            //     BoundaryTypes bcTypes = elemBcTypes[scvIdx];
            //     problemBoundaryTypes_(bcTypes, element.template subEntity<dim>(scvIdx));
            //     int bfIdx = fvGeometry.boundaryFaceIndex(fIdx, faceVertexIdx);
            //     for (int i = 0; i < numEq; i++)
            //     {
            //         if(bcTypes.isDirichlet(i) && refineAtDirichletBC_)
            //             return true;
            //         if(bcTypes.isNeumann(i) && refineAtFluxBC_)
            //         {
            //             PrimaryVariables fluxes(0.0);
            //             problem_.solDependentNeumann(fluxes, element, fvGeometry,
            //                                      intersection, scvIdx, bfIdx,
            //                                      elemVolVars);
            //             if (fluxes.infinity_norm() > eps_)
            //                 return true;
            //         }
            //     }
            // }
            return false;
        }
        else
        {
            auto bcTypes = problem_.boundaryTypes(element, scvf);
            for (int i = 0; i < numEq; i++)
            {
                if(bcTypes.isDirichlet(i) && refineAtDirichletBC_)
                    return true;

                else if(bcTypes.isNeumann(i) && refineAtFluxBC_)
                {
                    auto fluxes = problem_.neumann(element, fvGeometry, elemVolVars, scvf);
                    if (fluxes.infinity_norm() > eps_)
                        return true;
                }
            }
        }
        return false;
    }

    // only actually call the problem method when it exists
    template <class T = TypeTag>
    BoundaryTypes problemBoundaryTypes_(const Element &element,
                                        const typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox),
                                        SubControlVolume>::type& scv) const
    {
        return problem_.boundaryTypes(element, scv);
    }
    template <class T = TypeTag>
    BoundaryTypes problemBoundaryTypes_(const Element &element,
                                        const typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox),
                                        SubControlVolume>::type& scv) const
    {
        return BoundaryTypes();
    }

public:
    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     */
    void calculateIndicator()
    {
        if (!enableInitializationIndicator_)
            return;

        // first adapt for boundary conditions and sources to get a good initial solution
        if (nextMaxLevel_ == maxAllowedLevel_)
            adaptionIndicator_.calculateIndicator();

        // prepare an indicator for refinement
        indicatorVector_.assign(problem_.gridView().size(0), 0);

        // 1) calculate Indicator -> min, maxvalues
        // loop over all leaf elements
        for (const auto& element : elements(problem_.gridView()))
        {
            int globalIdxI = problem_.elementMapper().index(element);
            int level = element.level();
            using std::max;
            using std::min;
            maxLevel_ = max(level, maxLevel_);

            if (level < minAllowedLevel_)
            {
                nextMaxLevel_ = min(max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                indicatorVector_[globalIdxI] = refineCell;
                continue;
            }

            if (!(refineAtSource_ || refineAtFluxBC_ || refineAtDirichletBC_))
                continue;

            // get the fvGeometry and elementVolVars needed for the bc and source interfaces
            auto fvGeometry = localView(problem_.model().globalFvGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(problem_.model().curGlobalVolVars());
            elemVolVars.bindElement(element, fvGeometry, problem_.model().curSol());

            // Check if we have to refine around a source term
            if (indicatorVector_[globalIdxI] != refineCell && refineAtSource_)
            {
                if(hasSource_(element, fvGeometry, elemVolVars))
                {
                    nextMaxLevel_ = min(max(level + 1, nextMaxLevel_), maxAllowedLevel_);
                    indicatorVector_[globalIdxI] = refineCell;
                    continue;
                }
            }

            // get the element boundary types
            ElementBoundaryTypes elemBcTypes;
            elemBcTypes.update(problem_, element, fvGeometry);

            // Check if we have to refine at the boundary
            if (indicatorVector_[globalIdxI] != refineCell && (refineAtDirichletBC_ || refineAtFluxBC_))
            {
                // Calculate the boundary indicator for all boundary intersections
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    if (scvf.boundary())
                    {
                        if(hasRefineBC_(element, fvGeometry, elemVolVars, scvf, elemBcTypes))
                        {
                            nextMaxLevel_ = min(max(level + 1, nextMaxLevel_), maxAllowedLevel_);
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
    bool refine(const Element& element) const
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
    bool coarsen(const Element& element) const
    { return false; }

    int maxLevel() const
    { return maxLevel_; }

    /*! \brief Initializes the adaption indicator class */
    void init() {}

    /*! \brief If the model needs to be initialized after adaption.
     *         We always need to initialize since the hasRefineBC_ method needs information from
               an up-to-date model.
     */
    constexpr bool initializeModel() const
    { return true; }

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
    std::vector<int> indicatorVector_;
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
    bool refine(const Element& element) const
    {
        return false;
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool coarsen(const Element& element) const
    {
        return false;
    }

    constexpr bool initializeModel() const
    {
        return false;
    }

    /*! \brief Initializes the adaption indicator class*/
    void init() {}

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
