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
#ifndef DUMUX_GRIDADAPTIONINDICATOR2PLOCAL_HH
#define DUMUX_GRIDADAPTIONINDICATOR2PLOCAL_HH

#include <dune/common/version.hh>
#include <dumux/decoupled/common/impetproperties.hh>
#include <dumux/decoupled/2p/2pproperties.hh>

/**
 * \file
 * \brief  Class defining a standard, saturation dependent indicator for grid adaption
 */
namespace Dumux
{
/*!\ingroup IMPES
 * \brief  Class defining a standard, saturation dependent indicator for grid adaption
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridAdaptionIndicator2PLocal
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        sw = Indices::saturationW,
        sn = Indices::saturationNw
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

public:
    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     * This standard indicator is based on the saturation gradient.
     */
    void calculateIndicator()
    {
        // prepare an indicator for refinement
        if(indicatorVector_.size() != problem_.variables().cellDataGlobal().size())
        {
            indicatorVector_.resize(problem_.variables().cellDataGlobal().size());
        }
        indicatorVector_ = -1e100;

        Scalar globalMax = -1e100;
        Scalar globalMin = 1e100;

        Scalar maxLocalDelta = 0;

        ElementIterator eEndIt = problem_.gridView().template end<0>();
        // 1) calculate Indicator -> min, maxvalues
        // loop over all leaf-elements
        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eEndIt;
                ++eIt)
        {
            // calculate minimum and maximum saturation
            // index of the current leaf-elements
            int globalIdxI = problem_.variables().index(*eIt);

            if (refineAtSource_)
            {
            PrimaryVariables source(0.0);
            problem_.source(source, *eIt);
            for (int i = 0; i < 2; i++)
            {
                if (std::abs(source[i]) > 1e-10)
                {
                    indicatorVector_[globalIdxI] = 10;
                    break;
                }
            }
            }

            if (indicatorVector_[globalIdxI] == 10)
            {
                continue;
            }

            Scalar satI = 0.0;
            switch (saturationType_)
            {
            case sw:
                satI = problem_.variables().cellData(globalIdxI).saturation(wPhaseIdx);
            break;
            case sn:
                satI = problem_.variables().cellData(globalIdxI).saturation(nPhaseIdx);
                break;
            }

            globalMin = std::min(satI, globalMin);
            globalMax = std::max(satI, globalMax);

            // calculate refinement indicator in all cells
            IntersectionIterator isItend = problem_.gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItend; ++isIt)
            {
                if (indicatorVector_[globalIdxI] == 10)
                {
                    break;
                }

                const typename IntersectionIterator::Intersection &intersection = *isIt;
                // exit, if it is not a neighbor
                if (isIt->boundary())
                {
                    BoundaryTypes bcTypes;
                    problem_.boundaryTypes(bcTypes, *isIt);

                    for (int i = 0; i < 2; i++)
                    {
                        if (bcTypes.isNeumann(i))
                        {
                            PrimaryVariables flux(0.0);
                            problem_.neumann(flux, *isIt);

                            bool fluxBound = false;
                            for (int j = 0; j < 2; j++)
                            {
                            if (std::abs(flux[j]) > 1e-10)
                            {
                                if (refineAtFluxBC_)
                                {
                                    indicatorVector_[globalIdxI] = 10;
                                    fluxBound = true;
                                    break;
                                }
                            }
                            }
                            if (fluxBound)
                                break;
                        }
                        else if (bcTypes.isDirichlet(i))
                        {
                            if (refineAtDirichletBC_)
                            {
                                indicatorVector_[globalIdxI] = 10;
                                break;
                            }
                        }
                    }
                }
                else
                {
                    // get element pointer from neighbor
                    auto outside = intersection.outside();
                    int globalIdxJ = problem_.variables().index(outside);

                    // treat each intersection only from one side
                    if (eIt->level() > outside.level()
                            || (eIt->level() == outside.level() && globalIdxI < globalIdxJ))
                    {
                        Scalar satJ = 0.;
                        switch (saturationType_)
                        {
                        case sw:
                            satJ = problem_.variables().cellData(globalIdxJ).saturation(wPhaseIdx);
                            break;
                        case sn:
                            satJ = problem_.variables().cellData(globalIdxJ).saturation(nPhaseIdx);
                            break;
                        }

                        Scalar localdelta = std::abs(satI - satJ);
                        indicatorVector_[globalIdxI][0] = std::max(indicatorVector_[globalIdxI][0], localdelta);
                        indicatorVector_[globalIdxJ][0] = std::max(indicatorVector_[globalIdxJ][0], localdelta);

                        maxLocalDelta = std::max(maxLocalDelta, localdelta);
                    }
                }
            }
        }

        Scalar globaldelta = globalMax - globalMin;

        if (maxLocalDelta > 0.)
        {
            indicatorVector_ /= maxLocalDelta;
        }
        if (globaldelta > 0.)
        {
            maxLocalDelta /= globaldelta;
        }

        refineBound_ = refinetol_*maxLocalDelta;
        coarsenBound_ = coarsentol_*maxLocalDelta;
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
        return (indicatorVector_[problem_.elementMapper().index(element)] > refineBound_);
#else
        return (indicatorVector_[problem_.elementMapper().map(element)] > refineBound_);
#endif
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
        return (indicatorVector_[problem_.elementMapper().index(element)] < coarsenBound_);
#else
        return (indicatorVector_[problem_.elementMapper().map(element)] < coarsenBound_);
#endif
    }

    /*! \brief Initializes the adaption indicator class*/
    void init()
    {
        refineBound_ = 0.;
        coarsenBound_ = 0.;
    };

    /*! \brief Function for changing the indicatorVector values for refinement
     *
     *  \param idx Index of cell which may be refined
     */
    void setIndicatorToRefine(int idx)
    {
        indicatorVector_[idx] = coarsenBound_+1;
    }

    /*! \brief Function for changing the indicatorVector values for coarsening
     *
     *  \param idx Index of cell which may be coarsen
     */
    void setIndicatorToCoarse(int idx)
    {
        indicatorVector_[idx] = refineBound_ - 1;
    }

    /*! \brief Constructs a GridAdaptionIndicator instance
     *
     *  This standard indicator is based on the saturation gradient.
     *  It checks the local gradient compared to the maximum global gradient.
     *  The indicator is compared locally to a refinement/coarsening threshold to decide whether
     *  a cell should be marked for refinement or coarsening or should not be adapted.
     *
     * \param problem The problem object
     */
    GridAdaptionIndicator2PLocal (Problem& problem):
        problem_(problem)
    {
        refinetol_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, GridAdapt, RefineTolerance);
        coarsentol_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, GridAdapt, CoarsenTolerance);
        refineAtDirichletBC_ = GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, RefineAtDirichletBC);
        refineAtFluxBC_ = GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, RefineAtFluxBC);
        refineAtSource_ = GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, RefineAtSource);
    }

private:
    Problem& problem_;
    Scalar refinetol_;
    Scalar coarsentol_;
    Scalar refineBound_;
    Scalar coarsenBound_;
    ScalarSolutionType indicatorVector_;
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
    bool refineAtDirichletBC_;
    bool refineAtFluxBC_;
    bool refineAtSource_;
};
}

#endif
