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
 * \ingroup SequentialTwoPModel
 * \brief  Class defining a standard, saturation dependent indicator for grid adaption.
 */
#ifndef DUMUX_GRIDADAPTIONINDICATOR2PLOCAL_HH
#define DUMUX_GRIDADAPTIONINDICATOR2PLOCAL_HH

#include <dumux/porousmediumflow/sequential/impetproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief  Class defining a standard, saturation dependent indicator for grid adaption.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridAdaptionIndicator2PLocal
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GridView::Traits::template Codim<0>::Entity;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ScalarSolutionType = typename SolutionTypes::ScalarSolution;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

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
    /*!
     * \brief Calculates the indicator used for refinement/coarsening for each grid cell.
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

        // 1) calculate Indicator -> min, maxvalues
        // loop over all leaf-elements
        for (const auto& element : elements(problem_.gridView()))
        {
            // calculate minimum and maximum saturation
            // index of the current leaf-elements
            int globalIdxI = problem_.variables().index(element);

            if (refineAtSource_)
            {
            PrimaryVariables source(0.0);
            problem_.source(source, element);
            for (int i = 0; i < 2; i++)
            {
                using std::abs;
                if (abs(source[i]) > 1e-10)
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

            using std::min;
            globalMin = min(satI, globalMin);
            using std::max;
            globalMax = max(satI, globalMax);

            // calculate refinement indicator in all cells
            for (const auto& intersection : intersections(problem_.gridView(), element))
            {
                if (indicatorVector_[globalIdxI] == 10)
                {
                    break;
                }

                // exit, if it is not a neighbor
                if (intersection.boundary())
                {
                    BoundaryTypes bcTypes;
                    problem_.boundaryTypes(bcTypes, intersection);

                    for (int i = 0; i < 2; i++)
                    {
                        if (bcTypes.isNeumann(i))
                        {
                            PrimaryVariables flux(0.0);
                            problem_.neumann(flux, intersection);

                            bool fluxBound = false;
                            for (int j = 0; j < 2; j++)
                            {
                            using std::abs;
                            if (abs(flux[j]) > 1e-10)
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
                    if (element.level() > outside.level()
                            || (element.level() == outside.level() && globalIdxI < globalIdxJ))
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

                        using std::abs;
                        Scalar localdelta = abs(satI - satJ);
                        using std::max;
                        indicatorVector_[globalIdxI][0] = max(indicatorVector_[globalIdxI][0], localdelta);
                        indicatorVector_[globalIdxJ][0] = max(indicatorVector_[globalIdxJ][0], localdelta);

                        maxLocalDelta = max(maxLocalDelta, localdelta);
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

    /*!
     * \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool refine(const Element& element)
    {
        return (indicatorVector_[problem_.elementMapper().index(element)] > refineBound_);
    }

    /*!
     * \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool coarsen(const Element& element)
    {
        return (indicatorVector_[problem_.elementMapper().index(element)] < coarsenBound_);
    }

    /*!
     * \brief Initializes the adaption indicator class
     */
    void init()
    {
        refineBound_ = 0.;
        coarsenBound_ = 0.;
    };

    /*!
     * \brief Function for changing the indicatorVector values for refinement
     *
     *  \param idx Index of cell which may be refined
     */
    void setIndicatorToRefine(int idx)
    {
        indicatorVector_[idx] = coarsenBound_+1;
    }

    /*!
     * \brief Function for changing the indicatorVector values for coarsening
     *
     *  \param idx Index of cell which may be coarsen
     */
    void setIndicatorToCoarse(int idx)
    {
        indicatorVector_[idx] = refineBound_ - 1;
    }

    /*!
     * \brief Constructs a GridAdaptionIndicator instance
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
        refinetol_ = getParam<Scalar>("GridAdapt.RefineTolerance");
        coarsentol_ = getParam<Scalar>("GridAdapt.CoarsenTolerance");
        refineAtDirichletBC_ = getParam<bool>("GridAdapt.RefineAtDirichletBC");
        refineAtFluxBC_ = getParam<bool>("GridAdapt.RefineAtFluxBC");
        refineAtSource_ = getParam<bool>("GridAdapt.RefineAtSource");
    }

private:
    Problem& problem_;
    Scalar refinetol_;
    Scalar coarsentol_;
    Scalar refineBound_;
    Scalar coarsenBound_;
    ScalarSolutionType indicatorVector_;
    static const int saturationType_ = getPropValue<TypeTag, Properties::SaturationFormulation>();
    bool refineAtDirichletBC_; // switch for refinement at Dirichlet BC's
    bool refineAtFluxBC_; // switch for refinement at Neumann BC's
    bool refineAtSource_; // switch for refinement at sources
};
} // end namespace Dumux

#endif
