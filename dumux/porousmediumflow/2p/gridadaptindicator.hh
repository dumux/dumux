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
 * \ingroup TwoPModel
 * \brief Class defining a standard, saturation dependent indicator for grid adaptation.
 */

#ifndef DUMUX_TWOP_ADAPTION_INDICATOR_HH
#define DUMUX_TWOP_ADAPTION_INDICATOR_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief  Class defining a standard, saturation dependent indicator for grid adaptation.
 */
template<class TypeTag>
class TwoPGridAdaptIndicator
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    enum { saturationIdx = Indices::saturationIdx };

public:
    /*!
     * \brief The Constructor
     *
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     *
     *  Note: refineBound_, coarsenBound_ & maxSaturationDelta_ are chosen
     *        in a way such that the indicator returns false for all elements
     *        before having been calculated.
     */
    TwoPGridAdaptIndicator(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : gridGeometry_(gridGeometry)
    , refineBound_(std::numeric_limits<Scalar>::max())
    , coarsenBound_(std::numeric_limits<Scalar>::lowest())
    , maxSaturationDelta_(gridGeometry_->gridView().size(0), 0.0)
    , minLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MinLevel", 0))
    , maxLevel_(getParamFromGroup<std::size_t>(paramGroup, "Adaptive.MaxLevel", 0))
    {}

    /*!
     * \brief Function to set the minimum allowed level.
     *
     * \param minLevel The minimum level
     */
    void setMinLevel(std::size_t minLevel)
    {
        minLevel_ = minLevel;
    }

    /*!
     * \brief Function to set the maximum allowed level.
     *
     *\param maxLevel The maximum level
     */
    void setMaxLevel(std::size_t maxLevel)
    {
        maxLevel_ = maxLevel;
    }

    /*!
     * \brief Function to set the minumum/maximum allowed levels.
     *
     * \param minLevel The minimum level
     * \param maxLevel The maximum level
     */
    void setLevels(std::size_t minLevel, std::size_t maxLevel)
    {
        minLevel_ = minLevel;
        maxLevel_ = maxLevel;
    }

    /*!
     * \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     * \param sol The solution vector
     * \param refineTol The refinement tolerance
     * \param coarsenTol The coarsening tolerance
     *
     *  This standard two-phase indicator is based on the saturation gradient.
     */
    void calculate(const SolutionVector& sol,
                   Scalar refineTol = 0.05,
                   Scalar coarsenTol = 0.001)
    {
        //! Reset the indicator to a state that returns false for all elements
        refineBound_ = std::numeric_limits<Scalar>::max();
        coarsenBound_ = std::numeric_limits<Scalar>::lowest();
        maxSaturationDelta_.assign(gridGeometry_->gridView().size(0), 0.0);

        //! maxLevel_ must be higher than minLevel_ to allow for refinement
        if (minLevel_ >= maxLevel_)
            return;

        //! Check for inadmissible tolerance combination
        if (coarsenTol > refineTol)
            DUNE_THROW(Dune::InvalidStateException, "Refine tolerance must be higher than coarsen tolerance");

        //! Variables to hold the max/mon saturation values on the leaf
        Scalar globalMax = std::numeric_limits<Scalar>::lowest();
        Scalar globalMin = std::numeric_limits<Scalar>::max();

        //! Calculate minimum and maximum saturation
        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            //! Index of the current leaf-element
            const auto globalIdxI = gridGeometry_->elementMapper().index(element);

            //! Obtain the saturation at the center of the element
            const auto geometry = element.geometry();
            const auto elemSol = elementSolution(element, sol, *gridGeometry_);
            const Scalar satI = evalSolution(element, geometry, *gridGeometry_, elemSol, geometry.center())[saturationIdx];

            //! Maybe update the global minimum/maximum
            using std::min;
            using std::max;
            globalMin = min(satI, globalMin);
            globalMax = max(satI, globalMax);

            //! Calculate maximum delta in saturation for this cell
            for (const auto& intersection : intersections(gridGeometry_->gridView(), element))
            {
                //! Only consider internal intersections
                if (intersection.neighbor())
                {
                    //! Access neighbor
                    const auto outside = intersection.outside();
                    const auto globalIdxJ = gridGeometry_->elementMapper().index(outside);

                    //! Visit intersection only once
                    if (element.level() > outside.level() || (element.level() == outside.level() && globalIdxI < globalIdxJ))
                    {
                        //! Obtain saturation in the neighbor
                        const auto outsideGeometry = outside.geometry();
                        const auto elemSolJ = elementSolution(outside, sol, *gridGeometry_);
                        const Scalar satJ = evalSolution(outside, outsideGeometry, *gridGeometry_, elemSolJ, outsideGeometry.center())[saturationIdx];

                        using std::abs;
                        Scalar localdelta = abs(satI - satJ);
                        maxSaturationDelta_[globalIdxI] = max(maxSaturationDelta_[globalIdxI], localdelta);
                        maxSaturationDelta_[globalIdxJ] = max(maxSaturationDelta_[globalIdxJ], localdelta);
                    }
                }
            }
        }

        //! Compute the maximum delta in saturation
        const auto globalDelta = globalMax - globalMin;

        //! Compute the refinement/coarsening bounds
        refineBound_ = refineTol*globalDelta;
        coarsenBound_ = coarsenTol*globalDelta;

// TODO: fix adaptive simulations in parallel
//#if HAVE_MPI
//    // communicate updated values
//    using DataHandle = VectorExchange<ElementMapper, ScalarSolutionType>;
//    DataHandle dataHandle(problem_.elementMapper(), maxSaturationDelta_);
//    problem_.gridView().template communicate<DataHandle>(dataHandle,
//                                                         Dune::InteriorBorder_All_Interface,
//                                                         Dune::ForwardCommunication);
//
//    using std::max;
//    refineBound_ = problem_.gridView().comm().max(refineBound_);
//    coarsenBound_ = problem_.gridView().comm().max(coarsenBound_);
//
//#endif

        //! check if neighbors have to be refined too
        for (const auto& element : elements(gridGeometry_->gridView(), Dune::Partitions::interior))
            if (this->operator()(element) > 0)
                checkNeighborsRefine_(element);
    }

    /*!
     * \brief function call operator to return mark
     *
     * \return  1 if an element should be refined
     *         -1 if an element should be coarsened
     *          0 otherwise
     *
     * \param element A grid element
     */
    int operator() (const Element& element) const
    {
        if (element.hasFather()
            && maxSaturationDelta_[gridGeometry_->elementMapper().index(element)] < coarsenBound_)
        {
            return -1;
        }
        else if (element.level() < maxLevel_
                 && maxSaturationDelta_[gridGeometry_->elementMapper().index(element)] > refineBound_)
        {
            return 1;
        }
        else
            return 0;
    }

private:
    /*!
     * \brief Method ensuring the refinement ratio of 2:1
     *
     * For any given element, a loop over the neighbors checks if the
     * entities refinement would require that any of the neighbors has
     * to be refined, too. This is done recursively over all levels of the grid.
     *
     * \param element Element of interest that is to be refined
     * \param level level of the refined element: it is at least 1
     * \return true if everything was successful
     */
    bool checkNeighborsRefine_(const Element &element, std::size_t level = 1)
    {
        for(const auto& intersection : intersections(gridGeometry_->gridView(), element))
        {
            if(!intersection.neighbor())
                continue;

            // obtain outside element
            const auto outside = intersection.outside();

            // only mark non-ghost elements
            if (outside.partitionType() == Dune::GhostEntity)
                continue;

            if (outside.level() < maxLevel_ && outside.level() < element.level())
            {
                // ensure refinement for outside element
                maxSaturationDelta_[gridGeometry_->elementMapper().index(outside)] = std::numeric_limits<Scalar>::max();
                if(level < maxLevel_)
                    checkNeighborsRefine_(outside, ++level);
            }
        }

        return true;
    }

    std::shared_ptr<const GridGeometry> gridGeometry_;

    Scalar refineBound_;
    Scalar coarsenBound_;
    std::vector< Scalar > maxSaturationDelta_;
    std::size_t minLevel_;
    std::size_t maxLevel_;
};

} // end namespace Dumux

#endif
