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
 * \ingroup EmbeddedCoupling
 * \brief Data associated with a point source
 */

#ifndef DUMUX_POINTSOURCEDATA_HH
#define DUMUX_POINTSOURCEDATA_HH

#include <vector>
#include <unordered_map>
#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief A point source data class used for integration in multidimension models
 * \note The point source and related data are connected via an identifier (id)
 */
template<class MDTraits>
class PointSourceData
{
    using Scalar = typename MDTraits::Scalar;
    using ShapeValues = typename std::vector<Dune::FieldVector<Scalar, 1> >;

    // obtain the type tags of the sub problems
    using BulkTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using LowDimTypeTag = typename MDTraits::template SubDomain<1>::TypeTag;

    using BulkPrimaryVariables = GetPropType<BulkTypeTag, Properties::PrimaryVariables>;
    using LowDimPrimaryVariables = GetPropType<LowDimTypeTag, Properties::PrimaryVariables>;

    using BulkSolutionVector = GetPropType<BulkTypeTag, Properties::SolutionVector>;
    using LowDimSolutionVector = GetPropType<LowDimTypeTag, Properties::SolutionVector>;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    enum {
        bulkIsBox = GetPropType<BulkTypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box,
        lowDimIsBox = GetPropType<LowDimTypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box
    };

public:
    void addBulkInterpolation(const ShapeValues& shapeValues,
                              const std::vector<GridIndex<bulkIdx>>& cornerIndices,
                              GridIndex<bulkIdx> eIdx)
    {
        bulkElementIdx_ = eIdx;
        bulkCornerIndices_ = cornerIndices;
        bulkShapeValues_ = shapeValues;
    }

    void addLowDimInterpolation(const ShapeValues& shapeValues,
                                const std::vector<GridIndex<lowDimIdx>>& cornerIndices,
                                GridIndex<lowDimIdx> eIdx)
    {
        lowDimElementIdx_ = eIdx;
        lowDimCornerIndices_ = cornerIndices;
        lowDimShapeValues_ = shapeValues;
    }

    void addBulkInterpolation(GridIndex<bulkIdx> eIdx)
    {
        assert(!bulkIsBox);
        bulkElementIdx_ = eIdx;
    }

    void addLowDimInterpolation(GridIndex<lowDimIdx> eIdx)
    {
        assert(!lowDimIsBox);
        lowDimElementIdx_ = eIdx;
    }

    BulkPrimaryVariables interpolateBulk(const BulkSolutionVector& sol)
    {
        BulkPrimaryVariables bulkPriVars(0.0);
        if (bulkIsBox)
        {
            for (int i = 0; i < bulkCornerIndices_.size(); ++i)
                for (int priVarIdx = 0; priVarIdx < bulkPriVars.size(); ++priVarIdx)
                    bulkPriVars[priVarIdx] += sol[bulkCornerIndices_[i]][priVarIdx]*bulkShapeValues_[i];
        }
        else
        {
            bulkPriVars = sol[bulkElementIdx()];
        }
        return bulkPriVars;
    }

    LowDimPrimaryVariables interpolateLowDim(const LowDimSolutionVector& sol)
    {
        LowDimPrimaryVariables lowDimPriVars(0.0);
        if (lowDimIsBox)
        {
            for (int i = 0; i < lowDimCornerIndices_.size(); ++i)
                for (int priVarIdx = 0; priVarIdx < lowDimPriVars.size(); ++priVarIdx)
                    lowDimPriVars[priVarIdx] += sol[lowDimCornerIndices_[i]][priVarIdx]*lowDimShapeValues_[i];
        }
        else
        {
            lowDimPriVars = sol[lowDimElementIdx()];
        }
        return lowDimPriVars;
    }

    GridIndex<lowDimIdx> lowDimElementIdx() const
    { return lowDimElementIdx_; }

    GridIndex<bulkIdx> bulkElementIdx() const
    { return bulkElementIdx_; }

private:
    ShapeValues bulkShapeValues_, lowDimShapeValues_;
    std::vector<GridIndex<bulkIdx>> bulkCornerIndices_;
    std::vector<GridIndex<lowDimIdx>> lowDimCornerIndices_;
    GridIndex<lowDimIdx> lowDimElementIdx_;
    GridIndex<bulkIdx> bulkElementIdx_;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief A point source data class used for integration in multidimension models
 * \note The point source and related data are connected via an identifier (id)
 * When explicitly computing the circle average, i.e. the pressure for
 * the source term is computed as an integral over the circle describing
 * the surface of the one-dimensional tube. This exact determination of the bulk pressure
 * is necessary for convergence studies.
 */
template<class MDTraits>
class PointSourceDataCircleAverage : public PointSourceData<MDTraits>
{
    using ParentType = PointSourceData<MDTraits>;
    using Scalar = typename MDTraits::Scalar;
    using ShapeValues = typename std::vector<Dune::FieldVector<Scalar, 1> >;

    // obtain the type tags of the sub problems
    using BulkTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using LowDimTypeTag = typename MDTraits::template SubDomain<1>::TypeTag;

    using BulkPrimaryVariables = GetPropType<BulkTypeTag, Properties::PrimaryVariables>;
    using LowDimPrimaryVariables = GetPropType<LowDimTypeTag, Properties::PrimaryVariables>;

    using BulkSolutionVector = GetPropType<BulkTypeTag, Properties::SolutionVector>;
    using LowDimSolutionVector = GetPropType<LowDimTypeTag, Properties::SolutionVector>;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    enum {
        bulkIsBox = GetPropType<BulkTypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box,
        lowDimIsBox = GetPropType<LowDimTypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box
    };

public:
    PointSourceDataCircleAverage() : enableBulkCircleInterpolation_(false) {}

    BulkPrimaryVariables interpolateBulk(const BulkSolutionVector& sol)
    {
        // bulk interpolation on the circle is only enabled for source in the
        // lower dimensional domain if we use a circle distributed bulk sources
        // (see coupling manager circle sources)
        BulkPrimaryVariables bulkPriVars(sol[0]);
        bulkPriVars = 0.0;
        if (enableBulkCircleInterpolation_)
        {
            // compute the average of the bulk privars over the circle around the integration point
            // this computes $\bar{p} = \frac{1}{2\pi R} int_0^2*\pi p R \text{d}\theta.
            if (bulkIsBox)
            {
                assert(circleCornerIndices_.size() == circleShapeValues_.size());

                Scalar weightSum = 0.0;
                for (int j = 0; j < circleStencil_.size(); ++j)
                {
                    BulkPrimaryVariables priVars(0.0);
                    const auto& cornerIndices = circleCornerIndices_[circleStencil_[j]];
                    const auto& shapeValues = circleShapeValues_[circleStencil_[j]];
                    for (int i = 0; i < cornerIndices.size(); ++i)
                        for (int priVarIdx = 0; priVarIdx < priVars.size(); ++priVarIdx)
                            priVars[priVarIdx] += sol[cornerIndices[i]][priVarIdx]*shapeValues[i];
                    // multiply with weight and add
                    priVars *= circleIpWeight_[j];
                    weightSum += circleIpWeight_[j];
                    bulkPriVars += priVars;
                }
                bulkPriVars /= weightSum;
            }
            else
            {
                Scalar weightSum = 0.0;
                for (int j = 0; j < circleStencil_.size(); ++j)
                {
                    for (int priVarIdx = 0; priVarIdx < bulkPriVars.size(); ++priVarIdx)
                        bulkPriVars[priVarIdx] += sol[circleStencil_[j]][priVarIdx]*circleIpWeight_[j];

                    weightSum += circleIpWeight_[j];
                }
                bulkPriVars /= weightSum;
            }
            return bulkPriVars;
        }
        else
        {
            return ParentType::interpolateBulk(sol);
        }
    }

    void addCircleInterpolation(const std::unordered_map<GridIndex<bulkIdx>, std::vector<GridIndex<bulkIdx>> >& circleCornerIndices,
                                const std::unordered_map<GridIndex<bulkIdx>, ShapeValues>& circleShapeValues,
                                const std::vector<Scalar>& circleIpWeight,
                                const std::vector<GridIndex<bulkIdx>>& circleStencil)
    {
        circleCornerIndices_ = circleCornerIndices;
        circleShapeValues_ = circleShapeValues;
        circleIpWeight_ = circleIpWeight;
        circleStencil_ = circleStencil;
        enableBulkCircleInterpolation_ = true;

    }

    void addCircleInterpolation(const std::vector<Scalar>& circleIpWeight,
                                const std::vector<GridIndex<bulkIdx>>& circleStencil)
    {
        circleIpWeight_ = circleIpWeight;
        circleStencil_ = circleStencil;
        enableBulkCircleInterpolation_ = true;
    }

    const std::vector<GridIndex<bulkIdx>>& circleStencil()
    {
        return circleStencil_;
    }

private:
    std::unordered_map<GridIndex<bulkIdx>, std::vector<GridIndex<bulkIdx>> > circleCornerIndices_;
    std::unordered_map<GridIndex<bulkIdx>, ShapeValues> circleShapeValues_;
    std::vector<Scalar> circleIpWeight_;
    std::vector<GridIndex<bulkIdx>> circleStencil_;
    bool enableBulkCircleInterpolation_;
};

} // end namespace Dumux

#endif
