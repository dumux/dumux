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

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_POINTSOURCEDATA_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_POINTSOURCEDATA_HH

#include <vector>
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

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;
    template<std::size_t id> using SolutionVector = GetPropType<SubDomainTypeTag<id>, Properties::SolutionVector>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    template<std::size_t id>
    static constexpr bool isBox()
    { return GridGeometry<id>::discMethod == DiscretizationMethod::box; }

public:
    void addBulkInterpolation(const ShapeValues& shapeValues,
                              const std::vector<GridIndex<bulkIdx>>& cornerIndices,
                              GridIndex<bulkIdx> eIdx)
    {
        static_assert(isBox<bulkIdx>(), "This interface is only available for the box method.");
        bulkElementIdx_ = eIdx;
        bulkCornerIndices_ = cornerIndices;
        bulkShapeValues_ = shapeValues;
    }

    void addLowDimInterpolation(const ShapeValues& shapeValues,
                                const std::vector<GridIndex<lowDimIdx>>& cornerIndices,
                                GridIndex<lowDimIdx> eIdx)
    {
        static_assert(isBox<lowDimIdx>(), "This interface is only available for the box method.");
        lowDimElementIdx_ = eIdx;
        lowDimCornerIndices_ = cornerIndices;
        lowDimShapeValues_ = shapeValues;
    }

    void addBulkInterpolation(GridIndex<bulkIdx> eIdx)
    {
        static_assert(!isBox<bulkIdx>(), "This interface is not available for the box method.");
        bulkElementIdx_ = eIdx;
    }

    void addLowDimInterpolation(GridIndex<lowDimIdx> eIdx)
    {
        static_assert(!isBox<lowDimIdx>(), "This interface is not available for the box method.");
        lowDimElementIdx_ = eIdx;
    }

    PrimaryVariables<bulkIdx> interpolateBulk(const SolutionVector<bulkIdx>& sol) const
    {
        PrimaryVariables<bulkIdx> bulkPriVars(0.0);
        if (isBox<bulkIdx>())
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

    PrimaryVariables<lowDimIdx> interpolateLowDim(const SolutionVector<lowDimIdx>& sol) const
    {
        PrimaryVariables<lowDimIdx> lowDimPriVars(0.0);
        if (isBox<lowDimIdx>())
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
    GridIndex<bulkIdx> bulkElementIdx_;
    GridIndex<lowDimIdx> lowDimElementIdx_;
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

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;
    template<std::size_t id> using SolutionVector = GetPropType<SubDomainTypeTag<id>, Properties::SolutionVector>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    template<std::size_t id>
    static constexpr bool isBox()
    { return GridGeometry<id>::discMethod == DiscretizationMethod::box; }

public:
    PointSourceDataCircleAverage() : enableBulkCircleInterpolation_(false) {}

    PrimaryVariables<bulkIdx> interpolateBulk(const SolutionVector<bulkIdx>& sol) const
    {
        // bulk interpolation on the circle is only enabled for source in the
        // lower dimensional domain if we use a circle distributed bulk sources
        // (see coupling manager circle sources)
        PrimaryVariables<bulkIdx> bulkPriVars(sol[0]);
        bulkPriVars = 0.0;
        if (enableBulkCircleInterpolation_)
        {
            // compute the average of the bulk privars over the circle around the integration point
            // this computes $\bar{p} = \frac{1}{2\pi R} int_0^2*\pi p R \text{d}\theta.
            if (isBox<bulkIdx>())
            {
                assert(circleCornerIndices_.size() == circleShapeValues_.size());

                Scalar weightSum = 0.0;
                for (std::size_t j = 0; j < circleStencil_.size(); ++j)
                {
                    PrimaryVariables<bulkIdx> priVars(0.0);
                    const auto& cornerIndices = *(circleCornerIndices_[j]);
                    const auto& shapeValues = circleShapeValues_[j];
                    for (int i = 0; i < cornerIndices.size(); ++i)
                    {
                        const auto& localSol = sol[cornerIndices[i]];
                        const auto& shapeValue = shapeValues[i];
                        for (int priVarIdx = 0; priVarIdx < PrimaryVariables<bulkIdx>::size(); ++priVarIdx)
                            priVars[priVarIdx] += localSol[priVarIdx]*shapeValue;
                    }
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

    void addCircleInterpolation(const std::vector<const std::vector<GridIndex<bulkIdx>>*>& circleCornerIndices,
                                const std::vector<ShapeValues>& circleShapeValues,
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

    const std::vector<GridIndex<bulkIdx>>& circleStencil() const
    { return circleStencil_; }

private:
    std::vector<const std::vector<GridIndex<bulkIdx>>*> circleCornerIndices_;
    std::vector<ShapeValues> circleShapeValues_;
    std::vector<Scalar> circleIpWeight_;
    std::vector<GridIndex<bulkIdx>> circleStencil_;
    bool enableBulkCircleInterpolation_;
};

} // end namespace Dumux

#endif
