// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Data associated with a point source
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_POINTSOURCEDATA_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_POINTSOURCEDATA_HH

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>
#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief A point source data class used for integration in multidimensional models
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
    static constexpr bool isCVFE()
    { return DiscretizationMethods::isCVFE<typename GridGeometry<id>::DiscretizationMethod>; }

public:
    /*!
     * \brief Add the interpolation data in the bulk domain
     * \param shapeValues the values of the basis functions of the local dofs at the point source
     * \param dofIndices the global indices of the local dofs, in the order of the basis functions
     * \param eIdx the index of the element containing the point source
     */
    void addBulkInterpolation(const ShapeValues& shapeValues,
                              const std::vector<GridIndex<bulkIdx>>& dofIndices,
                              GridIndex<bulkIdx> eIdx)
    {
        static_assert(isCVFE<bulkIdx>(), "This interface is only available for control-volume finite element methods.");
        bulkElementIdx_ = eIdx;
        bulkDofIndices_ = dofIndices;
        bulkShapeValues_ = shapeValues;
    }

    /*!
     * \brief Add the interpolation data in the lower-dimensional domain
     * \param shapeValues the values of the basis functions of the local dofs at the point source
     * \param dofIndices the global indices of the local dofs, in the order of the basis functions
     * \param eIdx the index of the element containing the point source
     */
    void addLowDimInterpolation(const ShapeValues& shapeValues,
                                const std::vector<GridIndex<lowDimIdx>>& dofIndices,
                                GridIndex<lowDimIdx> eIdx)
    {
        static_assert(isCVFE<lowDimIdx>(), "This interface is only available for control-volume finite element methods.");
        lowDimElementIdx_ = eIdx;
        lowDimDofIndices_ = dofIndices;
        lowDimShapeValues_ = shapeValues;
    }

    void addBulkInterpolation(GridIndex<bulkIdx> eIdx)
    {
        static_assert(!isCVFE<bulkIdx>(), "This interface is not available for control-volume finite element methods.");
        bulkElementIdx_ = eIdx;
    }

    void addLowDimInterpolation(GridIndex<lowDimIdx> eIdx)
    {
        static_assert(!isCVFE<lowDimIdx>(), "This interface is not available for control-volume finite element methods.");
        lowDimElementIdx_ = eIdx;
    }

    PrimaryVariables<bulkIdx> interpolateBulk(const SolutionVector<bulkIdx>& sol) const
    {
        PrimaryVariables<bulkIdx> bulkPriVars(0.0);
        if (isCVFE<bulkIdx>())
        {
            for (int i = 0; i < bulkDofIndices_.size(); ++i)
                for (int priVarIdx = 0; priVarIdx < bulkPriVars.size(); ++priVarIdx)
                    bulkPriVars[priVarIdx] += sol[bulkDofIndices_[i]][priVarIdx]*bulkShapeValues_[i][0];
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
        if (isCVFE<lowDimIdx>())
        {
            for (int i = 0; i < lowDimDofIndices_.size(); ++i)
                for (int priVarIdx = 0; priVarIdx < lowDimPriVars.size(); ++priVarIdx)
                    lowDimPriVars[priVarIdx] += sol[lowDimDofIndices_[i]][priVarIdx]*lowDimShapeValues_[i][0];
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
    std::vector<GridIndex<bulkIdx>> bulkDofIndices_;
    std::vector<GridIndex<lowDimIdx>> lowDimDofIndices_;
    GridIndex<bulkIdx> bulkElementIdx_ = std::numeric_limits<GridIndex<bulkIdx>>::max();
    GridIndex<lowDimIdx> lowDimElementIdx_ = std::numeric_limits<GridIndex<lowDimIdx>>::max();
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief A point source data class used for integration in multidimensional models
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
    static constexpr bool isCVFE()
    { return DiscretizationMethods::isCVFE<typename GridGeometry<id>::DiscretizationMethod>; }

public:
    /*!
     * \brief The circle average of the bulk solution: a weight per bulk degree of freedom (the weights add
     *        up to one) and the bulk elements on the circle
     */
    struct CircleInterpolation
    {
        std::vector<GridIndex<bulkIdx>> dofs;
        std::vector<Scalar> weights;
        std::vector<GridIndex<bulkIdx>> stencil;
    };

    /*!
     * \brief The circle average from the interpolation data of the circle points (control-volume finite element methods)
     * \param circleDofIndices the bulk dofs of the element of each circle point
     * \param circleShapeValues the values of the basis functions at each circle point
     * \param circleIpWeight the weight of each circle point
     * \param circleStencil the bulk element of each circle point
     */
    static std::shared_ptr<const CircleInterpolation>
    makeCircleInterpolation(const std::vector<const std::vector<GridIndex<bulkIdx>>*>& circleDofIndices,
                            const std::vector<ShapeValues>& circleShapeValues,
                            const std::vector<Scalar>& circleIpWeight,
                            const std::vector<GridIndex<bulkIdx>>& circleStencil)
    {
        auto circle = std::make_shared<CircleInterpolation>();
        Scalar weightSum = 0.0;
        for (std::size_t j = 0; j < circleStencil.size(); ++j)
        {
            weightSum += circleIpWeight[j];
            const auto& dofIndices = *circleDofIndices[j];
            for (std::size_t i = 0; i < dofIndices.size(); ++i)
                addWeight_(circle->dofs, circle->weights, dofIndices[i], circleIpWeight[j]*circleShapeValues[j][i][0]);
        }
        finish_(*circle, circleStencil, weightSum);
        return circle;
    }

    /*!
     * \brief The circle average from the interpolation data of the circle points (cell-centered methods)
     * \param circleIpWeight the weight of each circle point
     * \param circleStencil the bulk element of each circle point
     */
    static std::shared_ptr<const CircleInterpolation>
    makeCircleInterpolation(const std::vector<Scalar>& circleIpWeight,
                            const std::vector<GridIndex<bulkIdx>>& circleStencil)
    {
        auto circle = std::make_shared<CircleInterpolation>();
        Scalar weightSum = 0.0;
        for (std::size_t j = 0; j < circleStencil.size(); ++j)
        {
            weightSum += circleIpWeight[j];
            addWeight_(circle->dofs, circle->weights, circleStencil[j], circleIpWeight[j]);
        }
        finish_(*circle, circleStencil, weightSum);
        return circle;
    }

    PrimaryVariables<bulkIdx> interpolateBulk(const SolutionVector<bulkIdx>& sol) const
    {
        // bulk interpolation on the circle is only enabled for source in the
        // lower dimensional domain if we use a circle distributed bulk sources
        // (see coupling manager circle sources)
        if (circle_)
        {
            // compute the average of the bulk privars over the circle around the integration point
            // this computes $\bar{p} = \frac{1}{2\pi R} int_0^2*\pi p R \text{d}\theta.
            PrimaryVariables<bulkIdx> bulkPriVars(0.0);
            for (std::size_t i = 0; i < circle_->dofs.size(); ++i)
                for (int priVarIdx = 0; priVarIdx < PrimaryVariables<bulkIdx>::size(); ++priVarIdx)
                    bulkPriVars[priVarIdx] += sol[circle_->dofs[i]][priVarIdx]*circle_->weights[i];
            return bulkPriVars;
        }
        else
            return ParentType::interpolateBulk(sol);
    }

    //! Use the circle average for the bulk values, shared with the other point sources on the circle
    void addCircleInterpolation(std::shared_ptr<const CircleInterpolation> circle)
    { circle_ = std::move(circle); }

    void addCircleInterpolation(const std::vector<const std::vector<GridIndex<bulkIdx>>*>& circleDofIndices,
                                const std::vector<ShapeValues>& circleShapeValues,
                                const std::vector<Scalar>& circleIpWeight,
                                const std::vector<GridIndex<bulkIdx>>& circleStencil)
    { circle_ = makeCircleInterpolation(circleDofIndices, circleShapeValues, circleIpWeight, circleStencil); }

    void addCircleInterpolation(const std::vector<Scalar>& circleIpWeight,
                                const std::vector<GridIndex<bulkIdx>>& circleStencil)
    { circle_ = makeCircleInterpolation(circleIpWeight, circleStencil); }

    //! The bulk elements on the circle, each once
    const std::vector<GridIndex<bulkIdx>>& circleStencil() const
    {
        static const std::vector<GridIndex<bulkIdx>> empty;
        return circle_ ? circle_->stencil : empty;
    }

private:
    static void addWeight_(std::vector<GridIndex<bulkIdx>>& dofs, std::vector<Scalar>& weights,
                           GridIndex<bulkIdx> dof, Scalar weight)
    {
        const auto it = std::find(dofs.begin(), dofs.end(), dof);
        if (it != dofs.end())
            weights[it - dofs.begin()] += weight;
        else
        {
            dofs.push_back(dof);
            weights.push_back(weight);
        }
    }

    static void finish_(CircleInterpolation& circle, const std::vector<GridIndex<bulkIdx>>& circleStencil, Scalar weightSum)
    {
        for (auto& w : circle.weights)
            w /= weightSum;
        circle.stencil = circleStencil;
        std::sort(circle.stencil.begin(), circle.stencil.end());
        circle.stencil.erase(std::unique(circle.stencil.begin(), circle.stencil.end()), circle.stencil.end());
    }

    std::shared_ptr<const CircleInterpolation> circle_;
};

} // end namespace Dumux

#endif
