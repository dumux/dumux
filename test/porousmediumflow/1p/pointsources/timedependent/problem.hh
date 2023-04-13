// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief A test problem for the one-phase model:
 * Water is injected in one single point in the middle of the domain.
 */

#ifndef DUMUX_1P_SINGULARITY_TIME_DEP_PROBLEM_HH
#define DUMUX_1P_SINGULARITY_TIME_DEP_PROBLEM_HH

#include "../timeindependent/problem.hh"

namespace Dumux {
/*!
 * \ingroup OnePTests
 * \brief  Test problem for the one-phase model:
 * Water is injected in a single point in the middle of the domain.
 *
 * The domain is box shaped. All sides have Dirichlet boundary conditions.
 *
 * In the middle of the domain, water is injected simulating a small diameter well.
 * The injection rate is time dependent.
 */
template <class TypeTag>
class OnePSingularityProblemTimeDependent : public OnePSingularityProblem<TypeTag>
{
    using ParentType = OnePSingularityProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    OnePSingularityProblemTimeDependent(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Applies a vector of point sources which
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in untis
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        // inject <time> kg/s water at position (0, 0), where <time> is the current simulation time
        auto function = [](const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolume& scv)
        { return PrimaryVariables(problem.getTime()); };

        pointSources.push_back(PointSource({0.0, 0.0}, function));
    }

    //! Set the current time at which we evaluate the source
    void setTime(Scalar time)
    { time_ = time; }

    //! Set the current time at which we evaluate the source
    Scalar getTime() const
    { return time_; }

private:
    Scalar time_ = 0.0;
};

} // end namespace Dumux

#endif
