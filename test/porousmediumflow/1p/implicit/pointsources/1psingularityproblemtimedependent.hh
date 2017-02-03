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
/*!
 * \file
 *
 * \brief A test problem for the one-phase model:
 * Water is injected in one single point in the middle of the domain.
 */
#ifndef DUMUX_1P_SINGULARITY_TIME_DEP_PROBLEM_HH
#define DUMUX_1P_SINGULARITY_TIME_DEP_PROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "1psingularityproblem.hh"

namespace Dumux
{
template <class TypeTag>
class OnePSingularityProblemTimeDependent;

namespace Properties
{
NEW_TYPE_TAG(OnePSingularityProblemTimeDependent, INHERITS_FROM(OnePSingularityCCProblem));

// Set the problem property
SET_TYPE_PROP(OnePSingularityProblemTimeDependent, Problem, OnePSingularityProblemTimeDependent<TypeTag>);

// point source
SET_TYPE_PROP(OnePSingularityProblemTimeDependent, PointSource, SolDependentPointSource<TypeTag>);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 * Water is injected in a single point in the middle of the domain.
 *
 * The domain is box shaped. All sides have Dirichlet boundary conditions.
 *
 * In the middle of the domain, water is injected simulating a small diameter well.
 * The injection rate is time dependent.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p_pointsources</tt> or
 * <tt>./test_cc1p_pointsources</tt>
 */
template <class TypeTag>
class OnePSingularityProblemTimeDependent : public OnePSingularityProblem<TypeTag>
{
    using ParentType = OnePSingularityProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    OnePSingularityProblemTimeDependent(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Applies a vector of point sources. The point sources
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
        // we use t+1 because this is an implicit Euler scheme
        auto function = [](const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolume& scv)
        { return PrimaryVariables(problem.timeManager().time() + problem.timeManager().timeStepSize()); };

        pointSources.push_back(PointSource({0.0, 0.0}, function));
    }
};

} //end namespace

#endif
