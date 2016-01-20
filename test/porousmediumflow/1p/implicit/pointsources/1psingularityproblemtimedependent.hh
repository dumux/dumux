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

#include <dumux/implicit/1p/1pmodel.hh>
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
SET_TYPE_PROP(OnePSingularityProblemTimeDependent, Problem, Dumux::OnePSingularityProblemTimeDependent<TypeTag>);

// point source
SET_TYPE_PROP(OnePSingularityProblemTimeDependent, PointSource, Dumux::TimeDependentPointSource<TypeTag>);
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
    typedef OnePSingularityProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;

    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

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
     * \param pointSources A vector of Dumux::PointSource s that contain
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
        auto function = [](const TimeManager& timeManager, const GlobalPosition& pos)
        { return PrimaryVariables(timeManager.time()); };

        pointSources.push_back(PointSource({0, 0}, function));
    }
};

} //end namespace

#endif
