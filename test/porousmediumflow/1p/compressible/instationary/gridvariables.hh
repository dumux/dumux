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
 * \ingroup OnePTests
 * \brief Wrapper around the current FVGridVariables to fulfill the layout
 *        of the new grid variables to test grid variables-based assembly.
 */
#ifndef DUMUX_COMPRESSIBLE_ONEP_TEST_GRID_VARIABLES_HH
#define DUMUX_COMPRESSIBLE_ONEP_TEST_GRID_VARIABLES_HH

#include <dumux/experimental/discretization/gridvariables.hh>

namespace Dumux::OnePCompressibleTest {

template<class GG, class BaseGridVariables, class SolutionVector>
class TestGridVariables
: public BaseGridVariables
, public Dumux::Experimental::GridVariables<GG, SolutionVector>
{
    using ExperimentalBase = Dumux::Experimental::GridVariables<GG, SolutionVector>;

public:
    // export some types to avoid ambiguity
    using GridGeometry = GG;
    using Scalar = typename BaseGridVariables::Scalar;

    template<class Problem>
    TestGridVariables(std::shared_ptr<Problem> problem,
                      std::shared_ptr<const GridGeometry> gridGeometry,
                      const SolutionVector& x)
    : BaseGridVariables(problem, gridGeometry)
    , ExperimentalBase(gridGeometry, x)
    {
        BaseGridVariables::init(x);
    }

    // update to a new solution
    void update(const SolutionVector& x)
    {
        BaseGridVariables::update(x);
        ExperimentalBase::update(x);
    }

    // overload some functions to avoid ambiguity
    decltype(auto) gridGeometry() const
    { return ExperimentalBase::gridGeometry(); }
};

} // end namespace Dumux::OnePCompressibleTest

#endif
