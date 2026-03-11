// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Wrapper around the current FVGridVariables to fulfill the layout
 *        of the new grid variables to test grid variables-based assembly.
 */
#ifndef DUMUX_COMPRESSIBLE_ONEP_TEST_GRID_VARIABLES_HH
#define DUMUX_COMPRESSIBLE_ONEP_TEST_GRID_VARIABLES_HH

#include <memory>

#include <dumux/experimental/common/variables.hh>

namespace Dumux::OnePCompressibleTest {

template<class GG, class BaseGridVariables, class SolutionVector>
class TestGridVariables
: public BaseGridVariables
, public Dumux::Experimental::Variables<SolutionVector>
{
    using ExperimentalBase = Dumux::Experimental::Variables<SolutionVector>;

public:
    // export some types to avoid ambiguity
    using GridGeometry = GG;
    using Scalar = typename BaseGridVariables::Scalar;

    template<class Problem>
    TestGridVariables(std::shared_ptr<Problem> problem,
                      std::shared_ptr<const GridGeometry> gridGeometry,
                      const SolutionVector& x)
    : BaseGridVariables(problem, gridGeometry)
    , ExperimentalBase(x)
    , gridGeometry_(gridGeometry)
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
    { return *gridGeometry_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
};

} // end namespace Dumux::OnePCompressibleTest

#endif
