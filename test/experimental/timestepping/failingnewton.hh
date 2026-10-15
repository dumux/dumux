// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Newton's method that reports a failure on a given call for the time stepping tests
 */
#ifndef DUMUX_TEST_TIMESTEPPING_FAILING_NEWTON_HH
#define DUMUX_TEST_TIMESTEPPING_FAILING_NEWTON_HH

#include <cstddef>

#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux::Test {

//! Newton's method that reports a failure on the given call of apply, after modifying the variables
template<class Assembler, class LinearSolver>
class FailingNewton : public NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;
public:
    using ParentType::ParentType;
    using typename ParentType::Variables;

    void failOnCall(const std::size_t call)
    { failOnCall_ = call; }

    bool apply(Variables& vars) override
    {
        const bool converged = ParentType::apply(vars);
        return (++calls_ == failOnCall_) ? false : converged;
    }

private:
    std::size_t calls_ = 0;
    std::size_t failOnCall_ = 0;
};

} // end namespace Dumux::Test

#endif
