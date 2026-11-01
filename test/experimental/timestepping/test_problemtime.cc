//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Test how the assemblers of the time integration methods pass the time to problems
 *
 * With CHECK_TIME_SETTER_REJECTION, the time is set on a problem with a time setter
 * other than `void setTime(Scalar t) const`, which must not compile.
 */
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/experimental/timestepping/implicitsolve.hh>

namespace Dumux::Test {

struct TimeIndependentProblem
{
    double value() const { return 1.0; }
};

struct TimeDependentProblem
{
    void setTime(double t) const { time_ = t; }
    double time() const { return time_; }
private:
    mutable double time_ = 0.0;
};

struct TimeDependentProblemWithUpdateTime : public TimeDependentProblem
{
    void updateTime(double t) { setTime(t); }
};

struct NonConstSetTimeProblem
{
    void setTime(double t) { time_ = t; }
private:
    double time_ = 0.0;
};

struct UpdateTimeProblem
{
    void updateTime(double t) { time_ = t; }
private:
    double time_ = 0.0;
};

struct SetTimeWithStepIndexProblem
{
    void setTime(double t, int) const { time_ = t; }
private:
    mutable double time_ = 0.0;
};

struct TimeLoopProblem
{
    void setTimeLoop(const double* time) { time_ = time; }
private:
    const double* time_ = nullptr;
};

struct TemplateSetTimeProblem
{
    template<class TimeLevel>
    void setTime(const TimeLevel&) {}
};

struct DerivedUpdateTimeProblem : public UpdateTimeProblem {};

} // end namespace Dumux::Test

int main()
{
    using namespace Dumux;
    using namespace Dumux::Experimental::Detail::TimeStepping;

    static_assert(!hasTimeSetter<Test::TimeIndependentProblem>());
    static_assert(hasTimeSetter<Test::TimeDependentProblem>());
    static_assert(hasTimeSetter<Test::TimeDependentProblemWithUpdateTime>());
    static_assert(hasTimeSetter<Test::NonConstSetTimeProblem>());
    static_assert(hasTimeSetter<Test::UpdateTimeProblem>());
    static_assert(hasTimeSetter<Test::SetTimeWithStepIndexProblem>());
    static_assert(hasTimeSetter<Test::TimeLoopProblem>());
    static_assert(hasTimeSetter<Test::TemplateSetTimeProblem>());
    static_assert(hasTimeSetter<Test::DerivedUpdateTimeProblem>());

    const Test::TimeIndependentProblem timeIndependentProblem;
    setProblemTime(timeIndependentProblem, 2.5);

    const Test::TimeDependentProblem timeDependentProblem;
    setProblemTime(timeDependentProblem, 2.5);
    if (timeDependentProblem.time() != 2.5)
        DUNE_THROW(Dune::Exception, "The time was not set: " << timeDependentProblem.time());

    const Test::TimeDependentProblemWithUpdateTime problemWithUpdateTime;
    setProblemTime(problemWithUpdateTime, 3.5);
    if (problemWithUpdateTime.time() != 3.5)
        DUNE_THROW(Dune::Exception, "The time was not set: " << problemWithUpdateTime.time());

#if CHECK_TIME_SETTER_REJECTION
    const Test::UpdateTimeProblem updateTimeProblem;
    setProblemTime(updateTimeProblem, 2.5);
#endif

    std::cout << "The time is passed to time-dependent problems" << std::endl;
    return 0;
}
