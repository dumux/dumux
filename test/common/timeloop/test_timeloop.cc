//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// time loop test
#include <config.h>

#include <iostream>
#include <cmath>
#include <chrono>
#include <thread>
#include <limits>
#include <optional>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/timeloop.hh>

void testTimeLoops(double tStart,
                   double tEnd,
                   double dt,
                   std::optional<double> dtInLoop = {})
{
    std::cout << "Testing with dt = " << dt
              << ", tStart = " << tStart
              << ", tEnd = " << tEnd
              <<std::endl;
    const auto timeSpan = tEnd - tStart;

    //! Standard time loop
    {
        std::cout << "------- Test time loop ----------" << std::endl;
        Dumux::TimeLoop<double> timeLoop(tStart, dt, tEnd);

        if (std::abs(timeLoop.time()-tStart) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Wrong time at start!");
        if (std::abs(timeLoop.endTime()-tEnd) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Wrong end time: " << timeLoop.endTime() << " should be " << tEnd);
        if (std::abs(timeLoop.timeStepSize()-dt) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Wrong time step size!");

        const auto smallerMaxDt = dt*0.03;
        timeLoop.setMaxTimeStepSize(smallerMaxDt);
        if (std::abs(timeLoop.timeStepSize()-smallerMaxDt) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Wrong time step size!");

        timeLoop.start(); // starts the timer
        timeLoop.stop(); // stops the timer
        timeLoop.resetTimer(); // resets the timer
        timeLoop.start(); // starts the timer again

        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        if (timeLoop.timeStepIndex() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Wrong time step index!");

        timeLoop.reportTimeStep();

        while (!timeLoop.finished())
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            timeLoop.advanceTimeStep();
            timeLoop.reportTimeStep();
            timeLoop.setMaxTimeStepSize(dtInLoop.value_or(timeLoop.timeStepSize()));
            timeLoop.setTimeStepSize(dtInLoop.value_or(timeLoop.timeStepSize()));
        }

        timeLoop.finalize();

        if (std::abs(timeLoop.time()-tEnd) > 1e-15*timeSpan)
            DUNE_THROW(Dune::InvalidStateException, "Ended with wrong end time!");
    }

    // check point timeLoop
    {
        std::cout << std::endl << "------- Test check point time loop ----------" << std::endl;
        Dumux::CheckPointTimeLoop<double> timeLoop(tStart, dt, tEnd);
        timeLoop.setPeriodicCheckPoint(0.03*dt, -0.03*dt);
        if (!timeLoop.isCheckPoint())
            DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

        timeLoop.start();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();
        timeLoop.setPeriodicCheckPoint(0.03*dt);

        if (std::signbit(timeLoop.timeStepSize()))
            DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());

        if (!timeLoop.isCheckPoint())
            DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

        timeLoop.setCheckPoint(0.02*dt); // can't be set because it's in the past
        timeLoop.setCheckPoint(0.08*dt);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();

        if (std::signbit(timeLoop.timeStepSize()))
            DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());

        if (!timeLoop.isCheckPoint())
            DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

        timeLoop.setEndTime(1e9*timeSpan);
        timeLoop.setTimeStepSize(1e7*timeSpan);

        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();

        if (std::signbit(timeLoop.timeStepSize()))
            DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());

        if (!timeLoop.isCheckPoint())
            DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

        const double scaling = 1.0e6;
        timeLoop.removeAllCheckPoints();
        timeLoop.setTimeStepSize(0.5e-6*timeSpan);
        timeLoop.setPeriodicCheckPoint(scaling*timeSpan, scaling*timeSpan);

        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();

        if (std::signbit(timeLoop.timeStepSize()))
            DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());

        timeLoop.setTimeStepSize(dtInLoop.value_or(timeLoop.timeStepSize())*scaling);
        while (!timeLoop.finished())
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            timeLoop.advanceTimeStep();
            timeLoop.reportTimeStep();
            timeLoop.setTimeStepSize(std::numeric_limits<double>::max());

            if (std::abs(timeLoop.time() - 2.0e6*timeSpan) < 1e-15 && !timeLoop.isCheckPoint())
                DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

            if (timeLoop.timeStepIndex() == 10)
                timeLoop.setPeriodicCheckPoint(1.0e8*timeSpan);

            if (std::signbit(timeLoop.timeStepSize()))
                DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());
        }

        timeLoop.finalize();

        if (std::abs(timeLoop.time()-1e9*timeSpan) > 1e9*1e-15*timeSpan)
            DUNE_THROW(Dune::InvalidStateException, "Ended with wrong end time!");

    }

    // check if time loops remembers time step before check point
    {
        std::cout << std::endl << "------- Test check point time loop ----------" << std::endl;
        Dumux::CheckPointTimeLoop<double> timeLoop(tStart, dt, tEnd);
        timeLoop.setCheckPoint(1.01*dt);
        timeLoop.start();
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();
        if (!(std::abs(timeLoop.timeStepSize()-dt) < 1e-14))
            DUNE_THROW(Dune::InvalidStateException, "Time Loop reduced time step size to " << timeLoop.timeStepSize()
                         << " after check point unnecessarily!");

        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();
    }

    // check if setting the checkpoint at the current time shows error message
    {
        // redirect std::cerr
        std::stringstream cerrBuffer;
        std::streambuf* cerrOriginal = std::cerr.rdbuf(cerrBuffer.rdbuf());

        Dumux::CheckPointTimeLoop<double> timeLoop(tStart, dt, tEnd);
        timeLoop.setCheckPoint(tStart);

        // get result and reset buffer
        const auto result = cerrBuffer.str();
        std::cerr.rdbuf(cerrOriginal);
        std::cout << "Setting check point at the current time. Printed error message: " << result << std::endl;

        if (result.empty())
            DUNE_THROW(Dune::Exception, "Setting a checkpoint at the current time should print a warning to the error stream");
        if (Dune::FloatCmp::eq(timeLoop.timeStepSize(), 0.0, 1e-10))
            DUNE_THROW(Dune::Exception, "Time Loop reduced time step size to 0!");
    }

    // check if setting timestep to zero shows error message
    {
        // redirect std::cerr
        std::stringstream cerrBuffer;
        std::streambuf* cerrOriginal = std::cerr.rdbuf(cerrBuffer.rdbuf());

        Dumux::TimeLoop<double> timeLoop(tStart, dt, tEnd);
        timeLoop.setTimeStepSize(0.0);

        // get result and reset buffer
        const auto result = cerrBuffer.str();
        std::cerr.rdbuf(cerrOriginal);
        std::cout << "Setting zero time step. Printed error message: " << result << std::endl;

        if (result.empty())
            DUNE_THROW(Dune::Exception, "Setting a zero timeStepSize should print a warning to the error stream");
    }
}

void testWithDurations()
{
    using namespace std::chrono_literals;
    const auto _test = [&] (auto&& timeLoop) {
        if (!Dune::FloatCmp::eq(timeLoop.timeStepSize(), 0.001, 1e-10))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected initial time step size");
        if (!Dune::FloatCmp::eq(timeLoop.endTime(), 3600.0, 1e-10))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected end time");
        std::cout << "Time loop constructed as expected" << std::endl;

        timeLoop.setTimeStepSize(0.1ms);
        if (!Dune::FloatCmp::eq(timeLoop.timeStepSize(), 0.0001, 1e-10))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected time step size");
        std::cout << "Setting dt from a duration successful" << std::endl;

        timeLoop.reset(1ms, 0.5h, 2h);
        if (!Dune::FloatCmp::eq(timeLoop.timeStepSize(), 1800.0, 1e-10))
            DUNE_THROW(Dune::InvalidStateException, "Resetting time step size failed");
        if (!Dune::FloatCmp::eq(timeLoop.endTime(), 7200.0, 1e-10))
            DUNE_THROW(Dune::InvalidStateException, "Resetting end time failed");
        std::cout << "Resetting from durations successful" << std::endl;

        timeLoop.setTime(0.5h);
        timeLoop.setMaxTimeStepSize(1h);
        if (!Dune::FloatCmp::eq(timeLoop.time(), 1800.0, 1e-10))
            DUNE_THROW(Dune::InvalidStateException, "Setting time failed");
        if (!Dune::FloatCmp::eq(timeLoop.maxTimeStepSize(), 3600.0, 1e-10))
            DUNE_THROW(Dune::InvalidStateException, "Setting max time step size failed");
    };

    _test(Dumux::TimeLoop<double>{0s, 1ms, 1h});
    _test(Dumux::CheckPointTimeLoop<double>{0s, 1ms, 1h});

    // with CTAD
    _test(Dumux::TimeLoop{0s, 1ms, 1h});
    _test(Dumux::CheckPointTimeLoop{0s, 1ms, 1h});

    // try setting check points from durations
    Dumux::CheckPointTimeLoop timeLoop{0s, 1ms, 1h};
    timeLoop.setMaxTimeStepSize(1h);
    timeLoop.setPeriodicCheckPoint(0.5ms);
    if (!Dune::FloatCmp::eq(timeLoop.maxTimeStepSize(), 0.0005, 1e-10))
        DUNE_THROW(Dune::InvalidStateException, "Setting periodic checkpoint failed");
}

// After many time steps of constant size, the time misses the check points by accumulated round-off.
// Setting the time step size after every time step prevents the time loop from stretching
// a time step to reach a check point exactly.
void testCheckPointsAfterManyTimeSteps()
{
    std::cout << std::endl << "------- Test check points after many time steps ----------" << std::endl;
    const double tEnd = 1.0;
    for (const int numSteps : {4000, 10000})
    {
        for (const bool periodic : {false, true})
        {
            const double dt = tEnd/numSteps;
            Dumux::CheckPointTimeLoop<double> timeLoop(0.0, dt, tEnd, /*verbose=*/false);
            if (periodic)
                timeLoop.setPeriodicCheckPoint(0.1*tEnd);
            else
                timeLoop.setCheckPoint({0.03*tEnd, 0.1*tEnd, 0.3*tEnd, tEnd});

            int numCheckPoints = 0;
            timeLoop.start();
            while (!timeLoop.finished())
            {
                if (timeLoop.timeStepIndex() == numSteps)
                    DUNE_THROW(Dune::InvalidStateException, "Time loop not finished after " << numSteps
                                << " time steps of size " << dt << " (t = " << timeLoop.time() << ")");

                timeLoop.advanceTimeStep();
                timeLoop.setTimeStepSize(dt);
                if (timeLoop.isCheckPoint())
                    ++numCheckPoints;
            }

            if (timeLoop.timeStepIndex() != numSteps)
                DUNE_THROW(Dune::InvalidStateException, "Time loop finished after " << timeLoop.timeStepIndex()
                            << " instead of " << numSteps << " time steps");
            if (!timeLoop.isCheckPoint())
                DUNE_THROW(Dune::InvalidStateException, "Missed the check point at the end time");

            const int expectedNumCheckPoints = periodic ? 10 : 4;
            if (numCheckPoints != expectedNumCheckPoints)
                DUNE_THROW(Dune::InvalidStateException, "Reached " << numCheckPoints << " instead of "
                            << expectedNumCheckPoints << " check points");

            std::cout << "Reached all " << numCheckPoints << (periodic ? " periodic" : "")
                      << " check points in " << numSteps << " time steps" << std::endl;
        }
    }
}

int main(int argc, char* argv[])
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // unit time
    testTimeLoops(0.0, 1.0, 0.1);

    // microseconds
    testTimeLoops(0.0, 1.0e-6, 1e-7);

    // large time scales
    testTimeLoops(0.0, 1.0e12, 1e9, {1e11});

    // large time scales but small initial time step
    testTimeLoops(0.0, 1.0e12, 0.1, {1e11});

    testWithDurations();

    testCheckPointsAfterManyTimeSteps();

    return 0;
}
