//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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
        std::cout << "Setting check point at the current time printed '" << result << "' to std::cerr" << std::endl;

        if (result.empty())
            DUNE_THROW(Dune::Exception, "Setting a checkpoint at the current time should print a warning to std::cerr");
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
        std::cout << "Setting zero time step printed '" << result << "' to std::cerr" << std::endl;

        if (result.empty())
            DUNE_THROW(Dune::Exception, "Setting a zero timeStepSize should print a warning to std::cerr");
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

    return 0;
}
