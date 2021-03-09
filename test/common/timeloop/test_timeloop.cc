// time loop test
#include <config.h>

#include <iostream>
#include <cmath>
#include <chrono>
#include <thread>
#include <limits>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/timeloop.hh>

int main(int argc, char* argv[])
{
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    //! Standard time loop
    double tStart = 0; double tEnd = 1; double dt = 0.1;

    //! Standard time loop
    {
        if (mpiHelper.rank() == 0) std::cout << "------- Test time loop ----------" << std::endl;
        Dumux::TimeLoop<double> timeLoop(tStart, dt, tEnd);

        if (std::abs(timeLoop.time()-tStart) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Wrong time at start!");
        if (std::abs(timeLoop.endTime()-tEnd) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Wrong end time: " << timeLoop.endTime() << " should be " << tEnd);
        if (std::abs(timeLoop.timeStepSize()-dt) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Wrong time step size!");

        timeLoop.setMaxTimeStepSize(0.03333333333333333);
        if (std::abs(timeLoop.timeStepSize()-0.03333333333333333) > 1e-15)
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
            timeLoop.setTimeStepSize(timeLoop.timeStepSize());
        }

        timeLoop.finalize();

        if (std::abs(timeLoop.time()-tEnd) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Ended with wrong end time!");
    }

    // check point timeLoop
    {
        if (mpiHelper.rank() == 0) std::cout << std::endl << "------- Test check point time loop ----------" << std::endl;
        Dumux::CheckPointTimeLoop<double> timeLoop(tStart, dt, tEnd);
        timeLoop.setPeriodicCheckPoint(0.03, -0.03);
        if (!timeLoop.isCheckPoint())
            DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

        timeLoop.start();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();
        timeLoop.setPeriodicCheckPoint(0.03);

        if (std::signbit(timeLoop.timeStepSize()))
            DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());

        if (!timeLoop.isCheckPoint())
            DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

        timeLoop.setCheckPoint(0.02); // can't be set because it's in the past
        timeLoop.setCheckPoint(0.08);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();

        if (std::signbit(timeLoop.timeStepSize()))
            DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());

        if (!timeLoop.isCheckPoint())
            DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

        timeLoop.setEndTime(1e9);
        timeLoop.setTimeStepSize(1e7);

        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();

        if (std::signbit(timeLoop.timeStepSize()))
            DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());

        if (!timeLoop.isCheckPoint())
            DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

        timeLoop.removeAllCheckPoints();
        timeLoop.setTimeStepSize(0.5e-6);
        timeLoop.setPeriodicCheckPoint(1.0e6, 1.0e6);

        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();

        if (std::signbit(timeLoop.timeStepSize()))
            DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());

        while (!timeLoop.finished())
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            timeLoop.advanceTimeStep();
            timeLoop.reportTimeStep();
            timeLoop.setTimeStepSize(std::numeric_limits<double>::max());

            if (std::abs(timeLoop.time() - 2.0e6) < 1e-15 && !timeLoop.isCheckPoint())
                DUNE_THROW(Dune::InvalidStateException, "This is supposed to be a check point! t = " << timeLoop.time());

            if (timeLoop.timeStepIndex() == 10)
                timeLoop.setPeriodicCheckPoint(1.0e8);

            if (std::signbit(timeLoop.timeStepSize()))
                DUNE_THROW(Dune::InvalidStateException, "Time step size is negative! dt = " << timeLoop.timeStepSize());
        }

        timeLoop.finalize();

        if (std::abs(timeLoop.time()-1e9) > 1e-15)
            DUNE_THROW(Dune::InvalidStateException, "Ended with wrong end time!");

    }

    // check if time loops remembers time step before check point
    {
        if (mpiHelper.rank() == 0) std::cout << std::endl << "------- Test check point time loop ----------" << std::endl;
        Dumux::CheckPointTimeLoop<double> timeLoop(tStart, dt, tEnd);
        timeLoop.setCheckPoint(0.101);
        timeLoop.start();
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();
        if (!(std::abs(timeLoop.timeStepSize()-0.1) < 1e-14))
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
        if (Dune::FloatCmp::eq(timeLoop.timeStepSize(), 0.0, 1e-10*tEnd))
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

    return 0;
}
