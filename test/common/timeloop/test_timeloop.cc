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

int main(int argc, char* argv[]) try
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

        timeLoop.setMaxTimeStepSize(0.033333333333333);
        if (std::abs(timeLoop.timeStepSize()-0.033333333333333) > 1e-15)
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

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
