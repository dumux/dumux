#include <iostream>
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>
#include <dune/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/timeloop.hh>

void testTimeLoopWithSmallCheckpoints() {
    std::cout << "\n------- Test Time Loop with Checkpoints ----------" << std::endl;

    double tStart = 1.0;
    double tEnd = 5.0;
    double dt = 1.0;
    double epsilon1 = 1e-15;
    double epsilon2 = 1e-16;
    double epsilon3 = 1e-17;

    Dumux::CheckPointTimeLoop<double> timeLoop(tStart, dt, tEnd);
    timeLoop.setCheckPoint(2.0 + epsilon1);
    timeLoop.setCheckPoint(3.0 + epsilon2);
    timeLoop.setCheckPoint(4.0 + epsilon3);

    timeLoop.start();

    while (!timeLoop.finished()) {
        timeLoop.advanceTimeStep();
        timeLoop.reportTimeStep();

        if (std::abs(timeLoop.timeStepSize()) < 1e-14) {
            DUNE_THROW(Dune::InvalidStateException, "Time step size is too small: " << timeLoop.timeStepSize());
        }
    }

    timeLoop.finalize();
}

int main(int argc, char* argv[]) {
    Dumux::initialize(argc, argv);
    testTimeLoopWithSmallCheckpoints();
    return 0;
}
