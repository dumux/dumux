# TimeLoop

In case of stationary PDEs, a time loop is not necessary. But in the case of non-stationary PDEs, an object of the class `TimeLoop` is instanciated which stores all parameters connected to the time management of the PDE such as time step width, current simulation time, total simulation time and so on. The base class of the time loop is implemented in [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/common/timeloop.hh).

### Key functionalities

* time():
    - Return the current simulation time, which is the time of the old/current time step before the current integration step.
* timeStepSize():
    - Return the suggested time step size.
* maxTimeStepSize():
    - Get the maximum possible time step size.
* advanceTimeStep():
    - Advance to the next time step.
* setTimeStepSize(dt):
    - Set the current time step size to a given value `dt`.
* finished():
    - Return true if the simulation is finished.
* reset(...):
    - Reset the time loop to provided starting time, end time, time step size etc.
* endTime():
    - Return the number of (simulated) seconds which the simulation runs.
* wallClockTime():
    - Return the current wall clock time (cpu time) spent on the simluation so far.
* timeStepIndex():
    - Return number of time steps which have been executed since the beginning of the simulation.
* reportTimeStep():
    - State info on cpu time.
* finalize():
    - Print final status and stops tracking the time.
* verbose():
    - Return if the time loop has verbose output or not.
* setVerbose():
    - Set the verbosity of the time loop.
