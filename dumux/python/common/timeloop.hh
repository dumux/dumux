#ifndef DUMUX_PYTHON_COMMON_TIMELOOP_HH
#define DUMUX_PYTHON_COMMON_TIMELOOP_HH

#include <dumux/common/timeloop.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dumux::Python {

template <class Scalar, class... options>
void registerTimeLoop(pybind11::handle scope,
                      pybind11::class_<CheckPointTimeLoop<Scalar>, options...> cls)
{
    using pybind11::operator""_a;

    using TimeLoop = CheckPointTimeLoop<Scalar>;
    cls.def(pybind11::init([](Scalar startTime, Scalar dt, Scalar endTime, bool verbose){
        return new TimeLoop(startTime, dt, endTime, verbose);
    }), "startTime"_a, "dt"_a, "endTime"_a, "verbose"_a=true);

    cls.def("time", &TimeLoop::time);
    cls.def("timeStepSize", &TimeLoop::timeStepSize);
    cls.def("start", &TimeLoop::start);
    cls.def("reset", &TimeLoop::reset, "startTime"_a, "dt"_a, "endTime"_a, "verbose"_a=true);
    cls.def("advanceTimeStep", &TimeLoop::advanceTimeStep);
    cls.def("setTimeStepSize", &TimeLoop::setTimeStepSize, "dt"_a);
    cls.def("setMaxTimeStepSize", &TimeLoop::setMaxTimeStepSize, "dt"_a);
    cls.def("finished", &TimeLoop::finished);
    cls.def("reportTimeStep", &TimeLoop::reportTimeStep);
    cls.def("finalize", [](TimeLoop& self){ self.finalize(); });
    cls.def("setPeriodicCheckPoint", &TimeLoop::setPeriodicCheckPoint, "interval"_a, "offset"_a=0.0);
    cls.def("isCheckPoint", &TimeLoop::isCheckPoint);
    cls.def("setCheckPoints", [](TimeLoop& self, const std::vector<double>& checkPoints) {
        self.setCheckPoint(checkPoints);
    });
}

template<class Scalar>
void registerTimeLoop(pybind11::handle scope, const char *clsName = "TimeLoop")
{
    pybind11::class_<CheckPointTimeLoop<Scalar>> cls(scope, clsName);
    registerTimeLoop(scope, cls);
}

} // namespace Dumux::Python

#endif
