// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief TODO: docme!
 */

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

    cls.def_property_readonly("time", &TimeLoop::time);
    cls.def_property_readonly("timeStepSize", &TimeLoop::timeStepSize);
    cls.def_property_readonly("finished", &TimeLoop::finished);
    cls.def_property_readonly("isCheckPoint", &TimeLoop::isCheckPoint);
    cls.def("start", &TimeLoop::start);
    cls.def(
        "reset",
        static_cast<void(TimeLoop::*)(Scalar, Scalar, Scalar, bool)>(&TimeLoop::reset),
        "startTime"_a, "dt"_a, "endTime"_a, "verbose"_a=true
    );
    cls.def("advanceTimeStep", &TimeLoop::advanceTimeStep);
    cls.def("setTimeStepSize", static_cast<void(TimeLoop::*)(Scalar)>(&TimeLoop::setTimeStepSize), "dt"_a);
    cls.def("setMaxTimeStepSize", &TimeLoop::template setMaxTimeStepSize<Scalar>, "dt"_a);
    cls.def("reportTimeStep", &TimeLoop::reportTimeStep);
    cls.def("finalize", [](TimeLoop& self){ self.finalize(); });
    cls.def(
        "setPeriodicCheckPoint",
        &TimeLoop::template setPeriodicCheckPoint<Scalar, Scalar>,
        "interval"_a, "offset"_a=0.0
    );
    cls.def("setCheckPoints", [](TimeLoop& self, const std::vector<Scalar>& checkPoints) {
        self.setCheckPoint(checkPoints.begin(), checkPoints.end());
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
