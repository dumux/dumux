// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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
