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

#ifndef DUMUX_PYTHON_IO_VTK_OUTPUTMODULE_HH
#define DUMUX_PYTHON_IO_VTK_OUTPUTMODULE_HH

#include <dumux/io/vtkoutputmodule.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/python/common/volumevariables.hh>
#include <dumux/io/velocityoutput.hh>

namespace Dumux::Python {

template <class GridVariables, class SolutionVector, class... options>
void registerVtkOutputModule(pybind11::handle scope,
                             pybind11::class_<VtkOutputModule<GridVariables, SolutionVector>, options...> cls)
{
    using pybind11::operator""_a;

    using VtkOutputModule = Dumux::VtkOutputModule<GridVariables, SolutionVector>;
    using VolumeVariables = typename VtkOutputModule::VolumeVariables;
    Dumux::Python::Impl::registerVolumeVariables<VolumeVariables>(scope);


    cls.def(pybind11::init([](const GridVariables& gridVariables,
                              const SolutionVector& sol,
                              const std::string& name){
        return new VtkOutputModule(gridVariables, sol, name);
    }));

    using Scalar = double;

    cls.def("addField", [](VtkOutputModule& self, const SolutionVector& sol, const std::string& name){
        self.addField(sol, name);
    });

    cls.def("write", [](VtkOutputModule& self, Scalar time){
        self.write(time);
    });

    cls.def("addVolumeVariable", [](VtkOutputModule& self,
                                    std::function<Scalar(const VolumeVariables&)>&& f,
                                    const std::string& name){
        self.addVolumeVariable(std::move(f), name);
    });

    using VelocityOutputType = Dumux::VelocityOutput<GridVariables>;
    cls.def("addVelocityOutput", [](VtkOutputModule& self, std::shared_ptr<VelocityOutputType> velocityOutput){
        self.addVelocityOutput(velocityOutput);
    });
};


template<class GridVariables, class SolutionVector>
void registerVtkOutputModule(pybind11::handle scope, const char *clsName = "VtkOutputModule")
{
    using VtkOutputModule = Dumux::VtkOutputModule<GridVariables, SolutionVector>;
    pybind11::class_<VtkOutputModule> cls(scope, clsName);
    registerVtkOutputModule(scope, cls);
}

} // namespace Dumux::Python

#endif
