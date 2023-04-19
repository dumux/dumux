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

#ifndef DUMUX_TEST_PYTHON_FVPROBLEM_BOUNDARY_LOOP_HH
#define DUMUX_TEST_PYTHON_FVPROBLEM_BOUNDARY_LOOP_HH

#include <string>
#include <memory>
#include <tuple>

#include <iostream>
#include <iomanip>
#include <dune/python/pybind11/iostream.h>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/python/pybind11/pybind11.h>

#include <dumux/common/boundarytypes.hh>
#include <dumux/python/common/boundarytypes.hh>

namespace Dumux::Python {

/////////////////////////////////////////////////////////
/// Some debugging/testing  stuff
///////////////////////////////////////////////////////////
template<class Problem>
void printProblemTest(const Problem& problem)
{
    std::size_t numNeumann = 0;
    std::size_t numDirichlet = 0;
    double totalSource = 0;
    const auto& gg = problem.gridGeometry();
    auto fvGeometry = localView(gg);
    for (const auto& element : elements(gg.gridView()))
    {
        fvGeometry.bindElement(element);

        for (const auto& scv : scvs(fvGeometry))
        {
            const auto boundaryTypes = problem.boundaryTypes(element, scv);

            if (boundaryTypes.hasDirichlet())
                numDirichlet += 1;
            else if (boundaryTypes.hasNeumann())
                numNeumann += 1;

            totalSource += problem.sourceAtPos(scv.center())[0]*scv.volume();
        }
    }

    std::cout << "[cpp] Found " << numNeumann << " Neumann faces and " << numDirichlet << " Dirichlet faces" << std::endl;
    std::cout << "[cpp] Total source " << totalSource << std::fixed << std::setprecision(3) << " kg/s" << std::endl;
}

template<class P>
class PrintProblemTest
{
public:
    using Problem = P;

    PrintProblemTest(std::shared_ptr<const Problem> problem)
    : problem_(problem)
    {}

    void print()
    {
        printProblemTest(*problem_);
    }

private:
    std::shared_ptr<const Problem> problem_;
};

template<class T, class... options>
void registerPrintProblemTest(pybind11::handle scope, pybind11::class_<T, options...> cls)
{
    using Problem = typename T::Problem;
    cls.def(pybind11::init([](std::shared_ptr<const Problem> problem){
        return std::make_unique<T>(problem);
    }));
    cls.def("print", [](T& self){
        pybind11::scoped_ostream_redirect stream(std::cout, pybind11::module::import("sys").attr("stdout"));
        self.print();
    });
}

} // end namespace Dumux::Python

#endif
