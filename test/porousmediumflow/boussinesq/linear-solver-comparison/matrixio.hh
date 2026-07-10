// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Filenames and MatrixMarket export helpers shared between the
 *        matrix-generation and matrix-solving stages of the Boussinesq
 *        linear-solver benchmark.
 *
 * The generation stage (main_generate_{pressure,vorticity}.cc) writes one
 * matrix/rhs pair per checkpoint time; the solving stage
 * (main_solve_{pressure,vorticity}.cc) reads them back in. Keeping the naming
 * convention in one place ensures both stages agree on it.
 */
#ifndef DUMUX_BOUSSINESQ_MATRIXIO_HH
#define DUMUX_BOUSSINESQ_MATRIXIO_HH

#include <cmath>
#include <string>

#include <dune/istl/matrixmarket.hh>

namespace Dumux {

// Checkpoint times are expected to be (near-)integers (e.g. 1, 2, 3, 4);
// they are rounded when turned into a filename.
inline std::string matrixFileName(const std::string& problemName, double time)
{ return problemName + "_t" + std::to_string(static_cast<long long>(std::llround(time))) + "_matrix.mm"; }

inline std::string rhsFileName(const std::string& problemName, double time)
{ return problemName + "_t" + std::to_string(static_cast<long long>(std::llround(time))) + "_rhs.mm"; }

template<class Matrix, class Vector>
void exportSystem(const std::string& problemName, double time, const Matrix& A, const Vector& b)
{
    Dune::storeMatrixMarket(A, matrixFileName(problemName, time));
    Dune::storeMatrixMarket(b, rhsFileName(problemName, time));
}

} // end namespace Dumux

#endif