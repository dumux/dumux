// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
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
 * \brief Solves a poisson problem discretized using TPFA.
 *        This exposes the minimal requirements for adding a new model.
 */
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/matrixpattern.hh>


int main (int argc, char *argv[])
{
    Dumux::initialize(argc, argv);

    Dumux::MatrixPattern pattern;
    pattern = Dumux::MatrixPattern{10, 10};
    pattern.resize(2, 2);

    pattern.add(0, 0);
    pattern.add(0, 1);
    pattern.add(1, 1);

    auto indexRange = indices(pattern);
    auto indexIterator = indexRange.begin();
    if (indexIterator->rowIndex != 0) DUNE_THROW(Dune::InvalidStateException, "Unexpected Index");
    if (indexIterator->colIndex != 0) DUNE_THROW(Dune::InvalidStateException, "Unexpected Index");

    ++indexIterator;
    if (indexIterator->rowIndex != 0) DUNE_THROW(Dune::InvalidStateException, "Unexpected Index");
    if (indexIterator->colIndex != 1) DUNE_THROW(Dune::InvalidStateException, "Unexpected Index");

    ++indexIterator;
    if (indexIterator->rowIndex != 1) DUNE_THROW(Dune::InvalidStateException, "Unexpected Index");
    if (indexIterator->colIndex != 1) DUNE_THROW(Dune::InvalidStateException, "Unexpected Index");

    if (++indexIterator != indexRange.end())
        DUNE_THROW(Dune::InvalidStateException, "Expected end of index range");

    if (pattern.nnz() != 3)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of nonzero entries");

    // ensure that adding the same index twice does not affect the nonzeroes
    pattern.add(0, 1);
    if (pattern.nnz() != 3)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of nonzero entries");

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
