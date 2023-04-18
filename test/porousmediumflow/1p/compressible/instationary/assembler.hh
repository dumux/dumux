// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * \ingroup OnePTests
 * \brief Dummy assembler that fulfills the new experimental assembly interface.
 */
#ifndef DUMUX_COMPRESSIBLE_ONEP_TEST_ASSEMBLER_HH
#define DUMUX_COMPRESSIBLE_ONEP_TEST_ASSEMBLER_HH

namespace Dumux::OnePCompressibleTest {

// Custom assembler to test assembly with grid variables,
// fulfilling the foreseen required interface
template<class Assembler>
class GridVarsAssembler : public Assembler
{
public:
    using Assembler::Assembler;
    using typename Assembler::GridVariables;
    using typename Assembler::ResidualType;

    using Variables = GridVariables;

    void assembleJacobianAndResidual(const GridVariables& gridVars)
    { Assembler::assembleJacobianAndResidual(gridVars.dofs()); }

    void assembleResidual(const GridVariables& gridVars)
    { Assembler::assembleResidual(gridVars.dofs()); }

    //! Remove residualNorm once this is fixed in the solvers !2113
    auto residualNorm(const GridVariables& gridVars)
    { return Assembler::residualNorm(gridVars.dofs()); }
};

} // end namespace Dumux::OnePCompressibleTest

#endif
