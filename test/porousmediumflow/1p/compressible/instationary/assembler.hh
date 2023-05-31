// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
};

} // end namespace Dumux::OnePCompressibleTest

#endif
