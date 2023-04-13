// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The model traits used in the tracer facet coupling test.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_TRACER_MODELTRAITS_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_TRACER_MODELTRAITS_HH

#include <dumux/porousmediumflow/tracer/model.hh>

namespace Dumux {

//! Custom model traits disabling diffusion
template<int nComp, bool useMol, bool enableCompDisp>
struct TracerTestModelTraits : public TracerModelTraits<nComp, useMol, enableCompDisp>
{
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableCompositionalDispersion() { return enableCompDisp; }
};

} // end namespace Dumux

#endif
