// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Preconditioners for mortar-coupling models.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_PRECONDITIONERS_HH
#define DUMUX_MULTIDOMAIN_MORTAR_PRECONDITIONERS_HH


#include <dune/common/exceptions.hh>
#include <dune/istl/preconditioner.hh>

namespace Dumux::Mortar {

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Preconditioner implementation that does nothing.
 */
template<typename X>
struct NoPreconditioner : public Dune::Preconditioner<X, X>
{
    explicit NoPreconditioner() {}

    void pre (X&, X&) override {}
    void apply (X& r, const X& x) override { r = x; }
    void post (X&) override {}

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }
};


} // end namespace Dumux::Mortar

#endif
