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
 * \brief Interface for subdomain solvers in mortar-coupling models.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_SOLVER_INTERFACE_HH
#define DUMUX_MULTIDOMAIN_MORTAR_SOLVER_INTERFACE_HH

#include <memory>

#include <dumux/io/grid/facetgridmanager.hh>

namespace Dumux::Mortar {

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Abstract base class for subdomain solvers.
 * \tparam MortarVector The type used to represent the solution in the mortar domain.
 * \tparam MortarGrid The grid type used to represent the mortar domain.
 * \tparam GG The subdomain grid geometry
 */
template<typename MortarVector,
         typename MortarGrid,
         typename GG>
class SubDomainSolver
{
 public:
    using GridGeometry = GG;
    using Trace = FacetGridManager<typename GG::GridView::Grid, MortarGrid>;
    using MortarSolutionVector = MortarVector;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;

    virtual ~SubDomainSolver() = default;
    explicit SubDomainSolver(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_{std::move(gridGeometry)}
    {}

    //! Solve the subdomain problem
    virtual void solve() = 0;

    //! Set the mortar boundary condition for the mortar with the given id
    virtual void setTraceVariables(std::size_t, MortarSolutionVector) = 0;

    //! Register a trace coupling to the mortar with the given id
    virtual void registerMortarTrace(std::shared_ptr<const Trace>, std::size_t) = 0;

    //! Assemble the variables on the trace overlapping with the given mortar domain
    virtual MortarSolutionVector assembleTraceVariables(std::size_t) const = 0;

    //! Return the underlying grid geometry
    const std::shared_ptr<const GridGeometry>& gridGeometry() const
    { return gridGeometry_; }

 private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
};

} // end namespace Dumux::Mortar

#endif
