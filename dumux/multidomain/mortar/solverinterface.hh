// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Interface for subdomain solvers in mortar-coupling models.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_SOLVER_INTERFACE_HH
#define DUMUX_MULTIDOMAIN_MORTAR_SOLVER_INTERFACE_HH

#include <cstddef>
#include <memory>

#include <dune/common/exceptions.hh>

#include "couplingmode.hh"
#include "trace.hh"

namespace Dumux::Mortar {

/*!
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
    using Trace = Mortar::Trace<GG, MortarGrid>;
    using MortarSolutionVector = MortarVector;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;

    virtual ~SubDomainSolver() = default;

    //! A solver for the subdomain with the given grid geometry
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

    /*!
     * \brief Set whether the subdomain is solved without the external data of its own
     *        problem, which turns its solution into the action of the linear operator on
     *        the mortar data.
     */
    virtual void setHomogeneous(bool homogeneous)
    {
        if (homogeneous)
            DUNE_THROW(Dune::NotImplemented, "This subdomain solver cannot drop the external data of its problem");
    }

    //! Set the mode in which mortar data enters this subdomain problem
    virtual void setCouplingMode(CouplingMode mode)
    {
        if (mode != CouplingMode::essential)
            DUNE_THROW(Dune::NotImplemented, "This subdomain solver only supports essential mortar coupling");
    }

    //! Return true if the entire boundary of this subdomain couples to mortar domains
    virtual bool isFloating() const
    { return false; }

    /*!
     * \brief The order of the trace data this solver accepts from the mortar domains: zero
     *        for one value per trace cell, one for one value per trace vertex.
     *
     * The model builds the projectors of this subdomain to produce data of this order.
     */
    virtual std::size_t traceDataOrder() const
    { return 0; }

    //! Return the underlying grid geometry
    const std::shared_ptr<const GridGeometry>& gridGeometry() const
    { return gridGeometry_; }

 private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
};

} // end namespace Dumux::Mortar

#endif
