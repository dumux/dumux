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
 * \ingroup Assembly
 * \brief The base class for the sub-control entity-local evaluation of
 *        the terms of equations in the context of finite-volume schemes
 */
#ifndef DUMUX_FV_OPERATORS_HH
#define DUMUX_FV_OPERATORS_HH

#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/localcontext.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Assembly
 * \brief Convenience alias to define the context finite-volume operators work on.
 * \tparam EV The element-(stencil)-local variables
 */
template<class EV>
using FVOperatorsContext = DefaultLocalContext<EV>;

/*!
 * \ingroup Assembly
 * \brief The base class for the sub-control entity-local evaluation of
 *        the terms of equations in the context of finite-volume schemes
 * \tparam EV The element-(stencil)-local variables
 */
template<class EV>
class FVOperators
{
    // context type on which to operate
    using LC = FVOperatorsContext<EV>;

    // The grid geometry on which the scheme operates
    using FVElementGeometry = typename LC::ElementGridGeometry;
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    // The variables (primary/secondary) on the grid
    using GridVariables = typename LC::ElementVariables::GridVariables;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;

public:
    //! export the local context on which this operates
    using LocalContext = LC;

    //! export the type used to store scalar values for all equations
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    // export the types of the flux/source/storage terms
    using FluxTerm = NumEqVector;
    using SourceTerm = NumEqVector;
    using StorageTerm = NumEqVector;

    /*!
     * \name Model specific interfaces
     * \note The following method are the model specific implementations of the
     *       operators appearing in the equation and should be overloaded by the
     *       implementations.
     */
    // \{

    /*!
     * \brief Compute the storage term of the equations for the given sub-control volume
     * \param problem The problem to be solved (could store additionally required quantities)
     * \param context The element-stencil-local required data
     * \param scv The sub-control volume
     * \note This must be overloaded by the implementation
     * \note In cell-centered schemes, depending on the assembler or the policy used
     *       therein, the context may contain variables (e.g. from neighboring cells
     *       within the stencil) with respect to which the storage term is not
     *       differentiated. This means that if a Newton solver is used to solve
     *       the system of equations, one ends up with a quasi-Newton scheme. If
     *       this is undesireable, make sure to set up an assembler that takes into
     *       account the dependencies of the storage term w.r.t variables in the context.
     */
     template<class Problem>
     static StorageTerm storage(const Problem& problem,
                                const LocalContext& context,
                                const SubControlVolume& scv)
    { DUNE_THROW(Dune::NotImplemented, "Storage operator not implemented!"); }

    /*!
     * \brief Compute the flux term of the equations for the given sub-control volume face
     * \param problem The problem to be solved (could store additionally required quantities)
     * \param context The element-stencil-local required data
     * \param scvf The sub-control volume face for which the flux term is to be computed
     * \note This must be overloaded by the implementation
     */
    template<class Problem>
    static FluxTerm flux(const Problem& problem,
                         const LocalContext& context,
                         const SubControlVolumeFace& scvf)
    { DUNE_THROW(Dune::NotImplemented, "This model does not implement a flux method!"); }

    /*!
     * \brief Compute the source term of the equations for the given sub-control volume
     * \param problem The problem to be solved (could store additionally required quantities)
     * \param context The element-stencil-local required data
     * \param scv The sub-control volume for which the source term is to be computed
     * \note This is a default implementation forwarding to interfaces in the problem
     * \note In cell-centered schemes, depending on the assembler or the policy used
     *       therein, the context may contain variables (e.g. from neighboring cells
     *       within the stencil) with respect to which the source term is not
     *       differentiated. This means that if a Newton solver is used to solve
     *       the system of equations, one ends up with a quasi-Newton scheme. If
     *       this is undesireable, make sure to set up an assembler that takes into
     *       account the dependencies of the storage term w.r.t variables in the context.
     */
     template<class Problem>
     static SourceTerm source(const Problem& problem,
                              const LocalContext& context,
                              const SubControlVolume& scv)
    {
        SourceTerm source(0.0);

        // add contributions from volume flux sources
        source += problem.source(context, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(context, scv);

        // multiply with scv volume
        const auto& elemVolVars = context.elementVariables().elemVolVars();
        source *= Extrusion::volume(scv)*elemVolVars[scv].extrusionFactor();

        return source;
    }
};

} // end namespace Dumux::Experimental

#endif
