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

#include <type_traits>
#include <dune/common/fvector.hh>

#include <dumux/common/typetraits/problem.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The base class for the sub-control entity-local evaluation of
 *        the terms of equations in the context of finite-volume schemes
 * \todo TODO: Should operators have a state? That is, be constructed and have non-static functions?
 */
template<class ElementVariables>
class FVOperators
{
    // The variables required for the evaluation of the equation
    using GridVariables = typename ElementVariables::GridVariables;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using ElementVolumeVariables = typename ElementVariables::ElementVolumeVariables;
    using ElementFluxVariablesCache = typename ElementVariables::ElementFluxVariablesCache;
    using Scalar = typename GridVariables::Scalar;

    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    static constexpr int numEq = PrimaryVariables::size();

    // The grid geometry on which the scheme operates
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    // user-input defined in the problem
    using Problem = std::decay_t<decltype(std::declval<GridVariables>().gridVolVars().problem())>;

public:
    //! export the type used to store scalar values for all equations
    //! TODO: How to obtain user-defined NumEqVector (ProblemTraits?)
    using NumEqVector = Dune::FieldVector<Scalar, numEq>;

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
     * \param scv The sub-control volume
     * \param volVars The primary & secondary variables evaluated for the scv
     * \note This must be overloaded by the implementation
     */
     static StorageTerm storage(const Problem& problem,
                                const SubControlVolume& scv,
                                const VolumeVariables& volVars)
    { DUNE_THROW(Dune::NotImplemented, "Storage operator not implemented!"); }

    /*!
     * \brief Compute the flux term of the equations for the given sub-control volume face
     * \param problem The problem to be solved (could store additionally required quantities)
     * \param element The grid element
     * \param fvGeometry The element-local view on the finite volume grid geometry
     * \param elemVolVars The element-local view on the grid volume variables
     * \param elemFluxVarsCache The element-local view on the grid flux variables cache
     * \param scvf The sub-control volume face for which the flux term is to be computed
     * \note This must be overloaded by the implementation
     */
    static FluxTerm flux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf)
    { DUNE_THROW(Dune::NotImplemented, "This model does not implement a flux method!"); }

    /*!
     * \brief Compute the source term of the equations for the given sub-control volume
     * \param problem The problem to be solved (could store additionally required quantities)
     * \param element The grid element
     * \param fvGeometry The element-local view on the finite volume grid geometry
     * \param elemVolVars The element-local view on the grid volume variables
     * \param scv The sub-control volume for which the source term is to be computed
     * \note This is a default implementation forwarding to interfaces in the problem
     */
     static SourceTerm source(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv)
    {
        SourceTerm source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        // multiply with scv volume
        source *= Extrusion::volume(scv)*elemVolVars[scv].extrusionFactor();

        return source;
    }
};

} // end namespace Dumux

#endif
