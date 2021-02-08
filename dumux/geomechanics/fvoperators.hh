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
 * \ingroup PorousmediumflowModels
 * \brief Sub-control entity-local evaluation of the operators
 *        within in the systems of equations of n-phase immiscible models.
 */
#ifndef DUMUX_FV_GEOMECHANICS_OPERATORS_HH
#define DUMUX_FV_GEOMECHANICS_OPERATORS_HH

#include <type_traits>

#include <dumux/assembly/fv/operators.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Sub-control entity-local evaluation of the operators
 *        within in the systems of equations of n-phase immiscible models.
 * \tparam ModelTraits defines model-related types and variables (e.g. number of phases)
 * \tparam ElementVariables the type of element-local view on the grid variables
 * \tparam StressType the constitutive model for stress computation at scv faces.
 */
template<class ModelTraits, class ElementVariables, class StressType>
class FVGeomechanicsOperators
: public FVOperators<ElementVariables>
{
    using ParentType = FVOperators<ElementVariables>;

    // The variables required for the evaluation of the equation
    using GridVariables = typename ElementVariables::GridVariables;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using ElementVolumeVariables = typename ElementVariables::ElementVolumeVariables;
    using ElementFluxVariablesCache = typename ElementVariables::ElementFluxVariablesCache;

    // The grid geometry on which the scheme operates
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using Problem = std::decay_t<decltype(std::declval<GridVariables>().gridVolVars().problem())>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int conti0EqIdx = ModelTraits::Indices::conti0EqIdx;

public:
    //! export the type used to store scalar values for all equations
    using typename ParentType::NumEqVector;

    // export the types of the flux/source/storage terms
    using typename ParentType::FluxTerm;
    using typename ParentType::SourceTerm;
    using typename ParentType::StorageTerm;

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
    { return StorageTerm(0.0); }

    // TODO: overload source here to add bulk gravity term!

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
    { return StressType::force(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache); }
};

} // end namespace Dumux

#endif
