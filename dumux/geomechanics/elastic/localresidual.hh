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
 * \ingroup Geomechanics
 * \brief Element-wise calculation of the local residual for problems
 *        using the elastic model considering linear elasticity.
 */
#ifndef DUMUX_GEOMECHANICS_ELASTIC_LOCAL_RESIDUAL_HH
#define DUMUX_GEOMECHANICS_ELASTIC_LOCAL_RESIDUAL_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/method.hh>
#include <dumux/assembly/felocalresidual.hh>

namespace Dumux {

//! Forward declaration of the implementation
template<class TypeTag, DiscretizationMethod dm>
class ElasticLocalResidualImpl;

/*!
 * \ingroup Geomechanics
 * \brief Element-wise calculation of the local residual for problems
 *        using the elastic model considering linear elasticity. This
 *        is a convenience alias the selects the right implementation
 *        depending on the chosen discretization scheme.
 * \note Currently the elastic model works for box discretization
 *       and finite element schemes
 */
template<class TypeTag>
using ElasticLocalResidual
= ElasticLocalResidualImpl<TypeTag, GetPropType<TypeTag, Properties::GridGeometry>::discMethod>;

/*!
 * \ingroup Geomechanics
 * \brief Implementation of the local residual for finite element schemes.
 */
template<class TypeTag>
class ElasticLocalResidualImpl<TypeTag, DiscretizationMethod::fem>
: public FELocalResidual<TypeTag>
{
    using ParentType = FELocalResidual<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SecondaryVariables = typename GridVariables::SecondaryVariables;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    // class assembling the stress tensor
    using StressType = GetPropType<TypeTag, Properties::StressType>;

public:
    using ParentType::ParentType;

    /*!
     * \brief Calculate the flux term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param feGeometry The finite-element geometry
     * \param elemSol The element solution vector
     * \param ipData The trial and ansatz space shape function values/gradients
     *               evaluated at the integration point
     * \param secVars The secondary variables evaluated at the integration point
     */
    template<class IpData, class ElementSolution>
    typename StressType::StressTensor computeFlux(const Problem& problem,
                                                  const Element& element,
                                                  const ElementGeometry& feGeometry,
                                                  const ElementSolution& elemSol,
                                                  const IpData& ipData,
                                                  const SecondaryVariables& secVars) const
    {
        // obtain force on the face from stress type
        return StressType::stressTensor(problem, element, feGeometry, elemSol, ipData, secVars);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param feGeometry The finite-element geometry
     * \param elemSol The element solution vector
     * \param ipData The trial and ansatz space shape function values/gradients
     *               evaluated at the integration point
     * \param secVars The secondary variables evaluated at the integration point
     */
    template<class IpData, class ElementSolution>
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const ElementGeometry& feGeometry,
                              const ElementSolution& elemSol,
                              const IpData& ipData,
                              const SecondaryVariables& secVars) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, feGeometry, elemSol, ipData, secVars);

        // TODO: add contribution from possible point sources

        // maybe add gravitational acceleration
        static const bool gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (gravity)
        {
            const auto& g = problem.spatialParams().gravity(ipData.ipGlobal());
            for (int dir = 0; dir < GridView::dimensionworld; ++dir)
                source[Indices::momentum(dir)] += secVars.solidDensity()*g[dir];
        }

        return source;
    }
};

/*!
 * \ingroup Geomechanics
 * \brief Implementation of the local residual for the box scheme.

 */
template<class TypeTag>
class ElasticLocalResidualImpl<TypeTag, DiscretizationMethod::box>
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;

    // class assembling the stress tensor
    using StressType = GetPropType<TypeTag, Properties::StressType>;

public:
    using ParentType::ParentType;

    /*!
     * \brief For the elastic model the storage term is zero since
     *        we neglect inertia forces.
     *
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluates the force in all grid directions acting
     *        on a face of a sub-control volume.
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux compuation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // obtain force on the face from stress type
        return StressType::force(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVolVars The volume variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     * \note This is the default implementation for geomechanical models adding to
     *       the user defined sources the source stemming from the gravitational acceleration.
     *
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        // maybe add gravitational acceleration
        static const bool gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (gravity)
        {
            const auto& g = problem.spatialParams().gravity(scv.center());
            for (int dir = 0; dir < GridView::dimensionworld; ++dir)
                source[Indices::momentum(dir)] += elemVolVars[scv].solidDensity()*g[dir];
        }

        return source;
    }
};

} // end namespace Dumux

#endif
