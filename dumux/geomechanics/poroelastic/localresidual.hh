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
 * \brief Element-wise calculation of the local residual
 *        for problems using the poroelastic model.
 */
#ifndef DUMUX_PRORELASTIC_LOCAL_RESIDUAL_HH
#define DUMUX_PRORELASTIC_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/geomechanics/elastic/localresidual.hh>

namespace Dumux {

//! Forward declaration of the implementation
template<class TypeTag, DiscretizationMethod dm>
class PoroElasticLocalResidualImpl;

/*!
 * \ingroup Geomechanics
 * \brief Element-wise calculation of the local residual for problems
 *        using the poro-elastic model. This is a convenience alias
 *        that selects the right implementation depending on the
 *        chosen discretization scheme.
 * \note Currently the poro-elastic model works for box
 *       discretization and finite element schemes
 */
template<class TypeTag>
using PoroElasticLocalResidual
= PoroElasticLocalResidualImpl<TypeTag, GetPropType<TypeTag, Properties::GridGeometry>::discMethod>;

/*!
 * \ingroup Geomechanics
 * \brief Implementation of the local residual for the box scheme.
 */
template<class TypeTag>
class PoroElasticLocalResidualImpl<TypeTag, DiscretizationMethod::box>
: public ElasticLocalResidual<TypeTag>
{
    using ParentType = ElasticLocalResidual<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

public:
    using ParentType::ParentType;

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
                              const SubControlVolume& scv) const
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
            // compute average density
            const auto& vv = elemVolVars[scv];
            const auto phi = vv.porosity();
            const auto rhoFluid = problem.effectiveFluidDensity(element, scv);
            const auto rhoAverage = phi*rhoFluid + (1.0 - phi)*vv.solidDensity();

            // add body force
            const auto& g = problem.spatialParams().gravity(scv.center());
            for (int dir = 0; dir < GridView::dimensionworld; ++dir)
                source[ Indices::momentum(dir) ] += rhoAverage*g[dir];
        }

        return source;
    }
};

/*!
 * \ingroup Geomechanics
 * \brief Implementation of the local residual for finite element schemes.
 */
template<class TypeTag>
class PoroElasticLocalResidualImpl<TypeTag, DiscretizationMethod::fem>
: public ElasticLocalResidual<TypeTag>
{
    using ParentType = ElasticLocalResidual<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FEElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SecondaryVariables = typename GetPropType<TypeTag, Properties::GridVariables>::SecondaryVariables;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
public:
    using ParentType::ParentType;

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
                              const FEElementGeometry& feGeometry,
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
            // compute average density
            const auto phi = secVars.porosity();
            const auto rhoFluid = problem.effectiveFluidDensity(element, ipData, elemSol);
            const auto rhoAverage = phi*rhoFluid + (1.0 - phi)*secVars.solidDensity();

            // add body force
            const auto& g = problem.spatialParams().gravity(ipData.ipGlobal());
            for (int dir = 0; dir < GridView::dimensionworld; ++dir)
                source[ Indices::momentum(dir) ] += rhoAverage*g[dir];
        }

        return source;
    }
};

} // end namespace Dumux

#endif
