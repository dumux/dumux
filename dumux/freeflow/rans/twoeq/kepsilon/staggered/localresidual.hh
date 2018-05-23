// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
  * \ingroup KepsilonModel
  * \copydoc Dumux::KepsilonResidualImpl
  */
#ifndef DUMUX_STAGGERED_KEPSILON_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_KEPSILON_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseLocalResidual, DiscretizationMethod discMethod>
class KEpsilonResidualImpl;

/*!
 * \ingroup KEpsilonModel
 * \brief Element-wise calculation of the residual for k-epsilon models using the staggered discretization
 */
template<class TypeTag, class BaseLocalResidual>
class KEpsilonResidualImpl<TypeTag, BaseLocalResidual, DiscretizationMethod::staggered>
: public BaseLocalResidual
{
    using ParentType = BaseLocalResidual;
    friend class StaggeredLocalResidual<TypeTag>;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using CellCenterResidual = CellCenterPrimaryVariables;

    static constexpr int turbulentKineticEnergyEqIdx = Indices::turbulentKineticEnergyEqIdx - ModelTraits::dim();
    static constexpr int dissipationEqIdx = Indices::dissipationEqIdx - ModelTraits::dim();

public:
    using ParentType::ParentType;

    // account for the offset of the cell center privars within the PrimaryVariables container
    static constexpr auto cellCenterOffset = ModelTraits::numEq() - CellCenterPrimaryVariables::dimension;
    static_assert(cellCenterOffset == ModelTraits::dim(), "cellCenterOffset must equal dim for staggered NavierStokes");

    //! Evaluate fluxes entering or leaving the cell center control volume.
    CellCenterPrimaryVariables computeStorageForCellCenter(const Problem& problem,
                                                           const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage = ParentType::computeStorageForCellCenter(problem, scv, volVars);

        storage[turbulentKineticEnergyEqIdx] = volVars.turbulentKineticEnergy();
        storage[dissipationEqIdx] = volVars.dissipation();

        return storage;
    }

    CellCenterPrimaryVariables computeSourceForCellCenter(const Problem& problem,
                                                          const Element &element,
                                                          const FVElementGeometry& fvGeometry,
                                                          const ElementVolumeVariables& elemVolVars,
                                                          const ElementFaceVariables& elemFaceVars,
                                                          const SubControlVolume &scv) const
    {
        CellCenterPrimaryVariables source = ParentType::computeSourceForCellCenter(problem, element, fvGeometry,
                                                                                   elemVolVars, elemFaceVars, scv);

        const auto& volVars = elemVolVars[scv];
        Scalar turbulentKineticEnergy = volVars.turbulentKineticEnergy();
        Scalar dissipation = volVars.dissipation();

        // production
        // turbulence production is equal to dissipation -> exclude both terms (according to local equilibrium hypothesis, see FLUENT)
        if (!volVars.isMatchingPoint())
        {
            source[turbulentKineticEnergyEqIdx] += 2.0 * volVars.kinematicEddyViscosity()
                                                   * volVars.stressTensorScalarProduct();
        }
        source[dissipationEqIdx] += volVars.cOneEpsilon()
                                    * dissipation / turbulentKineticEnergy
                                    * 2.0 * volVars.kinematicEddyViscosity()
                                    * volVars.stressTensorScalarProduct();

        // destruction
        // turbulence production is equal to dissipation -> exclude both terms (according to local equilibrium hypothesis, see FLUENT)
        if (!volVars.isMatchingPoint())
        {
            source[turbulentKineticEnergyEqIdx] -= dissipation;
        }
        source[dissipationEqIdx] -= volVars.cTwoEpsilon()
                                    * dissipation * dissipation
                                    / turbulentKineticEnergy;

        return source;
    }

protected:
     /*!
      * \brief Evaluate boundary conditions for a cell center dof
      */
    template<class ElementBoundaryTypes>
    void evalBoundaryForCellCenter_(CellCenterResidual& residual,
                                    const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const ElementFaceVariables& elemFaceVars,
                                    const ElementBoundaryTypes& elemBcTypes,
                                    const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        BaseLocalResidual::evalBoundaryForCellCenter_(residual, problem, element, fvGeometry,
                                                      elemVolVars, elemFaceVars, elemBcTypes, elemFluxVarsCache);
        for (auto&& scvf : scvfs(fvGeometry))
        {
            unsigned int elementID = problem.fvGridGeometry().elementMapper().index(element);
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideVolVars = elemVolVars[insideScv];

            // fixed value for the turbulent kinetic energy
            if (problem.inNearWallRegion(elementID))
            {
                residual[Indices::turbulentKineticEnergyEqIdx - cellCenterOffset]
                    = insideVolVars.turbulentKineticEnergy() - problem.turbulentKineticEnergyWallFunction(elementID);
            }

            // fixed value for the dissipation
            if (insideVolVars.inNearWallRegion() || insideVolVars.isMatchingPoint())
            {
                residual[Indices::dissipationEqIdx - cellCenterOffset]
                    = insideVolVars.dissipation() - problem.dissipationWallFunction(elementID);
            }
        }
    }
};
}

#endif
