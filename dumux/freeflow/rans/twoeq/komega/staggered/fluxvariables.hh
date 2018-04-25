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
  * \ingroup KOmegaModel
  * \copydoc Dumux::KOmegaFluxVariablesImpl
  */
#ifndef DUMUX_KOMEGA_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_KOMEGA_STAGGERED_FLUXVARIABLES_HH

#include <numeric>
#include <dumux/common/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/navierstokes/fluxvariables.hh>
#include <dumux/freeflow/rans/twoeq/komega/fluxvariables.hh>
#include <dumux/freeflow/rans/twoeq/komega/models.hh>

namespace Dumux
{

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class KOmegaFluxVariablesImpl;

/*!
  * \ingroup KOmegaModel
  * \brief The flux variables class for the k-omega model using the staggered grid discretization.
 */
template<class TypeTag>
class KOmegaFluxVariablesImpl<TypeTag, DiscretizationMethod::staggered>
: public NavierStokesFluxVariables<TypeTag>
{
    using ParentType = NavierStokesFluxVariables<TypeTag>;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;
    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;
    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    enum {
        turbulentKineticEnergyEqIdx = Indices::turbulentKineticEnergyEqIdx,
        dissipationEqIdx = Indices::dissipationEqIdx,
    };

public:

    /*!
    * \brief Computes the flux for the cell center residual.
    */
    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const ElementFaceVariables& elemFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux = ParentType::computeFluxForCellCenter(problem, element, fvGeometry,
                                                                               elemVolVars, elemFaceVars, scvf, fluxVarsCache);

        // calculate advective flux
        const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
        const bool isOutflowK = scvf.boundary() && bcTypes.isOutflow(turbulentKineticEnergyEqIdx);
        const bool isOutflowOmega = scvf.boundary() && bcTypes.isOutflow(dissipationEqIdx);
        auto upwindTermK = [](const auto& volVars)
        {
            return volVars.turbulentKineticEnergy();
        };
        auto upwindTermOmega = [](const auto& volVars)
        {
            return volVars.dissipation();
        };

        flux[turbulentKineticEnergyEqIdx - ModelTraits::dim()]
            = ParentType::advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTermK, isOutflowK);
        flux[dissipationEqIdx - ModelTraits::dim()]
            = ParentType::advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTermOmega, isOutflowOmega);
        Dune::dverb << " k_adv " << ParentType::advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTermK, isOutflowK);
        Dune::dverb << " w_adv " << ParentType::advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTermOmega, isOutflowOmega);

        // calculate diffusive flux
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // effective diffusion coefficients
        Scalar insideCoeff_k = insideVolVars.kinematicViscosity()
                                + ( insideVolVars.sigmaK() * insideVolVars.turbulentKineticEnergy() / insideVolVars.dissipation() );
        Scalar outsideCoeff_k = outsideVolVars.kinematicViscosity()
                                + ( outsideVolVars.sigmaK() * outsideVolVars.turbulentKineticEnergy() / outsideVolVars.dissipation() );
        Scalar insideCoeff_w = insideVolVars.kinematicViscosity()
                                + ( insideVolVars.sigmaOmega() * insideVolVars.turbulentKineticEnergy() / insideVolVars.dissipation() );
        Scalar outsideCoeff_w = outsideVolVars.kinematicViscosity()
                                + ( outsideVolVars.sigmaOmega() * outsideVolVars.turbulentKineticEnergy() / outsideVolVars.dissipation() );

        // scale by extrusion factor
        insideCoeff_k *= insideVolVars.extrusionFactor();
        outsideCoeff_k *= outsideVolVars.extrusionFactor();
        insideCoeff_w *= insideVolVars.extrusionFactor();
        outsideCoeff_w *= outsideVolVars.extrusionFactor();

        // average and distance
        Scalar coeff_k = (insideCoeff_k - outsideCoeff_k) * 0.5;
        Scalar coeff_w = (insideCoeff_w - outsideCoeff_w) * 0.5;
        Scalar distance = 0.0;

        // adapt boundary handling
        if (scvf.boundary())
        {
            distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            coeff_k = insideCoeff_k;
            coeff_w = insideCoeff_w;
        }
        else
        {
            distance = (outsideScv.dofPosition() - insideScv.dofPosition()).two_norm();
        }

        if (!isOutflowK)
        {
            flux[turbulentKineticEnergyEqIdx - ModelTraits::dim()]
                += coeff_k / distance
                   * (insideVolVars.turbulentKineticEnergy() - outsideVolVars.turbulentKineticEnergy())
                   * scvf.area();
            Dune::dverb << " k_diff " << coeff_k / distance
                                                 * (insideVolVars.turbulentKineticEnergy() - outsideVolVars.turbulentKineticEnergy())
                                                 * scvf.area();
        }
        if (!isOutflowOmega)
        {
            flux[dissipationEqIdx - ModelTraits::dim()]
                += coeff_w / distance
                   * (insideVolVars.dissipation() - outsideVolVars.dissipation())
                   * scvf.area();
            Dune::dverb << " w_diff " << coeff_w / distance
                                                 * (insideVolVars.dissipation() - outsideVolVars.dissipation())
                                                 * scvf.area();
        }
        return flux;
    }
};

} // end namespace

#endif
