// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SSTModel
 * \copydoc Dumux::SSTFluxVariablesImpl
 */
#ifndef DUMUX_SST_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_SST_STAGGERED_FLUXVARIABLES_HH

#include <numeric>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/navierstokes/fluxvariables.hh>
#include <dumux/freeflow/rans/twoeq/sst/fluxvariables.hh>

namespace Dumux {

/*!
  * \ingroup SSTModel
  * \brief The flux variables class for the SST model using the staggered grid discretization.
 */

// forward declaration
template<class TypeTag, class BaseFluxVariables, class DiscretizationMethod>
class SSTFluxVariablesImpl;

template<class TypeTag, class BaseFluxVariables>
class SSTFluxVariablesImpl<TypeTag, BaseFluxVariables, DiscretizationMethods::Staggered>
: public BaseFluxVariables
{
    using ParentType = BaseFluxVariables;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;

    static constexpr int turbulentKineticEnergyEqIdx = Indices::turbulentKineticEnergyEqIdx - ModelTraits::dim();
    static constexpr int dissipationEqIdx = Indices::dissipationEqIdx - ModelTraits::dim();

public:

    /*!
     * \brief Computes the flux for the cell center residual.
     */
    CellCenterPrimaryVariables computeMassFlux(const Problem& problem,
                                               const Element &element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               const SubControlVolumeFace &scvf,
                                               const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux = ParentType::computeMassFlux(problem, element, fvGeometry,
                                                                      elemVolVars, elemFaceVars, scvf, fluxVarsCache);

        // calculate advective flux
        auto upwindTermK = [](const auto& volVars)
        { return volVars.turbulentKineticEnergy() * volVars.density(); };
        auto upwindTermOmega = [](const auto& volVars)
        { return volVars.dissipation() * volVars.density(); };

        flux[turbulentKineticEnergyEqIdx]
            = ParentType::advectiveFluxForCellCenter(problem, fvGeometry, elemVolVars, elemFaceVars, scvf, upwindTermK);
        flux[dissipationEqIdx]
            = ParentType::advectiveFluxForCellCenter(problem, fvGeometry, elemVolVars, elemFaceVars, scvf, upwindTermOmega);

        // calculate diffusive flux
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        Scalar insideCoeff_k = 0.0, insideCoeff_w = 0.0, outsideCoeff_k = 0.0, outsideCoeff_w = 0.0;

        if(problem.sstModelVersion() == SSTModel::BSL)
        {
            insideCoeff_k = insideVolVars.viscosity()
                                    + ( insideVolVars.sigmaKBSL() * insideVolVars.dynamicEddyViscosity() );
            outsideCoeff_k = outsideVolVars.viscosity()
                                    + ( outsideVolVars.sigmaKBSL() * outsideVolVars.dynamicEddyViscosity() );
            insideCoeff_w = insideVolVars.viscosity()
                                    + ( insideVolVars.sigmaOmegaBSL() * insideVolVars.dynamicEddyViscosity() );
            outsideCoeff_w = outsideVolVars.viscosity()
                                    + ( outsideVolVars.sigmaOmegaBSL() * outsideVolVars.dynamicEddyViscosity() );
        }
        else if(problem.sstModelVersion() == SSTModel::SST)
        {
            insideCoeff_k = insideVolVars.viscosity()
                                    + ( insideVolVars.sigmaKSST() * insideVolVars.dynamicEddyViscosity() );
            outsideCoeff_k = outsideVolVars.viscosity()
                                    + ( outsideVolVars.sigmaKSST() * outsideVolVars.dynamicEddyViscosity() );
            insideCoeff_w = insideVolVars.viscosity()
                                    + ( insideVolVars.sigmaOmegaSST() * insideVolVars.dynamicEddyViscosity() );
            outsideCoeff_w = outsideVolVars.viscosity()
                                    + ( outsideVolVars.sigmaOmegaSST() * outsideVolVars.dynamicEddyViscosity() );
        }
        else
            DUNE_THROW(Dune::NotImplemented, "\nThis SST Model is not implemented.\n");


        // scale by extrusion factor
        insideCoeff_k *= insideVolVars.extrusionFactor();
        outsideCoeff_k *= outsideVolVars.extrusionFactor();
        insideCoeff_w *= insideVolVars.extrusionFactor();
        outsideCoeff_w *= outsideVolVars.extrusionFactor();

        Scalar distance = 0.0;
        Scalar coeff_k = 0.0;
        Scalar coeff_w = 0.0;
        if (scvf.boundary())
        {
            distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            coeff_k = insideCoeff_k;
            coeff_w = insideCoeff_w;
        }
        else
        {
            // average and distance
            coeff_k = arithmeticMean(insideCoeff_k, outsideCoeff_k,
                                    (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm(),
                                    (insideScv.dofPosition() - scvf.ipGlobal()).two_norm());
            coeff_w = arithmeticMean(insideCoeff_w, outsideCoeff_w,
                                    (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm(),
                                    (insideScv.dofPosition() - scvf.ipGlobal()).two_norm());
            distance = (outsideScv.dofPosition() - insideScv.dofPosition()).two_norm();
        }

        const auto bcTypes = problem.boundaryTypes(element, scvf);
        if (!(scvf.boundary() && (bcTypes.isOutflow(Indices::turbulentKineticEnergyEqIdx)
                                  || bcTypes.isSymmetry())))
        {
            flux[turbulentKineticEnergyEqIdx]
                += coeff_k / distance
                   * (insideVolVars.turbulentKineticEnergy() - outsideVolVars.turbulentKineticEnergy())
                   * Extrusion::area(fvGeometry, scvf);
        }
        if (!(scvf.boundary() && (bcTypes.isOutflow(Indices::dissipationEqIdx)
                                  || bcTypes.isSymmetry())))
        {
            flux[dissipationEqIdx]
                += coeff_w / distance
                   * (insideVolVars.dissipation() - outsideVolVars.dissipation())
                   * Extrusion::area(fvGeometry, scvf);
        }
        return flux;
    }

    /*!
     * \brief Returns the momentum flux over all staggered faces.
     */
    FacePrimaryVariables computeMomentumFlux(const Problem& problem,
                                             const Element& element,
                                             const SubControlVolumeFace& scvf,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementFaceVariables& elemFaceVars,
                                             const GridFluxVariablesCache& gridFluxVarsCache)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        return ParentType::computeFrontalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, gridFluxVarsCache)
               + ParentType::computeLateralMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, gridFluxVarsCache)
               + 2.0 / ModelTraits::dim() * insideVolVars.density() * insideVolVars.turbulentKineticEnergy()
                 * Extrusion::area(fvGeometry, scvf) * scvf.directionSign() * insideVolVars.extrusionFactor();
    }
};

} // end namespace

#endif
