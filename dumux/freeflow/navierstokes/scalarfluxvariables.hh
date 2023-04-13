// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_SCALAR_CONSERVATION_MODEL_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_SCALAR_CONSERVATION_MODEL_FLUXVARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/upwindscheme.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables base class for scalar quantities balanced in the Navier-Stokes model.
 */
template<class Problem,
         class ModelTraits,
         class FluxTypes,
         class ElementVolumeVariables,
         class ElementFluxVariablesCache,
         class UpwindScheme = UpwindScheme<typename ProblemTraits<Problem>::GridGeometry>>
class NavierStokesScalarConservationModelFluxVariables
: public FluxVariablesBase<Problem,
                           typename ProblemTraits<Problem>::GridGeometry::LocalView,
                           ElementVolumeVariables,
                           ElementFluxVariablesCache>
{
    using Scalar = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables::value_type;
    using GridGeometry = typename ProblemTraits<Problem>::GridGeometry;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Extrusion = Extrusion_t<typename ProblemTraits<Problem>::GridGeometry>;

public:

    struct AdvectionType
    {
        static Scalar flux(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           const int phaseIdx,
                           const ElementFluxVariablesCache& elemFluxVarsCache)
        {
            const auto velocity = problem.faceVelocity(element, fvGeometry, scvf);
            const Scalar volumeFlux = velocity*scvf.unitOuterNormal()*Extrusion::area(fvGeometry, scvf)*extrusionFactor_(elemVolVars, scvf);
            return volumeFlux;
        }
    };

    /*!
     * \brief Returns the advective flux computed by the respective law.
     */
    template<typename FunctionType>
    Scalar getAdvectiveFlux(const FunctionType& upwindTerm) const
    {
        if constexpr (ModelTraits::enableAdvection())
        {
            const auto& scvf = this->scvFace();
            const auto velocity = this->problem().faceVelocity(this->element(), this->fvGeometry(), scvf);
            const Scalar volumeFlux = velocity*scvf.unitOuterNormal()*Extrusion::area(this->fvGeometry(), scvf)*extrusionFactor_(this->elemVolVars(), scvf);
            return UpwindScheme::apply(*this, upwindTerm, volumeFlux, 0/*phaseIdx*/);
        }
        else
            return 0.0;
    }

    /*!
     * \brief Returns the conductive energy flux computed by the respective law.
     */
    Scalar heatConductionFlux() const
    {
        if constexpr (ModelTraits::enableEnergyBalance())
        {
            using HeatConductionType = typename FluxTypes::HeatConductionType;
            return HeatConductionType::flux(this->problem(),
                                            this->element(),
                                            this->fvGeometry(),
                                            this->elemVolVars(),
                                            this->scvFace(),
                                            this->elemFluxVarsCache());
        }
        else
            return 0.0;
    }

    /*!
     * \brief Returns the advective energy flux.
     */
    Scalar heatAdvectionFlux() const
    {
        if constexpr (ModelTraits::enableEnergyBalance())
        {
            const auto upwindTerm = [](const auto& volVars) { return volVars.density() * volVars.enthalpy(); };
            return getAdvectiveFlux(upwindTerm);
        }
        else
            return 0.0;
    }

    /*!
     * \brief Returns the total energy flux.
     */
    Scalar heatFlux() const
    {
        return heatConductionFlux() + heatAdvectionFlux();
    }

    /*!
     * \brief Adds the energy flux to a given flux vector.
     */
    template<class NumEqVector>
    void addHeatFlux(NumEqVector& flux) const
    {
        if constexpr (ModelTraits::enableEnergyBalance())
            flux[ModelTraits::Indices::energyEqIdx] = heatFlux();
    }

private:

    static Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
    }
};

} // end namespace Dumux

#endif
