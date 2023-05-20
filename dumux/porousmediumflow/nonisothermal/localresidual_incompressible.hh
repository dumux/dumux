// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NIModel
 * \brief Element-wise calculation of the local residual for non-isothermal
 *        fully implicit models.
 */

#ifndef DUMUX_ENERGY_LOCAL_RESIDUAL_INCOMPRESSIBLE_HH
#define DUMUX_ENERGY_LOCAL_RESIDUAL_INCOMPRESSIBLE_HH

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include "localresidual.hh"

namespace Dumux {

// forward declaration
template<class TypeTag, bool enableEneryBalance>
class EnergyLocalResidualIncompressibleImplementation;

template<class TypeTag>
using EnergyLocalResidualIncompressible = EnergyLocalResidualIncompressibleImplementation<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance()>;

/*!
 * \ingroup NIModel
 * \brief Element-wise calculation of the energy residual for non-isothermal problems.
 */
template<class TypeTag>
class EnergyLocalResidualIncompressibleImplementation<TypeTag, false>: public EnergyLocalResidualImplementation<TypeTag, false>
{}; //energy balance not enabled

/*!
 * \ingroup NIModel
 * \brief TODO docme!
 */
template<class TypeTag>
class EnergyLocalResidualIncompressibleImplementation<TypeTag, true>: public EnergyLocalResidualImplementation<TypeTag, true>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    enum { energyEqIdx = Indices::energyEqIdx };

public:
    /*!
     * \brief The advective phase energy fluxes for incompressible flow.
     *
     * Using specific internal energy $u$ instead of specific enthalpy \f$h\f$ for incompressible flow in convective flux
     * to account for otherwise neglected pressure work term (\f$\nabla p \cdot v\f$).
     *
     * Compressible formulation in EnergyLocalResidual (neglecting pressure work term (\f$\nabla p \cdot v\f$))
     * is
     \f{align*}{
     \frac{\partial}{\partial t} (\rho u) &= -\nabla \cdot (\rho v h) + \nabla \cdot (\lambda \nabla T)
     \f}
     *
     * Incompressible energy formulation is
     \f{align*}{
     \frac{\partial}{\partial t} (\rho u) = -\nabla \cdot (\rho v u) + \nabla \cdot (\lambda \nabla T)
     \f}
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     * \param phaseIdx The phase index
     */
    static void heatConvectionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars,
                                   int phaseIdx)
    {
        // internal energy used instead of enthalpy for incompressible flow
        auto upwindTerm = [phaseIdx](const auto& volVars)
        { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx)*volVars.internalEnergy(phaseIdx); };

        flux[energyEqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
    }
};

} // end namespace Dumux

#endif
