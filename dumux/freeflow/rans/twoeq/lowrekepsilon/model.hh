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
 * \ingroup LowReKEpsilonModel
 *
 * \brief A single-phase, isothermal low-Reynolds k-epsilon model
 *
 * \copydoc RANSModel
 *
 * The low-Reynolds k-epsilon models calculate the eddy viscosity with two additional PDEs,
 * one for the turbulent kinetic energy (k) and for the dissipation (\f$ \varepsilon \f$).
 * The model uses the one proposed by Chien \cite Chien1982a.
 * A good overview and additional models are given in Patel et al. \cite Patel1985a.
 *
 * The turbulent kinetic energy balance is identical with the one from the k-epsilon model,
 * but the dissipation includes a dampening function (\f$ D_\varepsilon \f$):
 * \f$ \varepsilon = \tilde{\varepsilon} + D_\varepsilon \f$:
 *
 * \f[
 *    \frac{\partial \left( k \right)}{\partial t}
 *    + \nabla \cdot \left( \textbf{v} k \right)
 *    - \nabla \cdot \left( \left( \nu + \frac{\nu_\text{t}}{\sigma_\text{k}} \right) \nabla k \right)
 *    - 2 \nu_\text{t} \textbf{S} \cdot \textbf{S}
 *    + \tilde{\varepsilon}
 *    + D_\varepsilon
 *    = 0
 * \f].
 *
 * The dissipation balance is changed by introducing additional functions
 * (\f$ E_\text{k}\f$, \f$ f_1 \f$, and \f$ f_2 \f$) to account for a dampening towards the wall:
 * \f[
 *   \frac{\partial \left( \tilde{\varepsilon} \right)}{\partial t}
 *   + \nabla \cdot \left( \textbf{v} \tilde{\varepsilon} \right)
 *   - \nabla \cdot \left( \left( \nu + \frac{\nu_\text{t}}{\sigma_{\varepsilon}} \right) \nabla \tilde{\varepsilon} \right)
 *   - C_{1\tilde{\varepsilon}} f_1 \frac{\tilde{\varepsilon}}{k} 2 \nu_\text{t} \textbf{S} \cdot \textbf{S}
 *   + C_{2\tilde{\varepsilon}} f_2 \frac{\tilde{\varepsilon}^2}{k}
 *   - E_\text{k}
 *   = 0
 * \f].
 *
 * The kinematic eddy viscosity \f$ \nu_\text{t} \f$ is dampened by \f$ f_\mu \f$:
 * \f[
 * \nu_\text{t} = C_\mu f_\mu \frac{k^2}{\tilde{\varepsilon}}
 * \f].
 *
 * The auxiliary and dampening functions are defined as:
 * \f[ D_\varepsilon = 2 \nu \nicefrac{k}{y^2} \f]
 * \f[ E_\text{k} = -2 \nu \frac{\tilde{\varepsilon}}{y^2} \exp \left( -0.5 y^+ \right) \f]
 * \f[ f_1 = 1 \f]
 * \f[ f_2 = 1 - 0.22 \exp \left( - \left( \frac{\mathit{Re}_\text{t}}{6} \right)^2 \right) \f]
 * \f[ f_\mu = 1 - \exp \left( -0.0115 y^+ \right) \f]
 * \f[ \mathit{Re}_\text{t} = \frac{k^2}{\nu \tilde{\varepsilon}} \f]
 * .
 *
 * Finally, the model is closed with the following constants:
 * \f[ \sigma_\text{k} = 1.00 \f]
 * \f[ \sigma_\varepsilon =1.30 \f]
 * \f[ C_{1\tilde{\varepsilon}} = 1.35 \f]
 * \f[ C_{2\tilde{\varepsilon}} = 1.80 \f]
 * \f[ C_\mu = 0.09 \f]
 */

#ifndef DUMUX_LOWREKEPSILON_MODEL_HH
#define DUMUX_LOWREKEPSILON_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>

#include "fluxvariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux
{
namespace Properties {

/*!
 * \ingroup LowReKEpsilonModel
 * \brief Traits for the low-Reynolds k-epsilon model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct LowReKEpsilonModelTraits : RANSModelTraits<dimension>
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions,
    //! one mass balance equation and two turbulent transport equations
    static constexpr int numEq() { return dim()+1+2; }

    //! The number of components
    static constexpr int numComponents() { return 1; }

    //! the indices
    using Indices = LowReKEpsilonIndices<dim(), numComponents()>;
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal low-Reynolds k-epsilon model
///////////////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal low-Reynolds k-epsilon model
NEW_TYPE_TAG(LowReKEpsilon, INHERITS_FROM(RANS));

//!< states some specifics of the isothermal low-Reynolds k-epsilon model
SET_PROP(LowReKEpsilon, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = LowReKEpsilonModelTraits<dim>;
};

//! The flux variables
SET_PROP(LowReKEpsilon, FluxVariables)
{
private:
    using BaseFluxVariables = NavierStokesFluxVariables<TypeTag>;
public:
    using type = LowReKEpsilonFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The local residual
SET_PROP(LowReKEpsilon, LocalResidual)
{
private:
    using BaseLocalResidual = NavierStokesResidual<TypeTag>;
public:
    using type = LowReKEpsilonResidual<TypeTag, BaseLocalResidual>;
};

//! Set the volume variables property
SET_PROP(LowReKEpsilon, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = LowReKEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific I/O fields
SET_PROP(LowReKEpsilon, IOFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
public:
    using type = LowReKEpsilonIOFields<FVGridGeometry>;
};

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal low-Reynolds k-epsilon model
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal low-Reynolds k-epsilon model
NEW_TYPE_TAG(LowReKEpsilonNI, INHERITS_FROM(RANSNI));

//! The model traits of the non-isothermal model
SET_PROP(LowReKEpsilonNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = LowReKEpsilonModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
SET_PROP(LowReKEpsilonNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = LowReKEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific non-isothermal I/O fields
SET_PROP(LowReKEpsilonNI, IOFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using IsothermalFields = LowReKEpsilonIOFields<FVGridGeometry>;
public:
    using type = FreeflowNonIsothermalIOFields<IsothermalFields, ModelTraits>;
};

// \}
}

} // end namespace

#endif
