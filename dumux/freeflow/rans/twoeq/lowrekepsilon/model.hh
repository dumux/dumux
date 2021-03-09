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
 *    \frac{\partial \left( \varrho k \right)}{\partial t}
 *    + \nabla \cdot \left( \textbf{v} \varhho k \right)
 *    - \nabla \cdot \left( \left( \mu + \frac{\mu_\text{t}}{\sigma_\text{k}} \right) \nabla k \right)
 *    - 2 \mu_\text{t} \textbf{S} \cdot \textbf{S}
 *    + \varrho \tilde{\varepsilon}
 *    + D_\varepsilon \varrho
 *    = 0
 * \f].
 *
 * The dissipation balance is changed by introducing additional functions
 * (\f$ E_\text{k}\f$, \f$ f_1 \f$, and \f$ f_2 \f$) to account for a dampening towards the wall:
 * \f[
 *   \frac{\partial \left( \varrho \tilde{\varepsilon} \right)}{\partial t}
 *   + \nabla \cdot \left( \textbf{v} \varrho \tilde{\varepsilon} \right)
 *   - \nabla \cdot \left( \left( \mu + \frac{\mu_\text{t}}{\sigma_{\varepsilon}} \right) \nabla \tilde{\varepsilon} \right)
 *   - C_{1\tilde{\varepsilon}} f_1 \frac{\tilde{\varepsilon}}{k} 2 \mu_\text{t} \textbf{S} \cdot \textbf{S}
 *   + C_{2\tilde{\varepsilon}} \varrho f_2 \frac{\tilde{\varepsilon}^2}{k}
 *   - E_\text{k} \varrho
 *   = 0
 * \f].
 *
 * The kinematic eddy viscosity \f$ \nu_\text{t} \f$ is dampened by \f$ f_\mu \f$:
 * \f[
 * \mu_\text{t} = \varrho C_\mu f_\mu \frac{k^2}{\tilde{\varepsilon}}
 * \f].
 *
 * The auxiliary and dampening functions are defined as:
 * \f[ D_\varepsilon = 2 \nu \frac{k}{y^2} \f]
 * \f[ E_\text{k} = -2 \nu \frac{\tilde{\varepsilon}}{y^2} \exp \left( -0.5 y^+ \right) \f]
 * \f[ f_1 = 1 \f]
 * \f[ f_2 = 1 - 0.22 \exp \left( - \left( \frac{\mathit{Re}_\text{t}}{6} \right)^2 \right) \f]
 * \f[ f_\mu = 1 - \exp \left( -0.0115 y^+ \right) \f]
 * \f[ \mathit{Re}_\text{t} = \frac{k^2}{\nu \tilde{\varepsilon}} \f].
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
#include <dumux/freeflow/rans/twoeq/indices.hh>
#include <dumux/freeflow/turbulencemodel.hh>

#include "problem.hh"
#include "fluxvariables.hh"
#include "localresidual.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux {
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
    static constexpr int numFluidComponents() { return 1; }

    //! the indices
    using Indices = RANSTwoEqIndices<dim(), numFluidComponents()>;

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::lowrekepsilon; }
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal low-Reynolds k-epsilon model
///////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal low-Reynolds k-epsilon model
struct LowReKEpsilon { using InheritsFrom = std::tuple<RANS>; };
} // end namespace TTag

//!< states some specifics of the isothermal low-Reynolds k-epsilon model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::LowReKEpsilon>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = LowReKEpsilonModelTraits<dim>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::LowReKEpsilon>
{
private:
    using BaseFluxVariables = NavierStokesFluxVariables<TypeTag>;
public:
    using type = LowReKEpsilonFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::LowReKEpsilon>
{
private:
    using BaseLocalResidual = NavierStokesResidual<TypeTag>;
public:
    using type = LowReKEpsilonResidual<TypeTag, BaseLocalResidual>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::LowReKEpsilon>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = LowReKEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::LowReKEpsilon> { using type = LowReKEpsilonIOFields; };

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal low-Reynolds k-epsilon model
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, non-isothermal low-Reynolds k-epsilon model
struct LowReKEpsilonNI { using InheritsFrom = std::tuple<LowReKEpsilon, RANSNI>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::LowReKEpsilonNI>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = LowReKEpsilonModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::LowReKEpsilonNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = LowReKEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific non-isothermal I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::LowReKEpsilonNI> { using type = FreeflowNonIsothermalIOFields<LowReKEpsilonIOFields, true/*turbulenceModel*/>; };

} // end properties
} // end namespace

#endif
