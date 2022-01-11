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
 * \ingroup TwoEqRansModel
 *
 * \brief A single-phase, isothermal 2-Eq. model
 *
 * \copydoc RANSModel
 *
 * Two additional PDEs, the turbulentKineticEnergy (k), and either the dissipation freqency (omega), or the dissipation (epsilon),
 * are used to calculate the eddy viscosity for this model.
 *
 * Turbulent Kinetic Energy balance:
 * \f[
 * \frac{\partial \left( \varrho k \right)}{\partial t}
 * + \nabla \cdot \left( \mathbf{v} \varrho k \right)
 * - \nabla \cdot \left[ \left( \mu +  \sigma_\textrm{k} \mu_\textrm{t} \right) \nabla k \right]
 * - Prod
 * + Dest
 * = 0
 * \f]
 *
 * Dissipation balance:
 * \f[
 * \frac{\partial \left( \varrho \omega,\varepsilon \right)}{\partial t}
 * + \nabla \cdot \left( \mathbf{v} \varrho \omega,\varepsilon \right)
 * - \nabla \cdot \left[ \left( \mu + \sigma_{\omega,\varepsilon} \mu_\textrm{t} \right) \nabla \omega,\varepsilon \right]
 * - Prod
 * + Dest
 * - CrossDisp
 * = 0
 * \f]
 *
 * The dynamic eddy viscosity \f$ \mu_\textrm{t} \f$ is calculated as follows:
 * \f[ \mu_\textrm{t} = \varrho \frac{k}{\tilde{\omega}} \f]
 * or
 * \f[ \mu_\textrm{t} = \varrho C_\mu \frac{k^2}{\tilde{\varepsilon}} \f].
 */
#ifndef DUMUX_TWOEQ_MODEL_HH
#define DUMUX_TWOEQ_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/turbulence/spatialparams.hh>

#include "fluxvariables.hh"
#include "localresidual.hh"
#include "iofields.hh"
#include "indices.hh"
#include "volumevariables.hh"

namespace Dumux::Properties {

/*!
 *\ingroup TwoEqModel
 * \brief Traits for the Two Eq Rans model
 *
 * \tparam dimension The dimension of the problem
 */
template<class BaseMassTraits>
struct TwoEqModelTraits : public BaseMassTraits
{
    using MassTraits = BaseMassTraits;

    //! There is one mass balance equation and two turbulent transport equations
    static constexpr int numEq() { return 1 + numTurbulenceEqs(); }

    //! The number of turbulence variables
    static constexpr int numTurbulenceEqs() { return 2; }

    //! The model does not include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }

    //! The indices
    using Indices = TwoEqIndices<MassTraits::numFluidComponents()>;

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel() { return TurbulenceModel::twoeq; }
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal twoeq rans model
///////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal two Eq model
struct RANSTwoEq { using InheritsFrom = std::tuple<NavierStokesMassOneP>; };
} // end namespace TTag

//! states some specifics of the isothermal two Eq model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::RANSTwoEq>
{
    using BaseMassTraits = NavierStokesMassOnePModelTraits;
    using type = TwoEqModelTraits<BaseMassTraits>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::RANSTwoEq>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using BaseFluxVars = FluxVariablesBase<Problem, FVElementGeometry, ElementVolumeVariables, ElementFluxVariablesCache>;
public:
    using type = TwoEqFluxVariables<TypeTag, BaseFluxVars>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::RANSTwoEq>
{ using type = TwoEqResidual<TypeTag, NavierStokesMassOnePLocalResidual<TypeTag>>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::RANSTwoEq>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesMassOnePVolumeVariablesTraits<PV, FSY, FST, MT>;
    using BaseNavierStokesMassOnePVolumeVariables = NavierStokesMassOnePVolumeVariables<Traits>;
    static constexpr int dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension;
public:
    using type = TwoEqVolumeVariables<Traits, BaseNavierStokesMassOnePVolumeVariables, dim>;
};

//! The spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RANSTwoEq>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = RANSSpatialParams<GridGeometry>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::RANSTwoEq> { using type = TwoEqIOFields; };

///////////////////////////////////////////////////////////////////////////
// default property values for the non-isothermal two Eq single phase model
///////////////////////////////////////////////////////////////////////////

// TODO: Set up non-isothermal models


} // Dumux::Properties

#endif
