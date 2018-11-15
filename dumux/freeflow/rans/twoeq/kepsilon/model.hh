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
 * \ingroup KEpsilonModel
 *
 * \brief A single-phase, isothermal k-epsilon model
 *
 * \copydoc RANSModel
 *
 * The k-epsilon models calculate the eddy viscosity with two additional PDEs,
 * one for the turbulent kinetic energy (k) and for the dissipation (\f$ \varepsilon \f$).
 * The model uses the one proposed by Launder and Sharma \cite Launder1994a.
 *
 * The turbulent kinetic energy balance is:
 * \f[
 *    \frac{\partial \left( k \right)}{\partial t}
 *    + \nabla \cdot \left( \textbf{v} k \right)
 *    - \nabla \cdot \left( \left( \nu + \frac{\nu_\text{t}}{\sigma_\text{k}} \right) \nabla k \right)
 *    - 2 \nu_\text{t} \textbf{S} \cdot \textbf{S}
 *    + \varepsilon
 *    = 0
 * \f].
 *
 * The dissipation balance is:
 * \f[
 *   \frac{\partial \left( \varepsilon \right)}{\partial t}
 *   + \nabla \cdot \left( \textbf{v} \varepsilon \right)
 *   - \nabla \cdot \left( \left( \nu + \frac{\nu_\text{t}}{\sigma_{\varepsilon}} \right) \nabla \varepsilon \right)
 *   - C_{1\varepsilon} \frac{\varepsilon}{k} 2 \nu_\text{t} \textbf{S} \cdot \textbf{S}
 *   + C_{2\varepsilon} \frac{\varepsilon^2}{k}
 *   = 0
 * \f].
 *
 * The kinematic eddy viscosity \f$ \nu_\text{t} \f$ is:
 * \f[
 * \nu_\text{t} = C_\mu \frac{k^2}{\tilde{\varepsilon}}
 * \f].
 *
 * Finally, the model is closed with the following constants:
 * \f[ \sigma_\text{k} = 1.00 \f]
 * \f[ \sigma_\varepsilon =1.30 \f]
 * \f[ C_{1\varepsilon} = 1.44 \f]
 * \f[ C_{2\varepsilon} = 1.92 \f]
 * \f[ C_\mu = 0.09 \f]
 */

#ifndef DUMUX_KEPSILON_MODEL_HH
#define DUMUX_KEPSILON_MODEL_HH

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
 * \ingroup KEpsilonModel
 * \brief Traits for the k-epsilon model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct KEpsilonModelTraits : RANSModelTraits<dimension>
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions,
    //! one mass balance equation and two turbulent transport equations
    static constexpr int numEq() { return dim()+1+2; }

    //! The number of components
    static constexpr int numComponents() { return 1; }

    //! the indices
    using Indices = KEpsilonIndices<dim(), numComponents()>;
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal k-epsilon model
///////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal k-epsilon model
struct KEpsilon { using InheritsFrom = std::tuple<RANS>; };
} // end namespace TTag

//!< states some specifics of the isothermal k-epsilon model
SET_PROP(KEpsilon, ModelTraits)
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::FVGridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = KEpsilonModelTraits<dim>;
};

//! The flux variables
SET_PROP(KEpsilon, FluxVariables)
{
private:
    using BaseFluxVariables = NavierStokesFluxVariables<TypeTag>;
public:
    using type = KEpsilonFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The local residual
SET_PROP(KEpsilon, LocalResidual)
{
private:
    using BaseLocalResidual = NavierStokesResidual<TypeTag>;
public:
    using type = KEpsilonResidual<TypeTag, BaseLocalResidual>;
};

//! Set the volume variables property
SET_PROP(KEpsilon, VolumeVariables)
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = KEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific I/O fields
SET_TYPE_PROP(KEpsilon, IOFields, KEpsilonIOFields);

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal k-epsilon model
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal k-epsilon model
struct KEpsilonNI { using InheritsFrom = std::tuple<RANSNI>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
SET_PROP(KEpsilonNI, ModelTraits)
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::FVGridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = KEpsilonModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
SET_PROP(KEpsilonNI, VolumeVariables)
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = KEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific non-isothermal I/O fields
SET_TYPE_PROP(KEpsilonNI, IOFields, FreeflowNonIsothermalIOFields<KEpsilonIOFields, true/*turbulenceModel*/>);

// \}
}

} // end namespace

#endif
