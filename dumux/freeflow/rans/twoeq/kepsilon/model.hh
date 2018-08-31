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

    //! return the names of the primary variables in cells
    template<class FluidSystem = void>
    static std::string primaryVariableNameCell(int pvIdx, int state = 0)
    {
        using ParentType = RANSModelTraits<dimension>;
        switch (pvIdx) {
            case 0:
                return ParentType::template primaryVariableNameCell<FluidSystem>(pvIdx, state);
            case 1:
                return "k";
            default:
                return "epsilon";
        }
    }
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal k-epsilon model
///////////////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal k-epsilon model
NEW_TYPE_TAG(KEpsilon, INHERITS_FROM(RANS));

//!< states some specifics of the isothermal k-epsilon model
SET_PROP(KEpsilon, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
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
    using type = KEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific I/O fields
SET_PROP(KEpsilon, IOFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
public:
    using type = KEpsilonIOFields<FVGridGeometry>;
};

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal k-epsilon model
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal k-epsilon model
NEW_TYPE_TAG(KEpsilonNI, INHERITS_FROM(RANSNI));

//! The model traits of the non-isothermal model
SET_PROP(KEpsilonNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = KEpsilonModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
SET_PROP(KEpsilonNI, VolumeVariables)
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
    using type = KEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific non-isothermal I/O fields
SET_PROP(KEpsilonNI, IOFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using IsothermalFields = KEpsilonIOFields<FVGridGeometry>;
public:
    using type = FreeflowNonIsothermalIOFields<IsothermalFields, ModelTraits>;
};

// \}
}

} // end namespace

#endif
