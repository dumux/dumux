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
 * \ingroup OneEqModel
 *
 * \brief A single-phase, isothermal one-equation turbulence model by Spalart-Allmaras
 *
 * \copydoc RANSModel
 *
 * This model, published by Spalart and Allmaras 1992 \cite Spalart1992a,
 * uses one additional PDE for a working variable \f$ \tilde{\nu} \f$.
 * This variable has the units of a viscosity and can be converted to the eddy
 * viscosity via a model function~(\f$ f_\text{v1} \f$):
 * \f[
 *  \nu_\text{t} = \tilde{\nu} f_\text{v1}
 * \f]
 *
 * Here, as proposed by Wilcox \cite Wilcox2008a and Versteeg \cite Versteeg2009a, the correction
 * term which account for the transition or trip, is dropped from the original equations,
 * such that the balance equation simplifies to:
 * \f[
 *   \frac{\partial \tilde{\nu}}{\partial t}
 *   + \nabla \cdot \left( \tilde{\nu} \textbf{v} \right)
 *   - c_\text{b1} \tilde{S} \tilde{\nu}
 *   - \frac{1}{\sigma_{\tilde{\nu}}} \nabla \cdot \left( \left[ \nu + \tilde{\nu} \right] \nabla \tilde{\nu} \right)
 *   - \frac{c_\text{b2}}{\sigma_{\tilde{\nu}}} \left| \nabla \tilde{\nu} \right|^2
 *   + c_\text{w1} f_\text{w} \frac{\tilde{\nu}^2}{y^2}
 *   = 0
 * \f]
 *
 * Here, a modified mean effective strain rate (\f$ \tilde{S} \f$) based on
 * the mean rotation rate tensor (\f$ \mathbf{\Omega} \f$) is used:
 * \f[
 *  \tilde{S} = \sqrt{2 \mathbf{\Omega} \cdot \mathbf{\Omega}} + \frac{\tilde{\nu}}{\kappa^2 y^2} f_\text{v2}
 * \f]
 * \f[
 *  \mathbf{\Omega} = \frac{1}{2} \left( \nabla \textbf{v}_\text{g}
 *                                  - \nabla \textbf{v}_\text{g}^\intercal \right)
 * \f]
 *
 * This balance equation is linked to the flow geometry by the distance to the closest wall ($y$).
 * Further, the model uses the following functions and expressions:
 * \f[ \chi = \frac{\tilde{\nu}}{\nu} \f]
 * \f[ f_\text{v1} = \frac{\chi^3}{\chi^3+c_\text{v1}^3} \f]
 * \f[ f_\text{v2} = 1 - \frac{\chi}{1+f_\text{v1}\chi} \f]
 * \f[ f_\text{w} = g_\text{w} \left( \frac{1+c_\text{w3}^6}{g^6_\text{w}+c_\text{w3}^6}
 *                             \right)^\frac{1}{6} \f]
 * \f[ g_\text{w} = r_\text{w} + c_\text{w2} (r_\text{w}^6 - r_\text{w}) \f]
 * \f[ r_\text{w} = \min \left[ \frac{\tilde{\nu}}{\tilde{S}\kappa^2 y^2},10\right] \f]
 * \f[ \sigma_{\tilde{\nu}} = \nicefrac{2}{3} \f]
 * \f[ c_\text{b1} = 0.1355 \f]
 * \f[ c_\text{b2} = 0.622 \f]
 * \f[ c_\text{v1} = 7.1 \f]
 * \f[ c_\text{w1} = \frac{c_\text{b1}}{\kappa^2}
 *                   + \frac{1+c_\text{b2}}{\sigma_{\tilde{\nu}}} \f]
 * \f[ c_\text{w2} = 0.3 \f]
 * \f[ c_\text{w3} = 2 \f]
 */

#ifndef DUMUX_ONEEQ_MODEL_HH
#define DUMUX_ONEEQ_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>
#include <dumux/freeflow/nonisothermal/iofields.hh>

#include "fluxvariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux
{
namespace Properties {

/*!
 * \ingroup OneEqModel
 * \brief Traits for the Spalart-Allmaras model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct OneEqModelTraits : RANSModelTraits<dimension>
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions,
    //! one mass balance equation and one turbulent transport equation
    static constexpr int numEq() { return dim()+1+1; }

    //! The number of components
    static constexpr int numComponents() { return 1; }

    //! the indices
    using Indices = OneEqIndices<dim(), numComponents()>;
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal Spalart-Allmaras model
///////////////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Spalart-Allmaras model
NEW_TYPE_TAG(OneEq, INHERITS_FROM(RANS));

//!< states some specifics of the isothermal Spalart-Allmaras model
SET_PROP(OneEq, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = OneEqModelTraits<dim>;
};

//! The flux variables
SET_PROP(OneEq, FluxVariables)
{
private:
    using BaseFluxVariables = NavierStokesFluxVariables<TypeTag>;
public:
    using type = OneEqFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The local residual
SET_PROP(OneEq, LocalResidual)
{
private:
    using BaseLocalResidual = NavierStokesResidual<TypeTag>;
public:
    using type = OneEqResidual<TypeTag, BaseLocalResidual>;
};

//! Set the volume variables property
SET_PROP(OneEq, VolumeVariables)
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
    using type = OneEqVolumeVariables<Traits, NSVolVars>;
};

//! The specific I/O fields
SET_TYPE_PROP(OneEq, IOFields, OneEqIOFields);

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal Spalart-Allmaras model
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Spalart-Allmaras model
NEW_TYPE_TAG(OneEqNI, INHERITS_FROM(RANSNI));

//! The model traits of the non-isothermal model
SET_PROP(OneEqNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = OneEqModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
SET_PROP(OneEqNI, VolumeVariables)
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
    using type = OneEqVolumeVariables<Traits, NSVolVars>;
};

//! The specific non-isothermal I/O fields
SET_TYPE_PROP(OneEqNI, IOFields, FreeflowNonIsothermalIOFields<OneEqIOFields, true/*turbulenceModel*/>);

// \}
}

} // end namespace

#endif
