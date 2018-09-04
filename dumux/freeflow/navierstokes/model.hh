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
 * \ingroup NavierStokesModel
 *
 * \brief A single-phase, isothermal Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\textup{T}}) = \nabla \cdot (\mu (\nabla \textbf{v} + \nabla \textbf{v}^{\textup{T}}))
   - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * By setting the property <code>EnableInertiaTerms</code> to <code>false</code> the Stokes
 * equation can be solved. In this case the term
 * \f[
 *    \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\textup{T}})
 * \f]
 * is neglected.
 *
 * The <B> mass balance equation </B>
 * \f[
       \frac{\partial \varrho}{\partial t} + \nabla \cdot (\varrho \textbf{v}) - q = 0
 * \f]
 *
 * closes the system.
 *
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_NAVIERSTOKES_MODEL_HH
#define DUMUX_NAVIERSTOKES_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/nonisothermal/model.hh>
#include <dumux/freeflow/nonisothermal/indices.hh>
#include <dumux/freeflow/nonisothermal/vtkoutputfields.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "fluxvariablescache.hh"
#include "indices.hh"
#include "vtkoutputfields.hh"

#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fourierslaw.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Traits for the Navier-Stokes model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct NavierStokesModelTraits
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions
    //! and one mass balance equation.
    static constexpr int numEq() { return dimension+1; }

    //! The number of phases is 1
    static constexpr int numPhases() { return 1; }

    //! The number of components is 1
    static constexpr int numComponents() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! The one-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return false; }

    //! The model is isothermal
    static constexpr bool enableEnergyBalance() { return false; }

    //! The model does not include a turbulence model
    static constexpr bool usesTurbulenceModel() { return false; }

    //! The type of turbulence model used (none)
    static constexpr auto turbulenceModel = TurbulenceModel::none;

    //! the indices
    using Indices = NavierStokesIndices<dim()>;

    //! return the names of the primary variables in cells
    template <class FluidSystem = void>
    static std::string primaryVariableNameCell(int pvIdx = 0, int state = 0)
    {
        return "p";
    }

    //! return the names of the primary variables on faces
    template <class FluidSystem = void>
    static std::string primaryVariableNameFace(int pvIdx = 0, int state = 0)
    {
        return "v";
    }
};

/*!
 * \ingroup NavierStokesModel
 * \brief Traits class for the volume variables of the Navier-Stokes model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam MT The model traits
 */
template<class PV,
         class FSY,
         class FST,
         class MT>
struct NavierStokesVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using ModelTraits = MT;
};

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Navier-Stokes model
NEW_TYPE_TAG(NavierStokes, INHERITS_FROM(FreeFlow));

//! The type tag for the corresponding non-isothermal model
NEW_TYPE_TAG(NavierStokesNI, INHERITS_FROM(NavierStokes));

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
SET_BOOL_PROP(NavierStokes, NormalizePressure, true); //!< Normalize the pressure term in the momentum balance by default

//!< states some specifics of the Navier-Stokes model
SET_PROP(NavierStokes, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr auto dim = GridView::dimension;
public:
    using type = NavierStokesModelTraits<dim>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(NavierStokes, FluidState){
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = Dumux::ImmiscibleFluidState<Scalar, FluidSystem>;
};

//! The local residual
SET_TYPE_PROP(NavierStokes, LocalResidual, NavierStokesResidual<TypeTag>);

//! Set the volume variables property
SET_PROP(NavierStokes, VolumeVariables)
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
public:
    using type = NavierStokesVolumeVariables<Traits>;
};

//! The flux variables
SET_TYPE_PROP(NavierStokes, FluxVariables, NavierStokesFluxVariables<TypeTag>);

//! The flux variables cache class, by default the one for free flow
SET_TYPE_PROP(NavierStokes, FluxVariablesCache, FreeFlowFluxVariablesCache<TypeTag>);

//! The specific vtk output fields
SET_PROP(NavierStokes, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
public:
     using type = NavierStokesVtkOutputFields<FVGridGeometry>;
};
//////////////////////////////////////////////////////////////////
// Property values for non-isothermal Navier-Stokes model
//////////////////////////////////////////////////////////////////

//! The model traits of the non-isothermal model
SET_PROP(NavierStokesNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr auto dim = GridView::dimension;
    using IsothermalTraits = NavierStokesModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! The specific non-isothermal vtk output fields
SET_PROP(NavierStokesNI, VtkOutputFields)
{
private:
     using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
     using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
     using IsothermalFields = NavierStokesVtkOutputFields<FVGridGeometry>;
public:
     using type = FreeflowNonIsothermalVtkOutputFields<IsothermalFields, ModelTraits>;
};

 // \}
}

} // end namespace

#endif // DUMUX_NAVIERSTOKES_MODEL_HH
