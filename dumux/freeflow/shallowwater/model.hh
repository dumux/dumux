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
 * \ingroup ShallowWaterModel
 *
 * \brief A two-dimensional shallow water equations model
 * The two-dimensonal shallow water equations (SWEs) can be written as
 *
 * \f[
 * \frac{\partial \mathbf{U}}{\partial t} +
 * \frac{\partial \mathbf{F}}{\partial x} + \\
 * \frac{\partial \mathbf{G}}{\partial y} - \mathbf{S_b} - \mathbf{S_f} = 0
 * \f]
 *
 * with  <B>U</B>, <B>F</B>, <B>G</B> are defined as
 *
 * \f[
 * \mathbf{U} = \begin{bmatrix} h \\ uh \\ vh \end{bmatrix},
 * \mathbf{F} = \begin{bmatrix} hu \\ hu^2  + \frac{1}{2} gh^2 \\ huv \end{bmatrix},
 * \mathbf{G} = \begin{bmatrix} hv \\ huv \\ hv^2  + \frac{1}{2} gh^2 \end{bmatrix}
 * \f]
 *
 * Z is the bedSurface, h the water depth, u the velocity in
 * x-direction and v the velocity in y-direction, g is the constant of gravity.
 *
 * The source terms for the bed friction <B>S_b</B> and bed slope
 * <B>S_f</B> are given as
 * \f[
 * \mathbf{S_b} = \begin{bmatrix} 0 \\ -gh \frac{\partial z}{\partial x}
 *                \\ -gh \frac{\partial z}{\partial y}\end{bmatrix},
 * \mathbf{S_f} = \begin{bmatrix} 0 \\ -ghS_{fx} \\ -ghS_{fy}\end{bmatrix}.
 * \f]
 *
 * A cell-centered finte volume method (cctpfa) is applied to solve the SWEs
 * in combination with a fully-implicit time discretization. For large time
 * step sizes (CFL > 1) this can lead to a strong semearing of sharp fronts.
 * This can be seen in the movement of short traveling waves (e.g. dam break
 * waves). Nevertheless the fully implicit time discretization showes often
 * no drawbacks in cases where no short waves are considered. Thus the model
 * can be a good choice for simulating flow in rivers and channels, where the
 * fully-implicit discretization allows large time steps and reduces the
 * overall computation time drastically.
 *
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_MODEL_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/flux/shallowwaterflux.hh>
#include <dumux/flux/fluxvariablescaching.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "indices.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup ShallowWaterModel
 * \brief Specifies a number properties of shallow water models.
 */
struct ShallowWaterModelTraits
{
    using Indices = ShallowWaterIndices;

    static constexpr int numEq() { return 3; }
    static constexpr int numPhases() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! Enable diffusion
    static constexpr bool enableDiffusion() { return false; }
};

/*!
 * \ingroup ShallowWaterModel
 * \brief Traits class for the volume variables of the shallow water model.
 *
 * \tparam PV The type used for primary variables
 * \tparam MT The model traits
 */
template<class PV,
         class MT>
struct ShallowWaterVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};


namespace Properties {

//! Type tag for shallow water equation model inherits from model properties
namespace TTag {
struct ShallowWater { using InheritsFrom = std::tuple<ModelProperties>; };
}// end namespace TTag

//////////////////////////////////////////////////////////////////
// Define properties
//////////////////////////////////////////////////////////////////

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterResidual<TypeTag>; };

template<class TypeTag>
struct FluxVariables<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterFluxVariables<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ShallowWater>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = ShallowWaterVolumeVariablesTraits<PV, MT>;
public:
    using type = ShallowWaterVolumeVariables<Traits>;
};

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::ShallowWater>
{ using type = FluxVariablesCaching::EmptyCache< GetPropType<TypeTag, Properties::Scalar> >; };

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::ShallowWater>
{ using type = FluxVariablesCaching::EmptyCacheFiller; };

template<class TypeTag>
struct IOFields<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterIOFields; };

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::ShallowWater>
{ using type = ShallowWaterFlux< GetPropType<TypeTag, Properties::NumEqVector> >; };

//template<class TypeTag> struct DiffusionType<TypeTag, TTag::ShallowWater> {using type = ShallowWaterExactRiemannSolver<TypeTag>;};

} // end properties
} // end namespace Dumux

#endif
