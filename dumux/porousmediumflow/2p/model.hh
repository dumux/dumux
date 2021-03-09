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
 * \ingroup TwoPModel
 * \brief Adaption of the fully implicit scheme to the two-phase flow model.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multi-phase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPFormulation::pwsn</tt> or <tt>TwoPFormulation::pnsw</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */

#ifndef DUMUX_TWOP_MODEL_HH
#define DUMUX_TWOP_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/spatialparams/fv.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include "formulation.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"
#include "saturationreconstruction.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPModel
 * \brief Specifies a number properties of two-phase models.
 */
template<TwoPFormulation formulation>
struct TwoPModelTraits
{
    using Indices = TwoPIndices;

    static constexpr TwoPFormulation priVarFormulation()
    { return formulation; }

    static constexpr int numEq() { return 2; }
    static constexpr int numFluidPhases() { return 2; }
    static constexpr int numFluidComponents() { return 2; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
};

/*!
 * \ingroup TwoPModel
 * \brief Traits class for the two-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam SSY The solid system type
 * \tparam SST The solid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 * \tparam SR The class used for reconstruction of
 *            nonwetting phase saturations in scvs
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT, class SR>
struct TwoPVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
    using SaturationReconstruction = SR;
};

// necessary for models derived from 2p
class TwoPIOFields;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the isothermal two-phase model
struct TwoP { using InheritsFrom = std::tuple<PorousMediumFlow>; };

//! The type tag for the non-isothermal two-phase model
struct TwoPNI { using InheritsFrom = std::tuple<TwoP>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal two-phase model
///////////////////////////////////////////////////////////////////////////
 //!< Set the default formulation to pwsn
template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoP>
{ static constexpr auto value = TwoPFormulation::p0s1; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoP> { using type = ImmiscibleLocalResidual<TypeTag>; };         //!< Use the immiscible local residual operator for the 2p model

//! The base model traits class
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::TwoP> { using type = TwoPModelTraits<getPropValue<TypeTag, Properties::Formulation>()>; };
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoP> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; }; //!< default the actually used traits to the base traits

//! Set the vtk output fields specific to the twop model
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoP> { using type = TwoPIOFields; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;

    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;

    using Traits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;
public:
    using type = TwoPVolumeVariables<Traits>;
};

//! The two-phase model uses the immiscible fluid state
template<class TypeTag>
struct FluidState<TypeTag, TTag::TwoP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

////////////////////////////////////////////////////////
// properties for the non-isothermal two-phase model
////////////////////////////////////////////////////////

//! The non-isothermal model traits class
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPNI> { using type = PorousMediumFlowNIModelTraits<GetPropType<TypeTag, Properties::BaseModelTraits>>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    static constexpr auto DM = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;

    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };

public:
    using type = TwoPVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

//! Set the vtk output fields specific to the non-isothermal twop model
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoPNI> { using type = EnergyIOFields<TwoPIOFields>; };

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPNI>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
