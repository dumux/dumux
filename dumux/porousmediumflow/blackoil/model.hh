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
 * \ingroup BlackOilModel
 * \brief A fully-implicit black-oil flow model using the VCVF discretization.
 *
 * The black-oil model is a three-phase, three-component model widely
 * used for oil reservoir simulation.  The phases are denoted by lower
 * index \f$\alpha \in \{ w, g, o \}\f$ ("water", "gas" and "oil") and
 * the components by upper index \f$\kappa \in \{ W, G, O \}\f$
 * ("Water", "Gas" and "Oil"). The model assumes partial miscibility:
 *
 * - Water and the gas phases are immisicible and are assumed to be
 *   only composed of the water and gas components respectively-
 * - The oil phase is assumed to be a mixture of the gas and the oil
 *  components.
 *
 * The densities of the phases are determined by so-called
 * <i>formation volume factors</i>:
 *
 * \f[
 * B_\alpha := \frac{\varrho_\alpha(1\,\text{bar})}{\varrho_\alpha(p_\alpha)}
 * \f]
 *
 * Since the gas and water phases are assumed to be immiscible, this
 * is sufficint to calculate their density. For the formation volume
 * factor of the the oil phase \f$B_o\f$ determines the density of
 * *saturated* oil, i.e. the density of the oil phase if some gas
 * phase is present.
 *
 * The composition of the oil phase is given by the <i>gas formation
 * factor</i> \f$R_s\f$, which defined as the volume of gas at
 * atmospheric pressure that is dissolved in saturated oil at a given
 * pressure:
 *
 * \f[
 * R_s := \frac{x_o^G(p)\,\varrho_{mol,o}(p)}{\varrho_g(1\,\text{bar})}\;.
 * \f]
 *
 * This allows to calculate all quantities required for the
 * mass-conservation equations for each component, i.e.
 *
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha  \right\}
 * - q^\kappa = 0 \;,
 * \f]
 * where \f$\mathrm{v}_\alpha\f$ is the filter velocity of the phase
 * \f$\alpha\f$.
 *
 * By default \f$\mathrm{v}_\alpha\f$ is determined by using the
 * standard multi-phase Darcy approach, i.e.
 * \f[ \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;, \f]
 * although the actual approach which is used can be specified via the
 * \c VelocityModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, VelocityModule, Ewoms::VcfvForchheimerVelocityModule<TypeTag>);
 * \endcode
 *
 * The primary variables used by this model are:
 * - The pressure of the phase with the lowest index
 * - The two saturations of the phases with the lowest indices
 */

#ifndef DUMUX_BLACK_OIL_MODEL_HH
#define DUMUX_BLACK_OIL_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/3p/model.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
// We do not use switchableprimaryvariables!
// #include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>

#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"
#include "localresidual.hh"

namespace Dumux {

/*!
 * \ingroup BlackOilModel
 * \brief Specifies a number properties of two-phase models.
 * \param useCS if we are using the contraint solver
 */
template<bool useCS, bool useMol>
struct BlackOilModelTraits
{
    using Indices = BlackOilIndices;

    static constexpr int numEq() { return 3; }
    static constexpr int numFluidPhases() { return 3; }
    static constexpr int numFluidComponents() { return 3; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; } // Black-oil model does not have diffusion coefficients...
    static constexpr bool enableEnergyBalance() { return false; } // We do not have a Energy equation!

    static constexpr bool useConstraintSolver() { return useCS; }
    static constexpr bool useMoles() { return useMol; }
};

namespace Properties {
// Create new type tags
namespace TTag {
//! The type tags for the pseudo-isothermal three-phase three-component model
struct BlackOil { using InheritsFrom = std::tuple<PorousMediumFlow>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the model traits
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::BlackOil>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static_assert(FluidSystem::numComponents == 3, "Only fluid systems with 3 components are supported by the 3p3c model!");
    static_assert(FluidSystem::numPhases == 3, "Only fluid systems with 3 phases are supported by the 3p3c model!");
public:
     using type = BlackOilModelTraits<getPropValue<TypeTag, Properties::UseConstraintSolver>(), getPropValue<TypeTag, Properties::UseMoles>()>;
};
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::BlackOil> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; };

//! Determines whether a constraint solver should be used explicitly
template<class TypeTag>
struct UseConstraintSolver<TypeTag, TTag::BlackOil> { static constexpr bool value = false; };

//! Set as default that _no_ component mass balance is replaced by the total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::BlackOil> { static constexpr int value = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents(); };
/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
template<class TypeTag>
struct FluidState<TypeTag, TTag::BlackOil>{
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    public:
        using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! The local residual function of the conservation equations
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BlackOil> { using type = BlackOilLocalResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::BlackOil>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;

    using BaseTraits = ThreePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
    template<class BaseTraits, class DT, class EDM>
    struct NCTraits : public BaseTraits
    {
        using DiffusionType = DT;
        using EffectiveDiffusivityModel = EDM;
    };

public:
    using type = BlackOilVolumeVariables<NCTraits<BaseTraits, DT, EDM>>;
};

//! The model after Millington (1961) is used for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::BlackOil> { using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::BlackOil> { using type = BlackOilIOFields; };

//! Use mole fractions in the balance equations by default
template<class TypeTag>
struct UseMoles<TypeTag, TTag::BlackOil> { static constexpr bool value = false; };

} // end namespace Properties
} // end namespace Dumux

#endif
