// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElasticLargeDef
 * \brief Large-deformation poroelastic model with volumetric-isochoric split.
 *
 * # Physical model
 *
 * Three-field mixed formulation \f$(\mathbf{d}, p_s, p_f)\f$ for a saturated
 * hyperelastic porous medium undergoing finite deformations.
 *
 * ## Kinematics
 *
 * \f[ \mathbf{F} = \mathbf{I} + \nabla_X \mathbf{d}, \quad
 *     J = \det\mathbf{F}, \quad
 *     \mathbf{C} = \mathbf{F}^\top\mathbf{F}. \f]
 *
 * ## Governing equations
 *
 * **Momentum** (displacement \f$\mathbf{d}\f$):
 * \f[ -\nabla_X \cdot \mathbf{P} = \mathbf{f}, \qquad
 *     \mathbf{P} = \mathbf{P}_\text{iso}(\mathbf{F}) + (p_s - \alpha_B p_f)\, J\, \mathbf{F}^{-\top}, \f]
 * where the isochoric part of the drained skeleton stress \f$\mathbf{P}_\text{iso}\f$ is supplied
 * by the problem via `isochoricFirstPiolaKirchhoffStressTensor(F, elemDisc, ipGlobal)` and
 * the body force \f$\mathbf{f}\f$ via the problem's `source`.
 *
 * **Solid pressure** (\f$p_s\f$, algebraic):
 * \f[ p_s^{\mathrm{eq}}(J) - p_s = 0 \f]
 * (volumetric constitutive law; \f$p_s^{\mathrm{eq}}\f$ supplied by the problem via
 * `volumetricPressure(J)`).
 *
 * **Fluid mass balance** (\f$p_f\f$, transient):
 * \f[ \frac{\partial}{\partial t}\!\left(\alpha_B J + S_p p_f\right)
 *     - \nabla_X \cdot \!\left(J\,\frac{\kappa}{\mu_f}\,\mathbf{C}^{-1} \nabla_X p_f\right) = q, \f]
 * with the Biot coefficient \f$\alpha_B\f$, the storage coefficient \f$S_p = 1/M\f$
 * (inverse Biot modulus, \f$S_p = 0\f$ for incompressible constituents) and the fluid
 * source \f$q\f$ given by the problem's `source`. In the small-strain limit
 * \f$J \approx 1 + \nabla \cdot \mathbf{d}\f$ this is Biot's storage equation.
 *
 * The spatial parameters provide `biotCoefficient(element, elemDisc, ipGlobal)`,
 * `storageCoefficient(element, elemDisc, ipGlobal)`,
 * `permeability(element, elemDisc, ipGlobal, J)` and
 * `fluidViscosity(element, elemDisc)`.
 *
 * The type tags `PoroElasticLargeDefMomentumModel`, `PoroElasticLargeDefSolidPressureModel` and
 * `PoroElasticLargeDefFluidPressureModel` are coupled via `PoroElasticLargeDefCouplingManager`
 * and assembled with `Experimental::MultiStageMultiDomainAssembler`.
 */
#ifndef DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_MODEL_HH
#define DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_MODEL_HH

#include <tuple>

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "localresidual.hh"
#include "variables.hh"

namespace Dumux {

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Index set for the momentum subdomain.
 */
template<int dim>
struct PoroElasticLargeDefMomentumIndices
{
    static constexpr int displacement(int i) { return i; }
    static constexpr int momentum(int i) { return i; }
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Model traits for the momentum subdomain (\f$\dim\f$ displacement DOFs per node).
 */
template<int dim>
struct PoroElasticLargeDefMomentumModelTraits
{
    using Indices = PoroElasticLargeDefMomentumIndices<dim>;
    static constexpr int numEq() { return dim; }
};

template<class PV, class MT>
struct PoroElasticLargeDefMomentumVVTraits
{ using PrimaryVariables = PV; using ModelTraits = MT; };

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Index set for the solid pressure subdomain.
 */
struct PoroElasticLargeDefSolidPressureIndices
{
    static constexpr int pressureIdx = 0;
    static constexpr int pressureEqIdx = 0;
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Model traits for the solid pressure subdomain (1 pressure DOF per node).
 */
struct PoroElasticLargeDefSolidPressureModelTraits
{
    using Indices = PoroElasticLargeDefSolidPressureIndices;
    static constexpr int numEq() { return 1; }
};

template<class PV, class MT>
struct PoroElasticLargeDefSolidPressureVVTraits
{ using PrimaryVariables = PV; using ModelTraits = MT; };

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Index set for the fluid pressure subdomain.
 */
struct PoroElasticLargeDefFluidPressureIndices
{
    static constexpr int pressureIdx = 0;
    static constexpr int pressureEqIdx = 0;
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Model traits for the fluid pressure subdomain (1 fluid pressure DOF per node).
 */
struct PoroElasticLargeDefFluidPressureModelTraits
{
    using Indices = PoroElasticLargeDefFluidPressureIndices;
    static constexpr int numEq() { return 1; }
};

template<class PV, class MT>
struct PoroElasticLargeDefFluidPressureVVTraits
{ using PrimaryVariables = PV; using ModelTraits = MT; };

} // end namespace Dumux

namespace Dumux::Properties::TTag {

struct PoroElasticLargeDefMomentumModel
{ using InheritsFrom = std::tuple<ModelProperties>; };

struct PoroElasticLargeDefSolidPressureModel
{ using InheritsFrom = std::tuple<ModelProperties>; };

struct PoroElasticLargeDefFluidPressureModel
{ using InheritsFrom = std::tuple<ModelProperties>; };

} // end namespace Dumux::Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PoroElasticLargeDefMomentumModel>
{
    using type = PoroElasticLargeDefMomentumModelTraits<
        GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PoroElasticLargeDefMomentumModel>
{ using type = PoroElasticLargeDefMomentumLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PoroElasticLargeDefMomentumModel>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using type = PoroElasticLargeDefMomentumVariables<
        PoroElasticLargeDefMomentumVVTraits<PV, MT>>;
};

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PoroElasticLargeDefSolidPressureModel>
{ using type = PoroElasticLargeDefSolidPressureModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PoroElasticLargeDefSolidPressureModel>
{ using type = PoroElasticLargeDefSolidPressureLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PoroElasticLargeDefSolidPressureModel>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using type = PoroElasticLargeDefSolidPressureVariables<
        PoroElasticLargeDefSolidPressureVVTraits<PV, MT>>;
};

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PoroElasticLargeDefFluidPressureModel>
{ using type = PoroElasticLargeDefFluidPressureModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PoroElasticLargeDefFluidPressureModel>
{ using type = PoroElasticLargeDefFluidPressureLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PoroElasticLargeDefFluidPressureModel>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using type = PoroElasticLargeDefFluidPressureVariables<
        PoroElasticLargeDefFluidPressureVVTraits<PV, MT>>;
};

} // end namespace Dumux::Properties

#endif
