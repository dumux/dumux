// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElastic
 * \brief A poroelastic geomechanical model
 *
 * This model describes the deformation of a porous medium using the theory of linear poroelasticity.
 * The momentum balance equation of a porous medium can be expressed by
 \f[
 \nabla\cdot\boldsymbol{\sigma_{\mathrm{eff}}} + \rho \mathbf{g} + \mathbf{f} = \rho\ddot{\mathbf{u}},
 \f]
 * where \f$ \boldsymbol{\sigma_{\mathrm{eff}}} \f$ is the effective stress tensor,
 * \f$ \rho = (1 - \phi) \rho_s + \phi \rho_f \f$ is the average density of solids and fluids within the porous medium,
 * \f$ \mathbf{f} \f$ in \f$ \mathrm{N/m^3} \f$  is the external force acting on the body per unit volume (e.g. magnetism),
 * and \f$ \mathbf{u} = \mathbf{x} - \mathbf{x}_{\mathrm{initial}} \f$ is the displacement,
 * defined as the difference in material points \f$ \mathbf{x} \f$ and \f$ \mathbf{x}_{\mathrm{initial}} \f$
 * in the deformed and undeformed (initial) state, respectively. The model assumes quasi-static conditions,
 * that is, the above momentum balance equation is solved under the assumption that the acceleration term
 * \f$ \rho\ddot{\mathbf{u}} \approx 0\f$.
 *
 * Using the concept of the effective stress, the effective stress tensor \f$ \boldsymbol{\sigma_{\mathrm{eff}}} \f$ is
 * determined by the stress tensor \f$ \boldsymbol{\sigma} \f$ , the effective pore pressure \f$ p_{\mathrm{eff}} \f$ and the Biot's coefficient \f$ \alpha \f$ :
 \f[
 \boldsymbol{\sigma_{\mathrm{eff}}} = \boldsymbol{\sigma} - \alpha p_{\mathrm{eff}} \mathbf{I}
 \f]
 *
 * Per default, Hookes' Law is used for expressing the stress tensor \f$ \boldsymbol{\sigma} \f$ as a function of the
 * displacement:
 \f[
 \boldsymbol{\sigma} = \lambda\mathrm{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2G \boldsymbol{\varepsilon},
 \f]
 * with
 \f[
 \boldsymbol{\varepsilon} = \frac{1}{2} \left[ \nabla\mathbf{u} + (\nabla\mathbf{u})^{\mathrm{T}} \right].
 \f]
 *
 * Primary variables are the displacements in each direction \f$ \mathbf{u} \f$.
 * Gravity can be enabled or disabled via a runtime parameter.
 *
 */
#ifndef DUMUX_GEOMECHANICS_POROELASTIC_MODEL_HH
#define DUMUX_GEOMECHANICS_POROELASTIC_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/geomechanics/elastic/indices.hh>
#include <dumux/geomechanics/elastic/model.hh>

#include <dumux/flux/hookeslaw.hh>
#include <dumux/flux/effectivestresslaw.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup PoroElastic
 * \brief Specifies a number properties of the poroelastic model
 */
template< int dim, int numSC, int numFP, int numFC >
struct PoroElasticModelTraits
{
    //! export the type encapsulating indices
    using Indices = ElasticIndices;
    //! the number of equations is equal to grid dimension
    static constexpr int numEq() { return dim; }
    //! This model does not consider fluid phases
    static constexpr int numFluidPhases() { return numFP; }
    //! This model does not consider fluid phases
    static constexpr int numFluidComponents() { return numFC; }
    //! We have one solid phase here
    static constexpr int numSolidComponents() { return numSC; }

    //! Energy balance not yet implemented
    static constexpr bool enableEnergyBalance() { return false; }
};

namespace Properties {

//! Type tag for the poro-elastic geomechanical model
// Create new type tags
namespace TTag {
struct PoroElastic { using InheritsFrom = std::tuple<Elastic>; };
} // end namespace TTag

//! Use the local residual of the poro-elastic model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PoroElastic> { using type = PoroElasticLocalResidual<TypeTag>; };

//! default vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PoroElastic> { using type = PoroElasticIOFields; };

//! The default model traits of the poro-elastic model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PoroElastic>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int numSC = GetPropType<TypeTag, Properties::SolidSystem>::numComponents;
    static constexpr int numFP = GetPropType<TypeTag, Properties::FluidSystem>::numPhases;
    static constexpr int numFC = GetPropType<TypeTag, Properties::FluidSystem>::numComponents;

public:
    using type = PoroElasticModelTraits<dim, numSC, numFP, numFC>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PoroElastic>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using DV = Dune::FieldVector<typename PV::value_type, dim>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;

    // we reuse the elastic volume variable traits here
    using Traits = ElasticVolumeVariablesTraits<PV, DV, MT, SST, SSY>;
public:
    using type = PoroElasticVolumeVariables<Traits>;
};

//! Per default, we use effective stresses on the basis of Hooke's Law
template<class TypeTag>
struct StressType<TypeTag, TTag::PoroElastic>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElasticStressType = HookesLaw< Scalar, GridGeometry >;
public:
    using type = EffectiveStressLaw< ElasticStressType, GridGeometry >;
};

} // namespace Properties
} // namespace Dumux

 #endif
