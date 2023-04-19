// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SSTModel
 *
 * \brief A single-phase, isothermal SST (Shear Stress Transport) -Eq. model
 *
 * \copydoc RANSModel
 *
 * Two additional PDEs, one for the turbulentKineticEnergy (k) and a second for the dissipation (omega)
 * are used to calculate the eddy viscosity for this model.
 * The model is taken from Menter, 1994 \cite Menter1994a.
 *
 * Turbulent Kinetic Energy balance:
 * \f[
 * \frac{\partial \varrho k}{\partial t}
 * + \nabla \cdot \left( \varrho \boldsymbol{u} k \right)
 * - 2\mu_t \boldsymbol{S}\cdot\boldsymbol{S}
 * + \beta^* \varrho\omega k
 * - \nabla\cdot \left[ \left( \mu + \sigma_k \mu_t \right)\nabla k\right]
 * = 0
 * \f]
 * and \f$ S_{ij} = \frac{1}{2} \left[ \frac{\partial}{\partial x_i} u_j + \frac{\partial}{\partial x_j} u_i \right] \f$
 * based on \f$ a_{ij} \cdot b_{ij} = \sum_{i,j} a_{ij} b_{ij} \f$.
 *
 * Dissipation(rate) balance:
 * \f[
 * \frac{\partial \varrho \omega}{\partial t}
 * + \nabla \cdot \left( \varrho \boldsymbol{u} \omega \right)
 * - \frac{\gamma}{\nu_t}\left(2\mu_t \boldsymbol{S}\cdot\boldsymbol{S}\right)
 * + \beta^* \varrho\omega^2
 * - \nabla\cdot \left[ \left( \mu + \sigma_k \mu_t \right)\nabla \omega\right]
 * - 2\varrho\left( 1-F_1\right) \sigma_{\omega 2} \frac{1}{\omega}\nabla k \nabla \omega
 * = 0
 * \f]
 *
 * The dynamic eddy viscosity \f$ \mu_\textrm{t} \f$ is calculated as follows:
 * \f[ \mu_t = \varrho \frac{a_1 k}{max\left( a_1 \omega; \Omega F_2\right)} \f]
 * and \f$ a_1 = 0.31 \f$
 * and \f$ \Omega = \sqrt{2\boldsymbol{\Omega}\cdot\boldsymbol{\Omega}} \text{ with } \boldsymbol{\Omega} = \frac{1}{2} \left( \nabla\boldsymbol{u} - \nabla^T\boldsymbol{u}\right) \f$
 * and \f$ F_2 = tanh\left( arg_2^2\right) \f$
 * and \f$ arg_2 = max\left( 2\frac{\sqrt{k}}{0.09\omega y}; \frac{500\nu}{y^2\omega} \right) \f$
 * where y is the distance to the closest wall and \f$ \nu \f$ is the kinematic viscosity.
 */

#ifndef DUMUX_SST_MODEL_HH
#define DUMUX_SST_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>
#include <dumux/freeflow/rans/twoeq/indices.hh>

#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/rans/problem.hh>

#include "fluxvariables.hh"
#include "localresidual.hh"
#include "volumevariables.hh"
#include "iofields.hh"

namespace Dumux {
namespace Properties {

/*!
 *\ingroup SSTModel
 * \brief Traits for the sst model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct SSTModelTraits : RANSModelTraits<dimension>
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions,
    //! one mass balance equation and two turbulent transport equations
    static constexpr int numEq() { return dim()+1+2; }

    //! The number of components
    static constexpr int numFluidComponents() { return 1; }

    //! The indices
    using Indices = RANSTwoEqIndices<dim(), numFluidComponents()>;

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::sst; }
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal SST single phase model
///////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal SST model
struct SST { using InheritsFrom = std::tuple<RANS>; };
} // end namespace TTag

//! states some specifics of the isothermal SST model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::SST>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = SSTModelTraits<dim>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::SST>
{
private:
    using BaseFluxVariables = NavierStokesFluxVariables<TypeTag>;
public:
    using type = SSTFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::SST>
{
private:
    using BaseLocalResidual = NavierStokesResidual<TypeTag>;
public:
    using type = SSTResidual<TypeTag, BaseLocalResidual>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::SST>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = SSTVolumeVariables<Traits, NSVolVars>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::SST> { using type = SSTIOFields; };

///////////////////////////////////////////////////////////////////////////
// default property values for the non-isothermal SST single phase model
///////////////////////////////////////////////////////////////////////////


// Create new type tags
namespace TTag {
//! The type tag for the single-phase, non-isothermal SST 2-Eq. model
struct SSTNI { using InheritsFrom = std::tuple<SST, RANSNI>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::SSTNI>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = SSTModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::SSTNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = SSTVolumeVariables<Traits, NSVolVars>;
};

//! The specific non-isothermal I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::SSTNI> { using type = FreeflowNonIsothermalIOFields<SSTIOFields, true/*turbulenceModel*/>; };

} // end properties
} // end namespace

#endif
