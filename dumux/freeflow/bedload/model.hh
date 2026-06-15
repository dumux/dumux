// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransportModel
 *
 * \brief A bedload transport model.
 *
 * The sediment balance equation, also known as Exner equation, balances the temporal change of the bedsurface \f$ z_b \f$
 * with the spatial change of the bedload transport rate \f$ \mathbf{q_b} \f$:
 *
 * \f[
 * \frac{\partial{z_b}}{\partial{t}} + \zeta \frac{\partial{q_{bx}}}{\partial{x}} + \zeta \frac{\partial{q_{by}}}{\partial{y}} = 0
 * \f]
 *
 * where \f$\zeta = \frac{1}{1-p}\f$ and \f$p\f$ is the porosity of the bed material.
 *
 * However, this model uses the mass based version of the sediment balance equation, which is more convenient,
 * when dealing with several grain classes. Using a first-order cell-centred finite volume discretisation in space
 * and a first order implicit Euler discretisation in time the mass based version of the sediment balance equation
 * reads as:
 *
 * \f[
 * m_i^\text{n+1} = m_i^\text{n} - \rho_{s,i} \Delta t \sum_{k=1}^{3} \mathbf{q}^\text{n+1}_{b,ik} \cdot \mathbf{n}_k l_k
 * \f]
 *
 * When using several grain classes the sediment balance equation must be solved for each grain class \f$ i \f$.
 * \f$ \text{n+1} \f$ and \f$ \text{n} \f$ denote the new and the old time levels, \f$ \rho_s \f$ the grain denisty.
 * The unit outward normal vector is designated as \f$ /mathbf{n} \f$, and the length of a face is referred to as \f$ l \f$.
 * The index \f$ k \f$ denotes the corresponding face. Consequently, \f$ q_{b,k} \f$ denote the bedload flux at the face \f$ k \f$.
 * The calculation of this flux is the primary objective within the finite volume approach and determines the stability
 * and accuracy of the scheme.
 *
 * The bedsurface is calculated using the mass \f$ m \f$ of all grain classes:
 *
 * \f[
 * z_b = z_{fix} + \zeta \sum_i \frac{m_i}{\rho_{\text{s},i}}
 * \f]
 *
 * where \f$z_{fix}\f$ is the bottom of the erodible layer.
 *
 */
#ifndef DUMUX_BEDLOAD_MODEL_HH
#define DUMUX_BEDLOAD_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/flux/fluxvariablescaching.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "indices.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup BedloadModel
 * \brief Traits for the bedload transport model.
 */
 template<int nGrainClasses>
struct BedloadModelTraits
{
    using Indices = BedloadIndices;

    static constexpr int numGrainClasses() { return nGrainClasses; }
    static constexpr int numEq() { return nGrainClasses; }
    static constexpr int numPhases() { return 1; }

    // Enable sediment
    static constexpr bool enableBedload() { return true; }
};

/*!
 * \ingroup BedloadModel
 * \brief Traits class for the volume variables of the bedload transport model.
 *
 * \tparam PV The type used for primary variables
 * \tparam MT The model traits
 */
template<class PV,
         class MT>
struct BedloadVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};


namespace Properties {

//! Type tag for the bedload transport model
// Create new type tags
namespace TTag {
struct Bedload{ using InheritsFrom = std::tuple<ModelProperties>; };
}// end namespace TTag

//////////////////////////////////////////////////////////////////
// Define properties
//////////////////////////////////////////////////////////////////

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Bedload>
{ using type = BedloadModelTraits<TypeTag::nGrainClasses>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Bedload>
{ using type = BedloadLocalResidual<TypeTag>; };

template<class TypeTag>
struct FluxVariables<TypeTag, TTag::Bedload>
{ using type = BedloadFluxVariables<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Bedload>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = BedloadVolumeVariablesTraits<PV, MT>;
public:
    using type = BedloadVolumeVariables<Traits>;
};

template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::Bedload>
{ using type = FluxVariablesCaching::EmptyCache< GetPropType<TypeTag, Properties::Scalar> >; };

template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::Bedload>
{ using type = FluxVariablesCaching::EmptyCacheFiller; };

template<class TypeTag>
struct IOFields<TypeTag, TTag::Bedload>
{ using type = BedloadIOFields; };

} // namespace Properties
} // namespace Dumux

#endif
