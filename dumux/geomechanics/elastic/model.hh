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
 * \ingroup Elastic
 * \brief Defines a type tag and some properties for the elastic geomechanical model
 */
#ifndef DUMUX_GEOMECHANICS_ELASTIC_MODEL_HH
#define DUMUX_GEOMECHANICS_ELASTIC_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/geomechanics/properties.hh>
#include <dumux/flux/hookeslaw.hh>

#include "indices.hh"
#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

/*!
 * \ingroup Elastic
 * \brief Specifies a number properties of the elastic model
 */
template< int dim, int numSolidComp >
struct ElasticModelTraits
{
    //! export the type encapsulating indices
    using Indices = ElasticIndices;
    //! the number of equations is equal to grid dimension
    static constexpr int numEq() { return dim; }
    //! This model does not consider fluid phases
    static constexpr int numFluidPhases() { return 0; }
    //! This model does not consider fluid phases
    static constexpr int numFluidComponents() { return 0; }
    //! We have one solid phase here
    static constexpr int numSolidComponents() { return numSolidComp; }

    //! Energy balance not yet implemented
    static constexpr bool enableEnergyBalance() { return false; }
};

/*!
 * \ingroup Elastic
 * \brief Traits class for the volume variables of the elastic model.
 *
 * \tparam PV The type used for primary variables
 * \tparam DV The type used for displacement vectors
 * \tparam MT The model traits
 * \tparam SST The solid state
 * \tparam SSY The solid system
 */
template<class PV, class DV, class MT, class SST, class SSY>
struct ElasticVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using DisplacementVector = DV;
    using ModelTraits = MT;
    using SolidState = SST;
    using SolidSystem = SSY;
};

namespace Properties {

//! Type tag for the elastic geomechanical model
// Create new type tags
namespace TTag {
struct Elastic { using InheritsFrom = std::tuple<Geomechanics>; };
} // end namespace TTag

//! Use the local residual of the elastic model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Elastic> { using type = ElasticLocalResidual<TypeTag>; };

//! The model traits of the elastic model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Elastic>
{
    using type = ElasticModelTraits< GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension,
                                     GetPropType<TypeTag, Properties::SolidSystem>::numComponents >;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Elastic>
{
private:
    static constexpr int dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension;
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using DV = Dune::FieldVector<typename PV::value_type, dim>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using Traits = ElasticVolumeVariablesTraits<PV, DV, MT, SST, SSY>;
public:
    using type = ElasticVolumeVariables<Traits>;
};

//! By default, we use hooke's law for stress evaluations
template<class TypeTag>
struct StressType<TypeTag, TTag::Elastic>
{
    using type = HookesLaw< GetPropType<TypeTag, Properties::Scalar>,
                            GetPropType<TypeTag, Properties::GridGeometry> >;
};

} // namespace Properties
} // namespace Dumux

 #endif
