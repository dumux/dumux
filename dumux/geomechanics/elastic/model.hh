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
 * \ingroup Properties
 * \ingroup Geomechanics
 * \ingroup Elastic
 * \brief Defines a type tag and some properties for the elastic geomechanical model
 */
#ifndef DUMUX_GEOMECHANICS_ELASTIC_MODEL_HH
#define DUMUX_GEOMECHANICS_ELASTIC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/geomechanics/properties.hh>
#include <dumux/discretization/hookeslaw.hh>

#include "indices.hh"
#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

/*!
 * \ingroup Geomechanics
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
 * \tparam MT The model traits
 * \tparam LP The class used for storing the lame parameters
 * \tparam SST The solid state
 * \tparam SSY The solid system
 */
template<class PV, class MT, class LP, class SST, class SSY>
struct ElasticVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
    using LameParams = LP;
    using SolidState = SST;
    using SolidSystem = SSY;
};

namespace Properties {

//! Type tag for the elastic geomechanical model
NEW_TYPE_TAG(Elastic, INHERITS_FROM(Geomechanics));

//! Use the local residual of the elastic model
SET_TYPE_PROP(Elastic, LocalResidual, ElasticLocalResidual<TypeTag>);

//! By default, we use hooke's law for stress evaluations
SET_TYPE_PROP(Elastic, ModelTraits, ElasticModelTraits< GET_PROP_TYPE(TypeTag, GridView)::dimension,
                                                        GET_PROP_TYPE(TypeTag, SolidSystem)::numComponents >);

//! Set the volume variables property
SET_PROP(Elastic, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using LP = typename GET_PROP_TYPE(TypeTag, SpatialParams)::LameParams;
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using Traits = ElasticVolumeVariablesTraits<PV, MT, LP, SST, SSY>;
public:
    using type = ElasticVolumeVariables<Traits>;
};

} // namespace Properties
} // namespace Dumux

 #endif
