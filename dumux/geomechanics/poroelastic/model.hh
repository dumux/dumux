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
 * \ingroup PoroElastic
 * \brief Defines a type tag and some properties for the poroelastic geomechanical model
 */
#ifndef DUMUX_GEOMECHANICS_POROELASTIC_MODEL_HH
#define DUMUX_GEOMECHANICS_POROELASTIC_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/geomechanics/elastic/indices.hh>
#include <dumux/geomechanics/elastic/model.hh>

#include <dumux/discretization/hookeslaw.hh>
#include <dumux/discretization/effectivestresslaw.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux {

/*!
 * \ingroup Geomechanics
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
NEW_TYPE_TAG(PoroElastic, INHERITS_FROM(Elastic));

//! Use the local residual of the poro-elastic model
SET_TYPE_PROP(PoroElastic, LocalResidual, PoroElasticLocalResidual<TypeTag>);

//! default vtk output fields specific to this model
SET_TYPE_PROP(PoroElastic, VtkOutputFields, PoroElasticVtkOutputFields);

//! The deault model traits of the poro-elastic model
SET_PROP(PoroElastic, ModelTraits)
{
private:
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
    static constexpr int numSC = GET_PROP_TYPE(TypeTag, SolidSystem)::numComponents;
    static constexpr int numFP = GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases;
    static constexpr int numFC = GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents;

public:
    using type = PoroElasticModelTraits<dim, numSC, numFP, numFC>;
};

//! Set the volume variables property
SET_PROP(PoroElastic, VolumeVariables)
{
private:
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using DV = Dune::FieldVector<typename PV::value_type, dim>;
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);

    // we reuse the elastic volume variable traits here
    using Traits = ElasticVolumeVariablesTraits<PV, DV, MT, SST, SSY>;
public:
    using type = PoroElasticVolumeVariables<Traits>;
};

//! Per default, we use effective stresses on the basis of Hooke's Law
SET_PROP(PoroElastic, StressType)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using ElasticStressType = HookesLaw< Scalar, FVGridGeometry >;
public:
    using type = EffectiveStressLaw< ElasticStressType, FVGridGeometry >;
};

} // namespace Properties
} // namespace Dumux

 #endif
