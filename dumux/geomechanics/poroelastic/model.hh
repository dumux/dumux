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
 * \ingroup PoroElastic
 * \brief Defines a type tag and some properties for the poroelastic geomechanical model
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

//! The deault model traits of the poro-elastic model
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
