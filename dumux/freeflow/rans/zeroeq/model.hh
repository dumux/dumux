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
 * \ingroup ZeroEqModel
 *
 * \brief A single-phase, isothermal Reynolds-Averaged Navier-Stokes 0-Eq. model
 *
 * \copydoc RANSModel
 *
 * These models calculate the eddy viscosity without solving additional PDEs,
 * only based on the wall distance and the velocity gradient.
 *
 * The following models are available:
 *  -# Prandtl's mixing length, e.g. \cite Oertel2012a
 *  -# Van-Driest modification, \cite vanDriest1956a and \cite Hanna1981a
 *  -# Baldwin-Lomax, \cite Baldwin1978a
 */
#ifndef DUMUX_ZEROEQ_MODEL_HH
#define DUMUX_ZEROEQ_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>
#include <dumux/freeflow/navierstokes/volumevariables.hh>

#include "problem.hh"
#include "volumevariables.hh"

namespace Dumux::Properties {

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal RANS 0-Eq. model
///////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes 0-Eq. model
struct ZeroEq { using InheritsFrom = std::tuple<RANS>; };
} // end namespace TTag

/*!
 * \ingroup ZeroEqModel
 * \brief Traits for the ZeroEq model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct ZeroEqModelTraits : RANSModelTraits<dimension>
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::zeroeq; }
};

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ZeroEq>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = ZeroEqModelTraits<dim>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ZeroEq>
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
    using type = ZeroEqVolumeVariables<Traits, NSVolVars>;
};

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal RANS 0-Eq. model
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, non-isothermal Reynolds-Averaged Navier-Stokes model
struct ZeroEqNI { using InheritsFrom = std::tuple<RANSNI>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ZeroEqNI>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = ZeroEqModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ZeroEqNI>
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
    using type = ZeroEqVolumeVariables<Traits, NSVolVars>;
};

} // end namespace Dumux::Properties

#endif
