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
 * \ingroup RANSModel
 *
 * \brief A single-phase, isothermal Reynolds-Averaged Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Reynolds-Averaged
 * Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}}) = \nabla \cdot (\mu_\textrm{eff} (\nabla \textbf{v} + \nabla \textbf{v}^{\text{T}}))
   - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * The effective viscosity is composed of the fluid and the eddy viscosity:
 * \f[
 *    \mu_\textrm{eff} = \mu + \mu_\textrm{t}
 * \f].
 */
#ifndef DUMUX_RANS_MODEL_HH
#define DUMUX_RANS_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

#include "iofields.hh"

namespace Dumux {

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Reynolds-Averaged Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes model
struct RANS { using InheritsFrom = std::tuple<NavierStokes>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
/*!
 * \ingroup RANSModel
 * \brief Traits for the Reynolds-averaged Navier-Stokes model
 *
 * \tparam dimension The dimension of the problem
 */
template<int dimension>
struct RANSModelTraits : NavierStokesModelTraits<dimension>
{
    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }
};

//! The model traits of the isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::RANS>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = RANSModelTraits<dim>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::RANS> { using type = RANSIOFields; };

//////////////////////////////////////////////////////////////////
// Property values for non-isothermal Reynolds-averaged Navier-Stokes model
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes model
struct RANSNI { using InheritsFrom = std::tuple<RANS>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::RANSNI>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;

    using IsothermalTraits = RANSModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! The specific non-isothermal I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::RANSNI> { using type = FreeflowNonIsothermalIOFields<RANSIOFields, true/*turbulenceModel*/>; };

//! Use Fourier's Law as default heat conduction type
template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::RANSNI> { using type = FouriersLaw<TypeTag>; };

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
