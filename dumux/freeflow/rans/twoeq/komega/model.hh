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
 * \ingroup KOmegaModel
 *
 * \brief A single-phase, isothermal k-omega 2-Eq. model
 *
 * \copydoc RANSModel
 *
 * Two additional PDEs, one for the turbulentKineticEnergy (k) and a second for the dissipation (omega)
 * are used to calculate the eddy viscosity for this model.
 * The model is taken from Wilcox, 2008 \cite Wilcox2008a.
 *
 * Turbulent Kinetic Energy balance:
 * \f[
 * \frac{\partial \left( k \right)}{\partial t}
 * + \nabla \cdot \left( \mathbf{v} k \right)
 * - \nabla \cdot \left[ \left( \nu +  \sigma_\textrm{k} \nu_\textrm{t} \right) \nabla k \right]
 * - P
 * + \beta_k^{*} k \omega
 * = 0
 * \f]
 * with \f$ P = 2 \nu_\textrm{t} \mathbf{S} \cdot \mathbf{S} \f$
 * and \f$ S_{ij} = \frac{1}{2} \left[ \frac{\partial}{\partial x_i} v_j + \frac{\partial}{\partial x_j} v_i \right] \f$
 * based on \f$ a_{ij} \cdot b_{ij} = \sum_{i,j} a_{ij} b_{ij} \f$.
 *
 * Dissipation balance:
 * \f[
 * \frac{\partial \left( \omega \right)}{\partial t}
 * + \nabla \cdot \left( \mathbf{v} \omega \right)
 * - \nabla \cdot \left[ \left( \nu + \sigma_{\omega} \nu_\textrm{t} \right) \nabla \omega \right]
 * - \alpha \frac{\omega}{k} P
 * + \beta_{\omega} \omega^2
 * - \frac{\sigma_d}{\omega} \nabla k \nabla \omega
 * = 0
 * \f]
 *
 * The kinematic eddy viscosity \f$ \nu_\textrm{t} \f$ is calculated as follows:
 * \f[ \nu_\textrm{t} = \frac{k}{\tilde{\omega}} \f]
 *
 * With a limited dissipation:
 * \f[ \tilde{\omega} = \textrm{max} \left\{ \omega, 0.875 \sqrt{\frac{P}{\nu_\textrm{t} \beta_\textrm{k}}} \right\} \f]
 *
 * And a cross-diffusion coefficient \f$ \sigma_\textrm{d} \f$
 * \f[
 *   \sigma_\text{d} =
 *   \begin{cases}
 *     0     & \mbox{, if } \; \nabla k \cdot \nabla \omega \le 0 \\
 *     0.125 & \mbox{, if } \; \nabla k \cdot \nabla \omega >   0
 *   \end{cases}.
 * \f]
 */

#ifndef DUMUX_KOMEGA_MODEL_HH
#define DUMUX_KOMEGA_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>

#include "fluxvariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{
namespace Properties {

/*!
 *\ingroup KOmegaModel
 * \brief Traits for the k-omega model
 */
template<int dimension>
struct KOmegaModelTraits
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions,
    //! one mass balance equation and two turbulent transport equations
    static constexpr int numEq() { return dim()+1+2; }

    //! The number of phases is always 1
    static constexpr int numPhases() { return 1; }

    //! The number of components
    static constexpr int numComponents() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! The one-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return true; }

    //! The model is isothermal
    static constexpr bool enableEnergyBalance() { return false; }

    //! The indices
    using Indices = KOmegaIndices<dim(), numComponents()>;

    //! The model includes a limiter to the production term
    static constexpr bool enableKOmegaProductionLimiter() { return true; }
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal k-omega single phase model
///////////////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal k-omega model
NEW_TYPE_TAG(KOmega, INHERITS_FROM(RANS));

//!< states some specifics of the isothermal k-omega model
SET_PROP(KOmega, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = KOmegaModelTraits<dim>;
};

//! The flux variables
SET_TYPE_PROP(KOmega, FluxVariables, KOmegaFluxVariables<TypeTag>);

//! The local residual
SET_TYPE_PROP(KOmega, LocalResidual, KOmegaResidual<TypeTag>);

//! Set the volume variables property
SET_PROP(KOmega, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = KOmegaVolumeVariables<Traits, NSVolVars>;
};

//! The specific vtk output fields
SET_PROP(KOmega, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
public:
    using type = KOmegaVtkOutputFields<FVGridGeometry>;
};

///////////////////////////////////////////////////////////////////////////
// default property values for the non-isothermal k-omega single phase model
///////////////////////////////////////////////////////////////////////////


//! The type tag for the single-phase, non-isothermal k-omega 2-Eq. model
NEW_TYPE_TAG(KOmegaNI, INHERITS_FROM(RANSNI));

//! Set the volume variables property
SET_PROP(KOmegaNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = KOmegaVolumeVariables<Traits, NSVolVars>;
};

// \}
}

} // end namespace

#endif // DUMUX_KOMEGA_MODEL_HH
