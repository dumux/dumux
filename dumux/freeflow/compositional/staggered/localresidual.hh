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
  * \ingroup FreeflowNCModel
  * \copydoc Dumux::FreeflowNCResidualImpl
  */
#ifndef DUMUX_FREEFLOW_NC_STAGGERED_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_NC_STAGGERED_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief Element-wise calculation of the multi-component free-flow residual for models using the staggered discretization
 */

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class FreeflowNCResidualImpl;

template<class TypeTag>
class FreeflowNCResidualImpl<TypeTag, DiscretizationMethod::staggered>
: public NavierStokesResidual<TypeTag>
{
    using ParentType = NavierStokesResidual<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr int numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

    using EnergyLocalResidual = typename ParentType::EnergyLocalResidual;

public:
    using ParentType::ParentType;

    //! Evaluate fluxes entering or leaving the cell center control volume.
    template<class VolumeVariables>
    CellCenterPrimaryVariables computeStorageForCellCenter(const Problem& problem,
                                                           const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage(0.0);

        const Scalar density = useMoles ? volVars.molarDensity() : volVars.density();

        // compute storage term of all components
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int eqIdx = compIdx;

            const Scalar massOrMoleFraction = useMoles ? volVars.moleFraction(compIdx) : volVars.massFraction(compIdx);
            const Scalar s =  density * massOrMoleFraction;

            if (eqIdx != ModelTraits::replaceCompEqIdx())
                storage[eqIdx] += s;
        }

        // in case one balance is substituted by the total mass balance
        if(ModelTraits::replaceCompEqIdx() < numComponents)
            storage[ModelTraits::replaceCompEqIdx()] = density;

        EnergyLocalResidual::fluidPhaseStorage(storage, volVars);

        return storage;
    }
};
} // end namespace Dumux

#endif
