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
  * \ingroup NavierStokesNCModel
  * \copydoc Dumux::NavierStokesNCResidualImpl
  */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_NC_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_NC_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag,  DiscretizationMethod discMethod>
class NavierStokesNCResidualImpl;

/*!
 * \ingroup NavierStokesNCModel
 * \brief Element-wise calculation of the multi-component Navier-Stokes residual for models using the staggered discretization
 */
template<class TypeTag>
class NavierStokesNCResidualImpl<TypeTag, DiscretizationMethod::staggered>
: public NavierStokesResidual<TypeTag>
{
    using ParentType = NavierStokesResidual<TypeTag>;
    friend class StaggeredLocalResidual<TypeTag>;
    friend ParentType;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using CellCenterResidual = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);

    enum {
        conti0EqIdx = Indices::conti0EqIdx,
        phaseIdx = Indices::phaseIdx,

        // The index of the component balance equation
        // that gets replaced with the total mass balance
        replaceCompEqIdx = Indices::replaceCompEqIdx
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    using ParentType::ParentType;

    //! Evaluate fluxes entering or leaving the cell center control volume.
    CellCenterPrimaryVariables computeStorageForCellCenter(const Problem& problem,
                                                           const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage(0.0);

        const Scalar density = useMoles ? volVars.molarDensity() : volVars.density();

        // compute storage term of all components
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int eqIdx = conti0EqIdx + compIdx;

            const Scalar massOrMoleFraction = useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx);
            const Scalar s =  density * massOrMoleFraction;

            if (eqIdx != replaceCompEqIdx)
                storage[eqIdx] += s;
        }

        // in case one balance is substituted by the total mass balance
        if(replaceCompEqIdx < numComponents)
            storage[replaceCompEqIdx] = density;

        this->computeStorageForCellCenterNonIsothermal_(std::integral_constant<bool, GET_PROP_VALUE(TypeTag, EnableEnergyBalance) >(),
                                                        problem, scv, volVars, storage);

        return storage;
    }

protected:

    /*!
     * \brief Sets a fixed Dirichlet value for a cell (such as pressure) at the boundary.
     *        This is a provisional alternative to setting the Dirichlet value on the boundary directly.
     */
    void setFixedCell_(CellCenterResidual& residual,
                       const Problem& problem,
                       const SubControlVolume& insideScv,
                       const ElementVolumeVariables& elemVolVars,
                       const BoundaryTypes& bcTypes) const
    {
        ParentType::setFixedCell_(residual, problem, insideScv, elemVolVars, bcTypes);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // get equation index
            const auto eqIdx = conti0EqIdx + compIdx;

            // set a fixed mole fraction for cells
            if(eqIdx != conti0EqIdx && bcTypes.isDirichletCell(eqIdx))
            {
                const auto& insideVolVars = elemVolVars[insideScv];
                const Scalar massOrMoleFraction = useMoles ? insideVolVars.moleFraction(phaseIdx, compIdx) : insideVolVars.massFraction(phaseIdx, compIdx);
                residual[eqIdx] = massOrMoleFraction - problem.dirichletAtPos(insideScv.center())[eqIdx];
            }
        }

    }
};
}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
