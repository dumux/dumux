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
  * \copydoc Dumux::NavierStokesNCFluxVariablesImpl
  */
#ifndef DUMUX_NAVIERSTOKES_NC_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_NC_STAGGERED_FLUXVARIABLES_HH

#include <numeric>
#include <dumux/common/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/navierstokes/fluxvariables.hh>

namespace Dumux
{

// forward declaration
template<class TypeTag, DiscretizationMethod Method>
class NavierStokesNCFluxVariablesImpl;

/*!
 * \ingroup NavierStokesNCModel
 * \brief The flux variables class for the multi-component Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag>
class NavierStokesNCFluxVariablesImpl<TypeTag, DiscretizationMethod::staggered>
: public NavierStokesFluxVariables<TypeTag>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);

    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    static constexpr auto numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    // The index of the component balance equation that gets replaced with the total mass balance
    static const int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    static const int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

    using ParentType = NavierStokesFluxVariables<TypeTag>;

    enum { conti0EqIdx = Indices::conti0EqIdx };

public:

    /*!
    * \brief Computes the flux for the cell center residual.
    */
    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const ElementFaceVariables& elemFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux(0.0);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // get equation index
            const auto eqIdx = conti0EqIdx + compIdx;

            bool isOutflow = false;
            if(scvf.boundary())
            {
                const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
                    if(bcTypes.isOutflow(eqIdx))
                        isOutflow = true;
            }

            auto upwindTerm = [compIdx](const auto& volVars)
            {
                const auto density = useMoles ? volVars.molarDensity() : volVars.density();
                const auto fraction =  useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx);
                return density * fraction;
            };

            flux[eqIdx] = ParentType::advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTerm, isOutflow);
        }

        // in case one balance is substituted by the total mass balance
        if (replaceCompEqIdx < numComponents)
        {
            flux[replaceCompEqIdx] = std::accumulate(flux.begin(), flux.end(), 0.0);
        }

        flux += MolecularDiffusionType::diffusiveFluxForCellCenter(problem, fvGeometry, elemVolVars, scvf);
        return flux;
    }
};

} // end namespace

#endif // DUMUX_NAVIERSTOKES_NC_STAGGERED_FLUXVARIABLES_HH
