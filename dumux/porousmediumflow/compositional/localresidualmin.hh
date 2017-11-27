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
 *
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model that also considers solid phases.
 */
#ifndef DUMUX_COMPOSITIONAL_LOCAL_RESIDUAL_MINERAL_HH
#define DUMUX_COMPOSITIONAL_LOCAL_RESIDUAL_MINERAL_HH

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(ReplaceCompEqIdx);
} // end namespace Properties

/*!
 * \ingroup Implicit
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 *
 */
template<class TypeTag>
class CompositionalMinLocalResidual: public GET_PROP_TYPE(TypeTag, CompositionalLocalResidual)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, CompositionalLocalResidual);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numSPhases = GET_PROP_VALUE(TypeTag, NumSPhases);
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    enum { conti0EqIdx = Indices::conti0EqIdx };

    //! The index of the component balance equation that gets replaced with the total mass balance
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < numComponents;

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    ResidualVector computeStorage(const Problem& problem,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars) const
    {
        ResidualVector storage(0.0);

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // compute storage term of all components within all fluid phases

        ParentType::computeStorage(problem, scv, volVars);
        for (int phaseIdx = numPhases; phaseIdx < numPhases + numSPhases; ++phaseIdx)
        {
            auto eqIdx = conti0EqIdx + numComponents-numPhases + phaseIdx;
            storage[eqIdx] += volVars.precipitateVolumeFraction(phaseIdx)
                             * massOrMoleDensity(volVars, phaseIdx);
        }

        return storage;
    }

protected:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }
};

} // end namespace Dumux

#endif
