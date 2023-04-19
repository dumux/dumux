// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MineralizationModel
 * \brief Element-wise calculation of the local residual for problems using a
 *        compositional model that also considers mineralization of solid phases.
 */

#ifndef DUMUX_COMPOSITIONAL_MINERALIZATION_LOCAL_RESIDUAL_HH
#define DUMUX_COMPOSITIONAL_MINERALIZATION_LOCAL_RESIDUAL_HH

#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>

namespace Dumux
{
/*!
 * \ingroup MineralizationModel
 * \brief Element-wise calculation of the local residual for problems
 *        using a one/two-phase n-component mineralization model.
 */
template<class TypeTag>
class MineralizationLocalResidual: public CompositionalLocalResidual<TypeTag>
{
    using ParentType = CompositionalLocalResidual<TypeTag>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int numSolidComps =  ModelTraits::numSolidComps();
    static constexpr int numInertSolidComps =  ModelTraits::numInertSolidComps();
    static constexpr int numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * We consider the volume-average here (e.g. phase mass inside a
     * sub control volume divided by the volume). The volume is multiplied
     * onto it afterwards in the local residual of the respective spatial
     * discretization scheme.
     *
     * \param problem The problem (Initial/Boundary conditions...) to be solved
     * \param scv The sub-control volume of the finite volume grid
     * \param volVars The volume variables (primary/secondary variables) in the scv
     * \return Amount per volume of the conserved quantities
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        auto storage = ParentType::computeStorage(problem, scv, volVars);

        const auto massOrMoleSolidDensity = [](const auto& volVars, const int sCompIdx)
        { return useMoles ? volVars.solidComponentMolarDensity(sCompIdx) : volVars.solidComponentDensity(sCompIdx); };

        // compute storage term of all solid components
        for (int sCompIdx = 0; sCompIdx < numSolidComps-numInertSolidComps; ++sCompIdx)
        {
            auto eqIdx = Indices::conti0EqIdx + numComponents + sCompIdx;
            storage[eqIdx] += volVars.solidVolumeFraction(sCompIdx)
                             * massOrMoleSolidDensity(volVars, sCompIdx);
        }

        return storage;
    }
};

} // end namespace Dumux

#endif
