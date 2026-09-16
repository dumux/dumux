// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the residual for problems
 *        using the n-phase immiscible fully implicit models.
 */
#ifndef DUMUX_IMMISCIBLE_LOCAL_RESIDUAL__HH
#define DUMUX_IMMISCIBLE_LOCAL_RESIDUAL__HH


#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/upwindscheme.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the residual for problems
 *        using the n-phase immiscible fully implicit models.
 * \tparam TypeTag The TypeTag
 */
template<class TypeTag>
class ImmiscibleLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;

    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using SubControlVolume = typename ElementDiscretization::SubControlVolume;
    using SubControlVolumeFace = typename ElementDiscretization::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridDiscretization>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int conti0EqIdx = ModelTraits::Indices::conti0EqIdx; //!< first index for the mass balance

    static_assert(!ModelTraits::enableEnergyBalance(), "The energy balance is not implemented yet");

public:
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage integral over a sub control volume
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub control volume
     * \param isPreviousTimeLevel If set to true, the storage term is evaluated on the previous time level
     */
    NumEqVector storageIntegral(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        const auto& vars = elemVars[scv];

        // partial time derivative of the phase mass, mass lumped
        NumEqVector storage(0.0);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            storage[conti0EqIdx + phaseIdx] = vars.porosity()
                                              * vars.density(phaseIdx)
                                              * vars.saturation(phaseIdx);

        storage *= Extrusion::volume(elemDisc, scv) * vars.extrusionFactor();

        return storage;
    }

    /*!
     * \brief Calculate the mass flux integral over a sub control volume face
     *
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     * \param scvf The sub control volume face
     */
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

        const auto& insideVars = elemVars[elemDisc.scv(scvf.insideScvIdx())];
        const auto& outsideVars = elemVars[elemDisc.scv(scvf.outsideScvIdx())];

        // The mobility is a control volume quantity here, so it is upwinded on the
        // integrated flux and the tensor carries the permeability alone.
        const auto tensor = [] (const auto& vars)
        {
            auto permeability = vars.permeability();
            permeability *= vars.extrusionFactor();
            return permeability;
        };

        NumEqVector flux(0.0);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar volumeFlux = 0.0;
            for (const auto& qpData : Dumux::CVFE::quadratureRule(elemDisc, scvf))
            {
                const auto& faceIpData = qpData.ipData();
                volumeFlux += qpData.weight()*(faceIpData.unitOuterNormal()*AdvectionType::fluxTerm(
                    problem, elemDisc.element(), elemDisc, elemVars, phaseIdx, faceIpData, tensor
                ));
            }

            // the physical quantity for which we perform upwinding
            const auto upwindTerm = [phaseIdx] (const auto& vars)
                                    { return vars.density(phaseIdx)*vars.mobility(phaseIdx); };

            flux[conti0EqIdx + phaseIdx] = volumeFlux*Dumux::Detail::upwindSchemeMultiplier(
                insideVars, outsideVars, upwindTerm, volumeFlux, phaseIdx, upwindWeight
            );
        }

        return flux;
    }
};

} // end namespace Dumux::Experimental

#endif
