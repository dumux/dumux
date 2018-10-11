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
 * \ingroup TracerModel
 * \brief Element-wise calculation of the local residual for problems
 *        using fully implicit tracer model.
 */
#ifndef DUMUX_TRANSPORT_LOCAL_RESIDUAL_HH
#define DUMUX_TRANSPORT_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup TracerModel
 * \brief Element-wise calculation of the local residual for problems
 *        using fully implicit tracer model.
 *
 */
template<class TypeTag>
class TransportLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    static constexpr int phaseIdx = 0;
    static constexpr int transportEqIdx = ModelTraits::Indices::transportEqIdx; //!< first index for the mass balance

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param problem TODO docme!
     * \param scv The sub control volume
     * \param volVars The primary and secondary varaibles on the scv
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        // formulation with mole balances
        storage[transportEqIdx] = volVars.porosity()
                                 * volVars.density(phaseIdx)
                                 * volVars.saturation(phaseIdx);

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param problem TODO docme!
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
     * \param scvf The sub control volume face
     * \param elemFluxVarsCache The cache related to flux compuation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        NumEqVector flux(0.0);

        // the physical quantities for which we perform upwinding
        auto upwindTerm = [&](const auto& volVars)
                          { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx); };

        flux[transportEqIdx] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);

        return flux;
    }

    /*!
     * \brief TODO docme!
     *
     * \param partialDerivatives TODO docme!
     * \param problem TODO docme!
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param curVolVars TODO docme!
     * \param scv The sub control volume.
     */
    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars,
                               const SubControlVolume& scv) const
    {
        const auto porosity = curVolVars.porosity();

        const auto d_storage = scv.volume()*porosity*curVolVars.density(phaseIdx)/this->timeLoop().timeStepSize();

        partialDerivatives[0][0] += d_storage;
    }
};

} // end namespace Dumux

#endif
