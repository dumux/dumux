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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the n-phase immiscible fully implicit models.
 */
#ifndef DUMUX_IMMISCIBLE_LOCAL_RESIDUAL_HH
#define DUMUX_IMMISCIBLE_LOCAL_RESIDUAL_HH

namespace Dumux
{
/*!
 * \ingroup OnePModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the n-phase immiscible fully implicit models.
 */
template<class TypeTag>
class ImmiscibleLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    // first index for the mass balance
    enum { conti0EqIdx = Indices::conti0EqIdx };

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

public:

    /*!
     * \brief Constructor. Gets the upwind weight.
     */
    ImmiscibleLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        upwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    }

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    PrimaryVariables computeStorage(const SubControlVolume& scv,
                                    const VolumeVariables& volVars) const
    {
        // partial time derivative of the phase mass
        PrimaryVariables storage;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            auto eqIdx = conti0EqIdx + phaseIdx;
            storage[eqIdx] = volVars.porosity()
                             * volVars.density(phaseIdx)
                             * volVars.saturation(phaseIdx);
        }

        return storage;
    }


    /*!
     * \brief Evaluate the mass flux over a face of a sub control volume
     * \param scvf The sub control volume face to compute the flux on
     */
    PrimaryVariables computeFlux(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 const BoundaryTypes& bcTypes,
                                 const FluxVariablesCache& fluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.initAndComputeFluxes(this->problem(),
                                      element,
                                      fvGeometry,
                                      elemVolVars,
                                      scvf,
                                      fluxVarsCache);

        // copy weight to local scope for use in lambda expression
        auto w = upwindWeight_;

        PrimaryVariables flux;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // the upwinding scheme
            auto upwindRule = [w, phaseIdx](const VolumeVariables& up, const VolumeVariables& dn)
                              { return (up.density(phaseIdx)*up.mobility(phaseIdx))*(w)
                                     + (dn.density(phaseIdx)*dn.mobility(phaseIdx))*(1-w); };

            auto eqIdx = conti0EqIdx + phaseIdx;
            flux[eqIdx] = fluxVars.advectiveFlux(phaseIdx, upwindRule);
        }
        return flux;
    }

private:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

    Scalar upwindWeight_;
};

}

#endif
