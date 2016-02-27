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
 *        using the one-phase fully implicit model.
 */
#ifndef DUMUX_1P_LOCAL_RESIDUAL_HH
#define DUMUX_1P_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup OnePModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase fully implicit model.
 */
template<class TypeTag>
class OnePLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    //index of the mass balance equation
    enum {
        conti0EqIdx = Indices::conti0EqIdx //index for the mass balance
    };

public:

    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    OnePLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
        mobilityUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MobilityUpwindWeight);
    }

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the OneP
     *        model.
     *
     * This function should not include the source and sink terms.
     *  \param scv The sub control volume
     *  \param volVars The current or previous volVars
     *
     * The volVars can be different to allow computing the implicit euler time derivative here
     */
    PrimaryVariables computeStorage(const SubControlVolume& scv,
                                    const VolumeVariables& volVars) const
    {
        PrimaryVariables storage(0);

        // partial time derivative of the wetting phase mass
        storage[conti0EqIdx] = volVars.density() * volVars.porosity();

        return storage;
    }


    /*!
     * \brief Evaluate the mass flux over a face of a sub-control
     *        volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param fIdx The index of the SCV face
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    PrimaryVariables computeFlux(const SubControlVolumeFace& scvFace) const
    {
        const auto& fluxVars = this->model_().fluxVars(scvFace);

        auto massWeight = massUpwindWeight_;
        auto mobWeight = mobilityUpwindWeight_;

        auto upwindRule = [massWeight, mobWeight](const VolumeVariables& up, const VolumeVariables& dn)
                          { return (up.density(0)*massWeight + dn.density(0)*(1-massWeight))
                                   *(up.mobility(0)*mobWeight + dn.density(0)*(1-mobWeight)); };

        auto flux = fluxVars.darcyFluxVars().computeFlux(0, upwindRule);

        return flux;
    }

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    Scalar temperature(const PrimaryVariables &priVars)
    { return this->problem_.temperature(); /* constant temperature */ }

private:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

    Scalar massUpwindWeight_;
    Scalar mobilityUpwindWeight_;
};

}

#endif
