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
 * \brief This file contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume.
 *
 * This means pressure, concentration and temperature gradients, phase
 * densities at the integration point, etc.
 */
#ifndef DUMUX_MPNC_FLUX_VARIABLES_HH
#define DUMUX_MPNC_FLUX_VARIABLES_HH

#include <dumux/common/spline.hh>

#include "diffusion/fluxvariables.hh"
#include "energy/fluxvariables.hh"
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/porousmediumflow/implicit/forchheimerfluxvariables.hh>

namespace Dumux
{

/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitFluxVariables
 * \brief This class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the MpNc model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class MPNCFluxVariables
    : public GET_PROP_TYPE(TypeTag, BaseFluxVariables)
{
    friend typename GET_PROP_TYPE(TypeTag, BaseFluxVariables); // be friends with base class
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    using BaseFluxVariables = typename GET_PROP_TYPE(TypeTag, BaseFluxVariables);

    enum {dim= GridView::dimension};
    enum {dimWorld= GridView::dimensionworld};
    enum {enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion)};
    enum {enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy)};
    enum {enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic)};
    enum {numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};

    using DimVector = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using FluxVariablesDiffusion = MPNCFluxVariablesDiffusion<TypeTag, enableDiffusion>;
    using FluxVariablesEnergy = MPNCFluxVariablesEnergy<TypeTag, enableEnergy, numEnergyEquations>;

public:
    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        elemVolVarsPtr_ = &elemVolVars;

        // velocities can be obtained from the Parent class.
        BaseFluxVariables::update(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary);

        // update the flux data of the energy module (i.e. isothermal
        // or non-isothermal)
        fluxVarsEnergy_.update(problem, element, fvGeometry, this->face(), *this, elemVolVars);

        // update the flux data of the diffusion module (i.e. with or
        // without diffusion)
        fluxVarsDiffusion_.update(problem, element, fvGeometry, this->face(), elemVolVars);

        extrusionFactor_ =
            (elemVolVars[this->face().i].extrusionFactor()
             + elemVolVars[this->face().j].extrusionFactor()) / 2;
    }

    /*!
     * \brief Returns a reference to the volume
     *        variables of the i-th sub-control volume of the current
     *        element.
     */
    const VolumeVariables &volVars(const unsigned int idx) const
    { return (*elemVolVarsPtr_)[idx]; }

    /*!
     * \brief Returns th extrusion factor for the sub-control volume face
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    ////////////////////////////////////////////////
    // forward calls to the diffusion module
    Scalar porousDiffCoeffL(const unsigned int compIdx) const
    { return fluxVarsDiffusion_.porousDiffCoeffL(compIdx); }

    Scalar porousDiffCoeffG(const unsigned int compIIdx,
                            const unsigned int compJIdx) const
    { return fluxVarsDiffusion_.porousDiffCoeffG(compIIdx, compJIdx); }

    const Scalar moleFraction(const unsigned int phaseIdx,
                              const unsigned int compIdx) const
    { return fluxVarsDiffusion_.moleFraction(phaseIdx, compIdx); }

    const GlobalPosition &moleFractionGrad(const unsigned int phaseIdx,
                                      const unsigned int compIdx) const
    { return fluxVarsDiffusion_.moleFractionGrad(phaseIdx, compIdx); }

    // end of forward calls to the diffusion module
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // forward calls to the temperature module
    const GlobalPosition &temperatureGrad() const
    { return fluxVarsEnergy_.temperatureGrad(); }

    const FluxVariablesEnergy &fluxVarsEnergy() const
    { return fluxVarsEnergy_; }
    // end of forward calls to the temperature module
    ////////////////////////////////////////////////

private:
    const ElementVolumeVariables* elemVolVarsPtr_;

    // The extrusion factor for the sub-control volume face
    Scalar extrusionFactor_;

    FluxVariablesDiffusion  fluxVarsDiffusion_;
    FluxVariablesEnergy     fluxVarsEnergy_;
};

} // end namespace

#endif
