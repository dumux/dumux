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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_SCALAR_CONSERVATION_MODEL_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_SCALAR_CONSERVATION_MODEL_FLUXVARIABLES_HH


#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/upwindscheme.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag, class UpwindScheme = UpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>>>
class NavierStokesScalarConservationModelFluxVariables
: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                           typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr bool enableAdvection = true; // TODO
    using Extrusion = Extrusion_t<GetPropType<TypeTag, Properties::GridGeometry>>;

public:

    /*!
     * \brief Returns the advective flux computed by the respective law.
     */
    template<typename FunctionType>
    Scalar advectiveFlux(const FunctionType& upwindTerm) const
    {
        if constexpr (enableAdvection)
        {
            const auto& scvf = this->scvFace();
            const auto velocity = this->problem().faceVelocity(this->element(), this->fvGeometry(), scvf);
            const Scalar volumeFlux = velocity*scvf.unitOuterNormal()*Extrusion::area(scvf)*extrusionFactor_(this->elemVolVars(), scvf);
            return UpwindScheme::apply(*this, upwindTerm, volumeFlux, 0/*phaseIdx*/);
        }
        else
            return 0.0;
    }

    /*!
     * \brief Returns the conductive flux computed by the respective law.
     * \note This overload is used in models considering local thermal equilibrium
     */
    Scalar heatConductionFlux() const
    {
        if constexpr (GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance())
        {
            using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;
            return HeatConductionType::flux(this->problem(),
                                            this->element(),
                                            this->fvGeometry(),
                                            this->elemVolVars(),
                                            this->scvFace(),
                                            this->elemFluxVarsCache());
        }
        else
            return 0.0;
    }

private:

    template<class ElementVolumeVariables, class SubControlVolumeFace>
    Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf) const
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
    }


};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
