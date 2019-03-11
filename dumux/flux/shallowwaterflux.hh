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
#ifndef DUMUX_FLUX_SHALLOW_WATER_FLUX_HH
#define DUMUX_FLUX_SHALLOW_WATER_FLUX_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>

namespace Dumux{

template<class TypeTag>
class ShallowWaterFlux
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    public:

    /*!
     * \ingroup Flux
     * \brief Prepares the riemann problem for the advective flux for
     *        the shallow water model. The actual model uses an exact
     *        Riemann solver after Torro and the reconstruction after
     *        Audusse and a flux limiter for small water depths.
     *
     * \todo The choice of the Riemann solver should be more flexible
     */

    using Cache = FluxVariablesCaching::EmptyAdvectionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;

    static NumEqVector flux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf)
    {

        //Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto& nxy = scvf.unitOuterNormal();

        auto riemannFlux = ShallowWater::riemannProblem(insideVolVars.waterDepth(),
                                                        outsideVolVars.waterDepth(),
                                                        insideVolVars.velocity(0),
                                                        outsideVolVars.velocity(0),
                                                        insideVolVars.velocity(1),
                                                        outsideVolVars.velocity(1),
                                                        insideVolVars.bedSurface(),
                                                        outsideVolVars.bedSurface(),
                                                        insideVolVars.gravity(),
                                                        nxy);

        NumEqVector localFlux(0.0);
        localFlux[0] = riemannFlux[0] * scvf.area();
        localFlux[1] = riemannFlux[1] * scvf.area();
        localFlux[2] = riemannFlux[2] * scvf.area();

        return localFlux;
    }
};

} // end namespace Dumux

#endif
