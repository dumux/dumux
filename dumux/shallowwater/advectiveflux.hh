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
#ifndef DUMUX_SHALLOW_WATER_ADVECTIVE_FLUX_HH
#define DUMUX_SHALLOW_WATER_ADVECTIVE_FLUX_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/shallowwater/numericalfluxes/exactriemannsolver.hh>
#include <dumux/shallowwater/numericalfluxes/fluxrotation.hh>
#include <dumux/shallowwater/numericalfluxes/letmodel.hh>

namespace Dumux
{

template<class TypeTag>
class SweAdvectiveFlux
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

  public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod myDiscretizationMethod = DiscretizationMethod::cctpfa;

    //! state the type for the corresponding cache
    //using Cache = SweAdvectiveFluxCache<TypeTag>;

    //! Compute the advective flux
    static NumEqVector flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf)
    {
        NumEqVector flux(0.0);

        //Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto& nxy = scvf.unitOuterNormal();

        Scalar cellStatesLeft[4] = {0.0};
        Scalar cellStatesRight[4] = {0.0};
        cellStatesLeft[0]  = insideVolVars.getH();
        cellStatesRight[0] = outsideVolVars.getH();

        cellStatesLeft[1]  = insideVolVars.getU();
        cellStatesRight[1] = outsideVolVars.getU();

        cellStatesLeft[2]  = insideVolVars.getV();
        cellStatesRight[2] = outsideVolVars.getV();

        cellStatesLeft[3]  = insideVolVars.getBottom();
        cellStatesRight[3] = outsideVolVars.getBottom();

        Scalar hllc_hl = cellStatesLeft[0];
        Scalar hllc_hr = cellStatesRight[0];
        Scalar thetal = cellStatesLeft[3] + cellStatesLeft[0];
        Scalar thetar = cellStatesRight[3] + cellStatesRight[0];

        auto ks_av = std::max(outsideVolVars.getKsH(), insideVolVars.getKsH());
        ks_av = std::max(ks_av,1.0E-9);
        ks_av = std::min(ks_av,0.1);

        Scalar mobility[3] = {1.0};
        letmobility(cellStatesLeft[0],cellStatesRight[0],ks_av,mobility);

        //-------------------------------------------------------------
        //
        // Inner boundary
        //
        //-------------------------------------------------------------
        //Hydrostatic reconstrucion after Audusse
        Scalar dzl = std::max(0.0,cellStatesRight[3] - cellStatesLeft[3]);
        cellStatesLeft[0] = std::max(0.0, hllc_hl - dzl);
        Scalar dzr = std::max(0.0,cellStatesLeft[3] - cellStatesRight[3]);
        cellStatesRight[0] = std::max(0.0, hllc_hr - dzr);

        //------------------ Flux rotation -------
        //make rotation for computing 1d HLLC flux
        stateRotation(nxy,cellStatesLeft);
        stateRotation(nxy,cellStatesRight);

        //---------------------------------
        //compute the HLLC flux
        //---------------------------------
        Scalar riemannFlux[3] = {0.0};
        computeExactRiemann(riemannFlux,cellStatesLeft[0],cellStatesRight[0],
                            cellStatesLeft[1],cellStatesRight[1],
                            cellStatesLeft[2],cellStatesRight[2],insideVolVars.getGravity());

        //---------------------------------
        //end of  HLLC flux
        //---------------------

        //redo rotation
        rotateFluxBack(nxy,riemannFlux);

        //Audusse
        Scalar hgzl = 0.5 * (cellStatesLeft[0] + hllc_hl) *(cellStatesLeft[0] - hllc_hl)  ;
        Scalar hgzr = 0.5 * (cellStatesRight[0] + hllc_hr) * (cellStatesRight[0] - hllc_hr);
        Scalar hdxzl = 0.0;
        Scalar hdxzr = 0.0;
        Scalar hdyzl = 0.0;
        Scalar hdyzr = 0.0;
        hdxzl = insideVolVars.getGravity() * nxy[0] * hgzl;
        hdyzl = insideVolVars.getGravity() * nxy[1] * hgzl;
        hdxzr = insideVolVars.getGravity() * nxy[0] * hgzr;
        hdyzr = insideVolVars.getGravity() * nxy[1] * hgzr;

        flux[0] = riemannFlux[0] * scvf.area() * mobility[0];
        flux[1] = (riemannFlux[1]  - hdxzl) * scvf.area() * mobility[1];
        flux[2] = (riemannFlux[2]  - hdyzl) * scvf.area() * mobility[2];

        return flux;
    }
};

} // end namespace Dumux

#endif
