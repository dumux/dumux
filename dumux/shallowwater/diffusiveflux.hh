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
 * \ingroup SweModel
 * \brief Defines the diffusive flux for the shallow water Model.
 */
#ifndef DUMUX_SHALLOW_WATER_DIFFUSIVE_FLUX_HH
#define DUMUX_SHALLOW_WATER_DIFFUSIVE_FLUX_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/frictionlaws/swefrictionlaws.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/shallowwater/numericalfluxes/fluxrotation.hh>
#include <dumux/shallowwater/numericalfluxes/letmodel.hh>

namespace Dumux
{

/*!
 * \file
 * \ingroup SweModel
 * \brief Diffusive flux for the shallow water Model.
 */
template<class TypeTag>
class SweDiffusiveFlux
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

    //! Compute the diffusive flux
    static NumEqVector flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf)
    {
        NumEqVector flux(0.0);

        //Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];
        const auto& nxy = scvf.unitOuterNormal();

        //For now assume a constant given background viscosity visc_h = 1.0e-6
        Scalar visc_h_bg = 1.0e-6;

        //For now assume a turbulence model based on just the background eddy viscosity
        Scalar turbulence_model = 2;

        //Get the center coordinates of left and right elements
        Scalar xl = insideScv.center()[0];
        Scalar xr = outsideScv.center()[0];
        Scalar yl = insideScv.center()[1];
        Scalar yr = outsideScv.center()[1];

        //Distance between the two cell centers
        Scalar dist = std::sqrt((xr-xl)*(xr-xl)+(yr-yl)*(yr-yl));

        //Get the velocities u and v inside and outside
        Scalar ul = insideVolVars.getU();
        Scalar ur = outsideVolVars.getU();
        Scalar vl = insideVolVars.getV();
        Scalar vr = outsideVolVars.getV();

        Scalar F_u = 0.0;
        Scalar F_v = 0.0;
        Scalar G_u = 0.0;
        Scalar G_v = 0.0;
        Scalar grad_u_s = 0.0;
        Scalar grad_v_s = 0.0;

        //initialise left and right viscosity
        Scalar viscosity_l = visc_h_bg;
        Scalar viscosity_r = visc_h_bg;

        //Elder scaling factors used for longitudinal and lateral diffusion
        Scalar f_eld_u = 1.0;
        Scalar f_eld_v = 1.0;

        //The water depths in the left and right elements
        Scalar hl = insideVolVars.getH();
        Scalar hr = outsideVolVars.getH();

        if (turbulence_model > 0){

          //Compute dx, dy
          //dx = xr - xl;
          //dy = yr - yl;

          //compute the gradients
          grad_u_s = (ur-ul)/dist;
          grad_v_s = (vr-vl)/dist;

          //constant viscosity
          if(turbulence_model == 1){
            viscosity_l = visc_h_bg;
            viscosity_r = visc_h_bg;
          }

          //Elder turbulence model
          if(turbulence_model == 2){
            viscosity_l = visc_h_bg;
            viscosity_r = visc_h_bg;

            Scalar cf_l = 0.0;
            Scalar cf_r = 0.0;

            Scalar ks_l;
            Scalar ks_r;

            //The roughness heights
            ks_l = problem.spatialParams().ks(element, insideScv);
            ks_r = problem.spatialParams().ks(element, outsideScv);
            //ks_l = 0.05;
            //ks_r = 0.05;
            //std::cout << "ks = " << ks_l << std::endl;

            //Use maximum or average ks?
            //Scalar ks = std::max(ks_l,ks_r);

            //The frictionlaws
            auto frictionlaw_l =  problem.spatialParams().frictionlaw(element, insideScv);
            auto frictionlaw_r =  problem.spatialParams().frictionlaw(element, outsideScv);
            //frictionlaw_l = 3;
            //frictionlaw_r = 3;

            //std::cout << "law = " << frictionlaw_l << std::endl;

            //Gravity
            Scalar grav_l = insideVolVars.getGravity();
            Scalar grav_r = outsideVolVars.getGravity();

            //The friction coefficients
            cf_l = computeUstarH(ks_l, hl, grav_l, frictionlaw_l);
            cf_r = computeUstarH(ks_r, hr, grav_r, frictionlaw_r);
            //cf_l = insideVolVars.getFrictionUstarH();
            //cf_r = outsideVolVars.getFrictionUstarH();

            //std::cout << "cf = " << cf_l << std::endl;

            //Velocity magnitude
            Scalar uvmag_l = std::sqrt(std::pow(ul,2.0) + std::pow(vl,2.0));
            Scalar uvmag_r = std::sqrt(std::pow(ur,2.0) + std::pow(vr,2.0));

            //Add Elder contribution nu_t = (kappa/6)*ustar*h = (kappa/6)*cf*|u|*h
            //kappa/ 6.0 --> 0.41/6 --> ~ 0.067 --> better named as c_eld
            // alpha_t ranges between 0.2 and 1.0 for most applications
            Scalar c_eld = 0.41 / 6.0;
            viscosity_l += c_eld * cf_l * uvmag_l * hl;
            viscosity_r += c_eld * cf_r * uvmag_r * hr;

            //The Elder model differentiates between diffusion in flow direction
            //and in cross-flow direction. This is realized using a factor f_eld,
            //which is computed using a simple average (left and right of the current edge)
            Scalar uav  = 0.5 * (ul+ur);
            Scalar vav  = 0.5 * (vl+vr);
            Scalar uvav = std::max(0.001, 0.5 * (uvmag_l+uvmag_r));
            f_eld_u = std::abs(uav)/uvav;
            f_eld_v = std::abs(vav)/uvav;

            //std::cout << "visc = " << viscosity_l << std::endl;
          }

          //compute the diffusive fluxes (should these depend on the lendth of the edge?)
          F_u = 0.5*(viscosity_l * hl + viscosity_r * hr) * grad_u_s * (1.0+nxy[0]*nxy[0]);
          F_v = 0.5*(viscosity_l * hl + viscosity_r * hr) * grad_v_s * nxy[0]*nxy[1];
          G_u = 0.5*(viscosity_l * hl + viscosity_r * hr) * grad_u_s * nxy[0]*nxy[1];
          G_v = 0.5*(viscosity_l * hl + viscosity_r * hr) * grad_v_s * (1.0+nxy[1]*nxy[1]);

        }

        //std::cout << "F_u " << F_u << " F_v " << F_v << " G_u " << G_u << " G_v " << G_v << std::endl;
        //The actual stress term / diffusion flux for this edge
        //note the use of the Elder coefficient, which is one when not using the Elder turbulence model
        Scalar diff_u = f_eld_u*(F_u + G_u);
        Scalar diff_v = f_eld_v*(F_v + G_v);

        Scalar mobility[3] = {1.0};
        letmobility(hl,hr,0.01,0.0001,mobility);

        //The contribution to the total flux
        flux[0] = 0.0;
        flux[1] = diff_u * scvf.area() * mobility[1];
        flux[2] = diff_v * scvf.area() * mobility[2];

        return flux;
    }
};

} // end namespace Dumux

#endif
