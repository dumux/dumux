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
 * \copydoc Dumux::NavierStokesResidualImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1P_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MASS_1P_LOCAL_RESIDUAL_HH

#include <dumux/assembly/cclocalresidual.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/freeflow/navierstokes/energy/localresidual.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux {


/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using the staggered discretization
 */
template<class TypeTag>
class NavierStokesMassOnePLocalResidual : public CCLocalResidual<TypeTag>
{
    using ParentType = CCLocalResidual<TypeTag>;
    friend class CCLocalResidual<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using EnergyLocalResidual = NavierStokesEnergyLocalResidual<ModelTraits>;
    //using PhasefieldLocalResidual = NavierStokesPhasefieldLocalResidual<ModelTraits>;

    static_assert(GridGeometry::discMethod == DiscretizationMethod::cctpfa);

    static constexpr int phi1Idx = Indices::phi1Idx;
    static constexpr int phi2Idx = Indices::phi2Idx;
    static constexpr int phi3Idx = Indices::phi3Idx;
    static constexpr int u1Idx = Indices::u1Idx;
    static constexpr int u2Idx = Indices::u2Idx;
    static constexpr int u3Idx = Indices::u3Idx;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage term of the equation
     *
     * \param problem The problem to solve
     * \param scv The sub-control volume over which we integrate the storage term
     * \param volVars The volume variables associated with the scv
     * \note has to be implemented by the model specific residual class
     *
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);
        const Scalar xi = getParam<Scalar>("Phasefield.xi");
        const Scalar delta = getParam<Scalar>("Phasefield.delta");
        const Scalar mineralDMolarDensity = getParam<Scalar>("Phasefield.MineralDMolarDensity");
        const Scalar mineralPMolarDensity = getParam<Scalar>("Phasefield.MineralPMolarDensity");
        //storage[Indices::conti0EqIdx] = volVars.density() * (volVars.phasefield(1) + delta)
        //    + volVars.phasefield(1);
        storage[Indices::conti0EqIdx] = 0.0;//volVars.density();

        //! The energy storage in the fluid phase
        EnergyLocalResidual::fluidPhaseStorage(storage, volVars);
        //PhasefieldLocalResidual::fluidPhaseStorage(storage, volVars);
        const auto& priVars = volVars.priVars();
        const static std::array<Scalar, 3> b_D = {-1.0, 1.0, 0.0};
        const static std::array<Scalar, 3> b_P = {0.0, 1.0, 1.0};
        storage[Indices::phasefield1EqIdx] = xi * xi * priVars[phi1Idx];
        storage[Indices::phasefield2EqIdx] = xi * xi * priVars[phi2Idx];
        storage[Indices::phasefield3EqIdx] = xi * xi * priVars[phi3Idx];
        storage[Indices::u1TransportEqIdx] = (priVars[phi1Idx] + delta) * priVars[u1Idx]
            + priVars[phi2Idx] * b_D[0] * mineralDMolarDensity
            + priVars[phi3Idx] * b_P[0] * mineralPMolarDensity;
            ;
        storage[Indices::u2TransportEqIdx] = (priVars[phi1Idx] + delta) * priVars[u2Idx]
            + priVars[phi2Idx] * b_D[1] * mineralDMolarDensity
            + priVars[phi3Idx] * b_P[1] * mineralPMolarDensity;
            ;
        storage[Indices::u3TransportEqIdx] = (priVars[phi1Idx] + delta) * priVars[u3Idx]
            + priVars[phi2Idx] * b_D[2] * mineralDMolarDensity
            + priVars[phi3Idx] * b_P[2] * mineralPMolarDensity;
            ;

        return storage;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source;
        const auto& priVars = elemVolVars[scv].priVars();
        const static Scalar sigma = getParam<Scalar>("Phasefield.sigma");
        const static Scalar xi = getParam<Scalar>("Phasefield.xi");
        const static Scalar react = getParam<Scalar>("Phasefield.Reaction");
        Scalar f_P = react * (priVars[u2Idx]*priVars[u3Idx] - 1.0);
        //source[Indices::phasefield1EqIdx] =
        //    // - gamma P'
        //    - 16.0 * sigma * (
        //    priVars[phi1Idx] * priVars[phi2Idx] * (priVars[phi2Idx] - priVars[phi1Idx])
        //    )
        //    - 4.0 * xi * priVars[phi1Idx] * priVars[phi2Idx] * f_P
        //    ;
        //source[Indices::phasefield2EqIdx] =
        //    // - gamma P'
        //    - 16.0 * sigma * (
        //    priVars[phi2Idx] * priVars[phi1Idx] * (priVars[phi1Idx] - priVars[phi2Idx])
        //    )
        //    + 4.0 * xi * priVars[phi1Idx] * priVars[phi2Idx] * f_P
        Scalar f_D = react * (priVars[u2Idx] / priVars[u1Idx] - 1.0);
        source[Indices::phasefield1EqIdx] =
            16.0/3.0 * sigma * (
                 -2*(priVars[phi1Idx] -3*priVars[phi1Idx]*priVars[phi1Idx]
                     +2*priVars[phi1Idx]*priVars[phi1Idx]*priVars[phi1Idx])
                 +  (priVars[phi2Idx] -3*priVars[phi2Idx]*priVars[phi2Idx]
                     +2*priVars[phi2Idx]*priVars[phi2Idx]*priVars[phi2Idx])
                 +  (priVars[phi3Idx] -3*priVars[phi3Idx]*priVars[phi3Idx]
                     +2*priVars[phi3Idx]*priVars[phi3Idx]*priVars[phi3Idx])
                 )
            -4.0 * xi * priVars[phi1Idx] * priVars[phi2Idx] * f_D
            -4.0 * xi * priVars[phi1Idx] * priVars[phi3Idx] * f_P;
            //-4.0 * xi * std::max(0.0, priVars[phi1Idx] * priVars[phi2Idx]) * f_D
            //-4.0 * xi * std::max(0.0, priVars[phi1Idx] * priVars[phi3Idx]) * f_P;
            ;
        source[Indices::phasefield2EqIdx] =
            16.0/3.0 * sigma * (
                    (priVars[phi1Idx] -3*priVars[phi1Idx]*priVars[phi1Idx]
                     +2*priVars[phi1Idx]*priVars[phi1Idx]*priVars[phi1Idx])
                -2* (priVars[phi2Idx] -3*priVars[phi2Idx]*priVars[phi2Idx]
                     +2*priVars[phi2Idx]*priVars[phi2Idx]*priVars[phi2Idx])
                +   (priVars[phi3Idx] -3*priVars[phi3Idx]*priVars[phi3Idx]
                     +2*priVars[phi3Idx]*priVars[phi3Idx]*priVars[phi3Idx])
                )
            +4.0 * xi * priVars[phi1Idx] * priVars[phi2Idx] * f_D;
            //+4.0 * xi * std::max(0.0, priVars[phi1Idx] * priVars[phi2Idx]) * f_D;
            ;
        source[Indices::phasefield3EqIdx] =
            16.0/3.0 * sigma * (
                    (priVars[phi1Idx] -3*priVars[phi1Idx]*priVars[phi1Idx]
                     +2*priVars[phi1Idx]*priVars[phi1Idx]*priVars[phi1Idx])
                +   (priVars[phi2Idx] -3*priVars[phi2Idx]*priVars[phi2Idx]
                     +2*priVars[phi2Idx]*priVars[phi2Idx]*priVars[phi2Idx])
                -2* (priVars[phi3Idx] -3*priVars[phi3Idx]*priVars[phi3Idx]
                     +2*priVars[phi3Idx]*priVars[phi3Idx]*priVars[phi3Idx])
                )
            +4.0 * xi * priVars[phi1Idx] * priVars[phi3Idx] * f_P;
            //+4.0 * xi * std::max(0.0, priVars[phi1Idx] * priVars[phi3Idx]) * f_P;
            ;

        return source;
    }


    /*!
     * \brief Evaluatex the mass flux over a face of a sub control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
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

        NumEqVector flux;

        // the physical quantities for which we perform upwinding
        auto upwindTerm = [](const auto& volVars) { return volVars.density(); };
        flux[Indices::conti0EqIdx] = fluxVars.advectiveFlux(upwindTerm);

        //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConvectionFlux(flux, fluxVars);

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        //! Add phasefield fluxes. For sharp-interface model the contribution is zero.
        //PhasefieldLocalResidual::phasefieldFlux(flux, fluxVars);
        diffusiveFlux(flux, fluxVars, fvGeometry, elemVolVars, scvf);
        advectiveFlux(flux, fluxVars, fvGeometry, elemVolVars, scvf);

        return flux;
    }

    void advectiveFlux(NumEqVector& flux, FluxVariables& fluxVars, FVElementGeometry fvGeometry,
            ElementVolumeVariables elemVolVars, SubControlVolumeFace scvf) const
    {
        const static Scalar delta = getParam<Scalar>("Phasefield.delta");

        auto upwindTerm1 = [](const auto& volVars) { return (volVars.priVar(phi1Idx) + delta) *
            volVars.priVar(u1Idx); };
        auto upwindTerm2 = [](const auto& volVars) { return (volVars.priVar(phi1Idx) + delta) *
            volVars.priVar(u2Idx); };
        auto upwindTerm3 = [](const auto& volVars) { return (volVars.priVar(phi1Idx) + delta) *
            volVars.priVar(u3Idx); };
        flux[Indices::u1TransportEqIdx] += fluxVars.advectiveFlux(upwindTerm1);
        flux[Indices::u2TransportEqIdx] += fluxVars.advectiveFlux(upwindTerm2);
        flux[Indices::u3TransportEqIdx] += fluxVars.advectiveFlux(upwindTerm3);
    }

    void diffusiveFlux(NumEqVector& flux, FluxVariables& fluxVars, FVElementGeometry fvGeometry,
            ElementVolumeVariables elemVolVars, SubControlVolumeFace scvf) const
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const static Scalar sigma = getParam<Scalar>("Phasefield.sigma");
        //const static Scalar sigma = problem.getSigma();
        const static Scalar xi = getParam<Scalar>("Phasefield.xi");
        const static Scalar delta = getParam<Scalar>("Phasefield.delta");
        const static Scalar D_u = getParam<Scalar>("Phasefield.DiffCoeff");
        const int numDiffusion = 6;
        const static std::array<Scalar, numDiffusion> diffCoeff = {
            xi*xi*sigma, xi*xi*sigma, xi*xi*sigma,// p1 - p3
        //    xi*xi*sigma, xi*xi*sigma,// p11, p12
            (insideVolVars.priVar(phi1Idx)+delta) * D_u,//u_A
            (insideVolVars.priVar(phi1Idx)+delta) * D_u,
            (insideVolVars.priVar(phi1Idx)+delta) * D_u
        //    (insideVolVars.priVar(Indices::p1Idx)+delta) * D_u,//u_3A
        //    (insideVolVars.priVar(Indices::p1Idx)+delta) * D_u,
        //    (insideVolVars.priVar(Indices::p1Idx)+delta) * D_u
        };

        for (int i = 0; i < numDiffusion; ++i)
        {
            int k = i + phi1Idx;
            const auto valInside = insideVolVars.priVar(k);
            const auto valOutside = outsideVolVars.priVar(k);

            Scalar tij;

            const Scalar ti = computeTpfaTransmissibility(scvf, insideScv, diffCoeff[i],
                                                          insideVolVars.extrusionFactor());
            //auto distanceVector = scvf.ipGlobal();
            //distanceVector -= insideScv.center();
            //distanceVector /= distanceVector.two_norm2();
            //const Scalar ti = diffCoeff[i] * insideVolVars.extrusionFactor() * (distanceVector *
            //        scvf.unitOuterNormal());

            if (scvf.boundary())
                tij = scvf.area()*ti;
            else
            {
                const auto outsideScvIdx = scvf.outsideScvIdx();
                const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
                const Scalar tj = -computeTpfaTransmissibility(scvf, outsideScv, diffCoeff[i],
                                                               outsideVolVars.extrusionFactor());
                //distanceVector = scvf.ipGlobal();
                //distanceVector -= outsideScv.center();
                //distanceVector /= distanceVector.two_norm2();
                //const Scalar tj = -diffCoeff[i] * outsideVolVars.extrusionFactor() * (distanceVector *
                //        scvf.unitOuterNormal());

                if(ti*tj < 0)
                    tij = 0;
                else
                    tij = scvf.area() * (ti * tj)/(ti + tj);
            }
        flux[k] += tij * (valInside - valOutside);
        }

    }

    private:

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif   // DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH
