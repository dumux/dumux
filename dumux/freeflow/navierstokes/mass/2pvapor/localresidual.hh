// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2PVAPOR_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2PVAPOR_LOCAL_RESIDUAL_HH

#include <dumux/freeflow/navierstokes/mass/2p/localresidual.hh>

namespace Dumux {

/*!
 * \brief Extends the 2p Cahn-Hilliard local residual with a vapor convection-diffusion equation.
 *
 * Added equation (vaporEqIdx):
 *   ∂c_v/∂t + ∇·(u c_v) − ∇·(D_eff ∇c_v) = ṁ
 *
 * where D_eff = D_gas · max(0, (1−φ)/2) confines diffusion to the gas phase.
 */
template<class TypeTag>
class NavierStokesMassTwoPVaporLocalResidual
: public NavierStokesMassTwoPLocalResidual<TypeTag>
{
    using ParentType = NavierStokesMassTwoPLocalResidual<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;
    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        auto storage = ParentType::computeStorage(problem, scv, volVars);
        storage[Indices::vaporEqIdx] = volVars.vaporConcentration();
        return storage;
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        auto flux = ParentType::computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        const auto velocity = problem.faceVelocity(element, fvGeometry, scvf);
        const Scalar volumeFlux = velocity * scvf.unitOuterNormal() * scvf.area();

        // advective flux: upwind vapor concentration
        flux[Indices::vaporEqIdx] = upwindVapor_(elemVolVars, scvf, volumeFlux) * volumeFlux;

        // diffusive flux: −D_eff ∇c_v · n, D_eff phase-weighted to gas region
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradCv(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
            gradCv.axpy(elemVolVars[localDof].vaporConcentration(),
                        fluxVarCache.gradN(localDof.index()));

        const auto phiIn  = elemVolVars[scvf.insideScvIdx()].phaseField();
        const auto phiOut = scvf.boundary()
            ? phiIn
            : elemVolVars[scvf.outsideScvIdx()].phaseField();
        // gas volume fraction: (1−φ)/2, clamped to [0,1]
        const auto gasFraction = std::max(Scalar(0), 0.5*(1.0 - 0.5*(phiIn + phiOut)));
        const Scalar dEff = problem.vaporDiffusivity() * gasFraction;

        flux[Indices::vaporEqIdx] -= dEff * (gradCv * scvf.unitOuterNormal()) * scvf.area();

        return flux;
    }

private:
    Scalar upwindVapor_(const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf,
                        const Scalar flux) const
    {
        using std::signbit;
        return signbit(flux)
            ? elemVolVars[scvf.outsideScvIdx()].vaporConcentration()
            : elemVolVars[scvf.insideScvIdx()].vaporConcentration();
    }
};

} // end namespace Dumux

#endif
