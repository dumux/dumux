// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassTwoPLocalResidual
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2P_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2P_LOCAL_RESIDUAL_HH

#include <type_traits>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/common/typetraits/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for two-phase flow.
 */
template<class TypeTag>
class NavierStokesMassTwoPLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
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
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage term of the equation
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);
        storage[Indices::conti0EqIdx] = volVars.pressure()*1e-15;
        storage[Indices::phaseFieldEqIdx] = volVars.phaseField();
        storage[Indices::chemicalPotentialEqIdx] = 0.0;
        return storage;
    }

    /*!
     * \brief Evaluate the mass flux over a face of a sub control volume.
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
        NumEqVector flux(0.0);

        // TODO variable extrusion factor?
        const auto velocity = problem.faceVelocity(element, fvGeometry, scvf);
        const Scalar volumeFlux = velocity*scvf.unitOuterNormal()
            *Extrusion::area(fvGeometry, scvf)*elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        flux[Indices::conti0EqIdx] = volumeFlux;

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradPhaseField(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradChemicalPotential(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            // v.axpy(a, w) means v += a*w
            gradPhaseField.axpy(
                volVars.phaseField(),
                fluxVarCache.gradN(scv.indexInElement())
            );
            gradChemicalPotential.axpy(
                volVars.chemicalPotential(),
                fluxVarCache.gradN(scv.indexInElement())
            );
        }

        // advection of the phase field
        flux[Indices::phaseFieldEqIdx] = upwindSchemeMultiplier_(
            elemVolVars, scvf, volumeFlux,
            [](const auto& volVars) { return volVars.phaseField(); }
        )*volumeFlux;

        flux[Indices::phaseFieldEqIdx] += -1.0*vtmv(
            scvf.unitOuterNormal(), problem.mobility(), gradChemicalPotential
        )*scvf.area();

        flux[Indices::chemicalPotentialEqIdx] = -1.0*vtmv(
            scvf.unitOuterNormal(), problem.surfaceTension(), gradPhaseField
        )*scvf.area();

        return flux;
    }

    /*!
     * \brief Compute the source term of the equation
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
     * \param scv The sub control volume
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        // add contributions from problem (e.g. double well potential)
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        return source;
    }

private:

    template<class ElemVolVars, class SubControlVolumeFace, class UpwindTermFunction, class Scalar>
    Scalar upwindSchemeMultiplier_(const ElemVolVars& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   const Scalar flux,
                                   const UpwindTermFunction& upwindTerm) const
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        constexpr Scalar upwindWeight = 1.0; // TODO: pass this from outside?

        using std::signbit;
        if (signbit(flux)) // if sign of flux is negative
            return (upwindWeight*upwindTerm(outsideVolVars)
                    + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
        else
            return (upwindWeight*upwindTerm(insideVolVars)
                    + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
    }

};

} // end namespace Dumux

#endif
