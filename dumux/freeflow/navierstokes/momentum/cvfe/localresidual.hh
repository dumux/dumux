// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesResidualImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/cvfelocalresidual.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/flux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using the CVFE discretizations
 */
template<class TypeTag>
class NavierStokesMomentumCVFELocalResidual
: public CVFELocalResidual<TypeTag>
{
    using ParentType = CVFELocalResidual<TypeTag>;

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
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr auto dim = GridView::dimension;

    using FluxContext = NavierStokesMomentumFluxContext<Problem, FVElementGeometry, ElementVolumeVariables, ElementFluxVariablesCache>;
    using FluxHelper = NavierStokesMomentumFluxCVFE<GridGeometry, NumEqVector>;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param fvGeometry The finite-volume geometry of the element
     * \param scv The sub-control volume over which we integrate the storage term
     * \param volVars The volume variables associated with the scv
     * \param isPreviousStorage If set to true, the storage term is evaluated on the previous time level.
     *
     */
    NumEqVector computeStorage(const Problem& problem,
                               const FVElementGeometry& fvGeometry,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars,
                               const bool isPreviousStorage) const
    {
        return problem.density(fvGeometry.element(), fvGeometry, scv, isPreviousStorage) * volVars.velocity();
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVolVars The volume variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     *
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source = ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);

        // add rho*g (note that gravity might be zero in case it's disabled in the problem)
        source +=  problem.density(element, fvGeometry, scv) * problem.gravity();

        // Axisymmetric problems in 2D feature an extra source terms arising from the transformation to cylindrical coordinates.
        // See Ferziger/Peric: Computational methods for Fluid Dynamics (2020)
        // https://doi.org/10.1007/978-3-319-99693-6
        // Chapter 9.9 and Eq. (9.81) and comment on finite volume methods
        if constexpr (dim == 2 && isRotationalExtrusion<Extrusion>)
        {
            // the radius with respect to the rotation axis
            const auto r = scv.center()[Extrusion::radialAxis] - fvGeometry.gridGeometry().bBoxMin()[Extrusion::radialAxis];

            // The velocity term is new with respect to Cartesian coordinates and handled below as a source term
            // it only enters the balance of the momentum balance in radial direction
            source[Extrusion::radialAxis] += -2.0*problem.effectiveViscosity(element, fvGeometry, scv)
                * elemVolVars[scv].velocity(Extrusion::radialAxis) / (r*r);

            // Pressure term (needed because we incorporate pressure in terms of a surface integral).
            // grad(p) becomes div(pI) + (p/r)*n_r in cylindrical coordinates. The second term
            // is new with respect to Cartesian coordinates and handled below as a source term.
            source[Extrusion::radialAxis] += problem.pressure(element, fvGeometry, scv)/r;
        }

        return source;
    }

        /*!
     * \brief Evaluates the mass flux over a face of a sub control volume.
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
        FluxContext context(problem, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
        FluxHelper fluxHelper;

        NumEqVector flux(0.0);
        flux += fluxHelper.advectiveMomentumFlux(context);
        flux += fluxHelper.diffusiveMomentumFlux(context);
        flux += fluxHelper.pressureContribution(context);
        return flux;
    }
};

} // end namespace Dumux

#endif
