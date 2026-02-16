// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMomentumCVFELocalResidual
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/boundaryflag.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/assembly/cvfelocalresidual.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/flux.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/felocalresidual.hh>

namespace Dumux {

namespace Detail {

//! helper struct detecting if a problem has new source interface
template<class P, class FVG, class EV, class IPD>
using SourceWithIpDataInterface = decltype(
    std::declval<P>().source(std::declval<FVG>(), std::declval<EV>(), std::declval<IPD>())
);

template<class P, class FVG, class EV, class IPD>
constexpr inline bool hasProblemSourceWithIpDataInterface()
{ return Dune::Std::is_detected<SourceWithIpDataInterface, P, FVG, EV, IPD>::value; }

} // end namespace Detail

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using CVFE discretizations
 */
template<class TypeTag>
class NavierStokesMomentumCVFELocalResidual
: public CVFELocalResidual<TypeTag>
{
    using ParentType = CVFELocalResidual<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVariablesCache = typename GridVariables::GridVolumeVariables;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using Variables = typename GridVariablesCache::VolumeVariables;

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
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr auto dim = GridView::dimension;

    using LocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using BaseIpData = Dumux::CVFE::InterpolationPointData<typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate, GlobalPosition>;
    using FluxContext = NavierStokesMomentumFluxContext<Problem, FVElementGeometry, ElementVariables, ElementFluxVariablesCache>;
    using FluxHelper = NavierStokesMomentumFluxCVFE<GridGeometry, NumEqVector>;
    using FeResidual = NavierStokesMomentumFELocalResidual<Scalar, NumEqVector, LocalBasis, Extrusion>;

    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;
    using FluxFunctionContext = NavierStokesMomentumFluxFunctionContext<Problem, FVElementGeometry, ElementVariables, FluxVariablesCache>;
    using FluxFunctionHelper = NavierStokesMomentumFluxFunctionCVFE<GridGeometry, NumEqVector>;

public:
    //! Use the parent type's constructor
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage term of the equation
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
                               const Variables& volVars,
                               const bool isPreviousStorage) const
    {
        return problem.density(fvGeometry.element(), fvGeometry, ipData(fvGeometry, scv), isPreviousStorage) * volVars.velocity();
    }

    /*!
     * \brief Calculate the storage integral
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub control volume
     * \param isPreviousTimeLevel If set to true, the storage term is evaluated on the previous time level.
     *
     */
    NumEqVector storageIntegral(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        const auto& volVars = elemVars[scv];
        // We apply mass lumping
        NumEqVector storage = this->asImp().problem().density(fvGeometry.element(), fvGeometry, ipData(fvGeometry, scv), isPreviousTimeLevel)
                            * volVars.velocity();

        storage *= Extrusion::volume(fvGeometry, scv) * volVars.extrusionFactor();

        return storage;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub-control volume over which we integrate the source term
     *
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVariables& elemVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source;

        if constexpr (Detail::hasProblemSourceWithIpDataInterface<Problem, FVElementGeometry, ElementVariables, BaseIpData>())
        {
            source = problem.source(fvGeometry, elemVars, ipData(fvGeometry, scv.center()));

            // ToDo: point source data with ipData
            // add contribution from possible point sources
            if (!problem.pointSourceMap().empty())
                source += problem.scvPointSources(element, fvGeometry, elemVars, scv);
        }
        else
            source = ParentType::computeSource(problem, element, fvGeometry, elemVars, scv);


        // add rho*g (note that gravity might be zero in case it's disabled in the problem)
        const auto& data = ipData(fvGeometry, scv);
        source +=  problem.density(element, fvGeometry, data) * problem.gravity();

        // Axisymmetric problems in 2D feature an extra source term arising from the transformation to cylindrical coordinates.
        // See Ferziger/Peric: Computational methods for Fluid Dynamics (2020)
        // https://doi.org/10.1007/978-3-319-99693-6
        // Chapter 9.9 and Eq. (9.81) and comment on finite volume methods
        if constexpr (dim == 2 && isRotationalExtrusion<Extrusion>)
        {
            // the radius with respect to the rotation axis
            const auto r = scv.center()[Extrusion::radialAxis] - fvGeometry.gridGeometry().bBoxMin()[Extrusion::radialAxis];

            // The velocity term is new with respect to Cartesian coordinates and handled below as a source term
            // It only enters the balance of the momentum balance in radial direction
            source[Extrusion::radialAxis] += -2.0*problem.effectiveViscosity(element, fvGeometry, data)
                * elemVars[scv].velocity(Extrusion::radialAxis) / (r*r);

            // Pressure term (needed because we incorporate pressure in terms of a surface integral).
            // grad(p) becomes div(pI) + (p/r)*n_r in cylindrical coordinates. The second term
            // is new with respect to Cartesian coordinates and handled below as a source term.
            source[Extrusion::radialAxis] += problem.pressure(element, fvGeometry, data)/r;
        }

        return source;
    }

    /*!
     * \brief Calculate the source integral
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub control volume
     *
     */
    NumEqVector sourceIntegral(const FVElementGeometry& fvGeometry,
                               const ElementVariables& elemVars,
                               const SubControlVolume& scv) const
    {
        static_assert(!(dim == 2 && isRotationalExtrusion<Extrusion>), "Rotational extrusion source terms are not implemented for integral interface.");

        const auto& problem = this->asImp().problem();

        NumEqVector source(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv))
        {
            source += qpData.weight() * (problem.source(fvGeometry, elemVars, qpData.ipData())
                                        + problem.density(fvGeometry.element(), fvGeometry, qpData.ipData()) * problem.gravity());
        }

        // ToDo: point source data with ipData
        // add contribution from possible point sources
        if (!problem.pointSourceMap().empty())
            source += Extrusion::volume(fvGeometry, scv) * problem.scvPointSources(fvGeometry.element(), fvGeometry, elemVars, scv);

        source *= elemVars[scv].extrusionFactor();

        return source;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a sub control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVars The variables for all local dofs of the element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVariables& elemVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxContext context(problem, fvGeometry, elemVars, elemFluxVarsCache, scvf);
        FluxHelper fluxHelper;

        NumEqVector flux(0.0);
        flux += fluxHelper.advectiveMomentumFlux(context);
        flux += fluxHelper.diffusiveMomentumFlux(context);
        flux += fluxHelper.pressureContribution(context);
        return flux;
    }

    /*!
     * \brief Calculates the flux integral over a sub control volume face.
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param elemFluxVarsCache the flux variable caches for the element's flux stencils
     * \param scvf The sub control volume face
     *
     */
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const ElementFluxVariablesCache& elemFluxVarsCache,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();

        NumEqVector flux(0.0);
        GlobalPosition velIntegral(0.0);
        FluxFunctionHelper fluxFunctionHelper;

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& fluxVarsCache = elemFluxVarsCache[qpData.ipData()];
            FluxFunctionContext context(this->problem(), fvGeometry, elemVars, fluxVarsCache);

            velIntegral += context.velocity() * qpData.weight();
            flux += qpData.weight() * ( fluxFunctionHelper.diffusiveMomentumFluxIntegrand(context, qpData.ipData())
                                      + fluxFunctionHelper.pressureFluxIntegrand(context, qpData.ipData()) );
        }
        flux += fluxFunctionHelper.advectiveMomentumFluxIntegral(problem, fvGeometry, elemVars, elemFluxVarsCache, scvf, velIntegral);

        flux *= elemVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor();

        return flux;
    }

    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVariables& prevElemVolVars,
                                     const ElementVariables& curElemVolVars) const
    {
        FeResidual::addStorageTerms(
            residual, problem, fvGeometry, prevElemVolVars, curElemVolVars, this->timeLoop().timeStepSize()
        );
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars,
                                           const ElementFluxVariablesCache& elemFluxVarsCache,
                                           const ElementBoundaryTypes &elemBcTypes) const
    {
        FeResidual::addFluxAndSourceTerms(
            residual, problem, fvGeometry, elemVars, elemFluxVarsCache, elemBcTypes
        );
    }
};

} // end namespace Dumux

#endif
