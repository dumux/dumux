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
#include <dumux/discretization/fem/integrationpointdata.hh>
#include <dumux/discretization/cvfe/integrationpointdata.hh>
#include <dumux/assembly/cvfelocalresidual.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/flux.hh>

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
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr auto dim = GridView::dimension;

    using LocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using BaseIpData = Dumux::CVFE::IntegrationPointData<typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate, GlobalPosition>;
    using IpData = FEIntegrationPointData<GlobalPosition, LocalBasis>;
    using BoundaryFlag = Dumux::BoundaryFlag<typename GridView::Grid>;
    using FaceIpData = FEFaceIntegrationPointData<GlobalPosition, LocalBasis, BoundaryFlag>;
    using FluxContext = NavierStokesMomentumFluxContext<Problem, FVElementGeometry, ElementVolumeVariables, ElementFluxVariablesCache>;
    using FluxFunctionContext = NavierStokesMomentumFluxFunctionContext<Problem, FVElementGeometry, ElementVolumeVariables, IpData>;
    using FluxHelper = NavierStokesMomentumFluxCVFE<GridGeometry, NumEqVector>;

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
                               const VolumeVariables& volVars,
                               const bool isPreviousStorage) const
    {
        return problem.density(fvGeometry.element(), fvGeometry, ipData(fvGeometry, scv), isPreviousStorage) * volVars.velocity();
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
        NumEqVector source;

        if constexpr (Detail::hasProblemSourceWithIpDataInterface<Problem, FVElementGeometry, ElementVolumeVariables, BaseIpData>())
        {
            source = problem.source(fvGeometry, elemVolVars, ipData(fvGeometry, scv.center()));

            // ToDo: point source data with ipData
            // add contribution from possible point sources
            if (!problem.pointSourceMap().empty())
                source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);
        }
        else
            source = ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);


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
                * elemVolVars[scv].velocity(Extrusion::radialAxis) / (r*r);

            // Pressure term (needed because we incorporate pressure in terms of a surface integral).
            // grad(p) becomes div(pI) + (p/r)*n_r in cylindrical coordinates. The second term
            // is new with respect to Cartesian coordinates and handled below as a source term.
            source[Extrusion::radialAxis] += problem.pressure(element, fvGeometry, data)/r;
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

    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& prevElemVolVars,
                                     const ElementVolumeVariables& curElemVolVars) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {
            // Make sure we don't iterate over quadrature points if there are no hybrid dofs
            if( nonCVLocalDofs(fvGeometry).empty() )
                return;

            static const auto intOrder
                = getParamFromGroup<int>(problem.paramGroup(), "Assembly.FEIntegrationOrderStorage", 4);

            const auto &localBasis = fvGeometry.feLocalBasis();
            using RangeType = typename LocalBasis::Traits::RangeType;
            std::vector<RangeType> integralShapeFunctions(localBasis.size(), RangeType(0.0));

            // We apply mass lumping such that we only need to calculate the integral of basis functions
            const auto& geometry = fvGeometry.elementGeometry();
            const auto& quadRule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), intOrder);
            for (const auto& quadPoint : quadRule)
            {
                const Scalar qWeight = quadPoint.weight()*Extrusion::integrationElement(geometry, quadPoint.position());
                // Obtain and store shape function values and gradients at the current quad point
                IpData ipData(geometry, quadPoint.position(), localBasis);

                // get density from the problem
                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                    integralShapeFunctions[localDof.index()] += ipData.shapeValue(localDof.index())*qWeight;
            }

            for (const auto& localDof : nonCVLocalDofs(fvGeometry))
            {
                const auto localDofIdx = localDof.index();
                const auto& data = ipData(fvGeometry, localDof);
                const auto curDensity = problem.density(element, fvGeometry, data, false);
                const auto prevDensity = problem.density(element, fvGeometry, data, true);
                const auto curVelocity = curElemVolVars[localDofIdx].velocity();
                const auto prevVelocity = prevElemVolVars[localDofIdx].velocity();
                auto timeDeriv = (curDensity*curVelocity - prevDensity*prevVelocity);
                timeDeriv /= this->timeLoop().timeStepSize();

                // add storage to residual
                for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    residual[localDofIdx][eqIdx] += integralShapeFunctions[localDofIdx]*timeDeriv[eqIdx];
            }
        }
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVolumeVariables& elemVolVars,
                                           const ElementFluxVariablesCache& elemFluxVarsCache,
                                           const ElementBoundaryTypes &elemBcTypes) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {
            // Make sure we don't iterate over quadrature points if there are no hybrid dofs
            if( nonCVLocalDofs(fvGeometry).empty() )
                return;

            static const bool enableUnsymmetrizedVelocityGradient
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
            static const auto intOrder
                = getParamFromGroup<int>(problem.paramGroup(), "Assembly.FEIntegrationOrderFluxAndSource", 4);

            const auto &localBasis = fvGeometry.feLocalBasis();

            const auto& geometry = fvGeometry.elementGeometry();
            const auto& quadRule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), intOrder);
            for (const auto& quadPoint : quadRule)
            {
                const Scalar qWeight = quadPoint.weight()*Extrusion::integrationElement(geometry, quadPoint.position());

                // Obtain and store shape function values and gradients at the current quad point
                IpData ipData(geometry, quadPoint.position(), localBasis);
                FluxFunctionContext context(problem, fvGeometry, elemVolVars, ipData);
                const auto& v = context.velocity();
                const auto& gradV = context.gradVelocity();

                // get viscosity from the problem
                const Scalar mu = problem.effectiveViscosity(element, fvGeometry, ipData);
                // get density from the problem
                const Scalar density = problem.density(element, fvGeometry, ipData);

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto localDofIdx = localDof.index();
                    NumEqVector fluxAndSourceTerm(0.0);
                    // add advection term
                    if (problem.enableInertiaTerms())
                        fluxAndSourceTerm += density*(v*ipData.gradN(localDofIdx))*v;

                    // add diffusion term
                    fluxAndSourceTerm -= enableUnsymmetrizedVelocityGradient ?
                                            mu*mv(gradV, ipData.gradN(localDofIdx))
                                            : mu*mv(gradV + getTransposed(gradV), ipData.gradN(localDofIdx));

                    // add pressure term
                    fluxAndSourceTerm += problem.pressure(element, fvGeometry, ipData) * ipData.gradN(localDofIdx);

                    // add gravity term rho*g (note that gravity might be zero in case it's disabled in the problem)
                    fluxAndSourceTerm += density * problem.gravity();

                    // finally add source and Neumann term and add everything to residual
                    const auto sourceAtIp = problem.source(fvGeometry, elemVolVars, ipData);

                    for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    {
                        fluxAndSourceTerm[eqIdx] += ipData.shapeValue(localDofIdx) * sourceAtIp[eqIdx];
                        residual[localDofIdx][eqIdx] += qWeight*fluxAndSourceTerm[eqIdx];
                    }
                }
            }

            if(elemBcTypes.hasNeumann())
                residual += evalNeumannSegments_(problem, element, fvGeometry, elemVolVars, elemFluxVarsCache, elemBcTypes);
        }
    }

private:
    ElementResidualVector evalNeumannSegments_(const Problem& problem,
                                               const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFluxVariablesCache& elemFluxVarsCache,
                                               const ElementBoundaryTypes &elemBcTypes) const
    {
        ElementResidualVector flux(0.0);

        static const auto intOrder
            = getParamFromGroup<int>(problem.paramGroup(), "Assembly.FEIntegrationOrderBoundary", 4);

        const auto &localBasis = fvGeometry.feLocalBasis();

        const auto& geometry = fvGeometry.elementGeometry();
        for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            if(!intersection.boundary())
                continue;

            const auto& bcTypes = elemBcTypes.get(fvGeometry, intersection);
            if(!bcTypes.hasNeumann())
                continue;

            // select quadrature rule for intersection faces (dim-1)
            auto isGeometry = intersection.geometry();
            const auto& faceRule = Dune::QuadratureRules<Scalar, dim-1>::rule(isGeometry.type(), intOrder);
            for (const auto& quadPoint : faceRule)
            {
                // position of quadrature point in local coordinates of inside element
                auto local = geometry.local(isGeometry.global(quadPoint.position()));

                // get quadrature rule weight for intersection
                Scalar qWeight = quadPoint.weight() * Extrusion::integrationElement(isGeometry, quadPoint.position());
                FaceIpData faceIpData(geometry, local, localBasis, intersection.centerUnitOuterNormal(), BoundaryFlag{ intersection });

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto& boundaryFlux = qWeight*problem.boundaryFlux(fvGeometry, elemVolVars, elemFluxVarsCache, faceIpData);
                    // only add fluxes to equations for which Neumann is set
                    for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                        if (bcTypes.isNeumann(eqIdx))
                            flux[localDof.index()][eqIdx] += faceIpData.shapeValue(localDof.index()) * boundaryFlux[eqIdx];
                }

            }

        }

        return flux;
    }
};

} // end namespace Dumux

#endif
