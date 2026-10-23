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
#include <dumux/common/concepts/variables_.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/common/typetraits/griddiscretization.hh>

#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

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
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVariablesCache = Concept::GridVariablesCache_t<GridVariables>;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using Variables = Concept::Variables_t<GridVariables>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using SubControlVolume = typename ElementDiscretization::SubControlVolume;
    using SubControlVolumeFace = typename ElementDiscretization::SubControlVolumeFace;
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
    using FluxHelper = NavierStokesMomentumFluxCVFE<GridGeometry, NumEqVector>;
    using FeResidual = NavierStokesMomentumFELocalResidualTerms<Scalar, NumEqVector, LocalBasis, Extrusion>;

    using FluxTerms = Experimental::CVFENavierStokesMomentumFluxAtIp<Scalar, GridGeometry>;

public:
    //! Use the parent type's constructor
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage term of the equation
     *
     * \param problem The problem to solve
     * \param elemDisc The finite-volume geometry of the element
     * \param scv The sub-control volume over which we integrate the storage term
     * \param vars The variables associated with the scv
     * \param isPreviousStorage If set to true, the storage term is evaluated on the previous time level.
     *
     */
    NumEqVector computeStorage(const Problem& problem,
                               const ElementDiscretization& elemDisc,
                               const SubControlVolume& scv,
                               const Variables& vars,
                               const bool isPreviousStorage) const
    {
        return problem.density(elemDisc.element(), elemDisc, ipData(elemDisc, scv), isPreviousStorage) * vars.velocity();
    }

    /*!
     * \brief Calculate the storage integral
     *
     * \param elemDisc The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub control volume
     * \param isPreviousTimeLevel If set to true, the storage term is evaluated on the previous time level.
     *
     */
    NumEqVector storageIntegral(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        const auto& vars = elemVars[scv];
        // We apply mass lumping
        NumEqVector storage = this->asImp().problem().density(elemDisc.element(), elemDisc, ipData(elemDisc, scv), isPreviousTimeLevel)
                            * vars.velocity();

        storage *= Extrusion::volume(elemDisc, scv) * vars.extrusionFactor();

        return storage;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param elemDisc The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub-control volume over which we integrate the source term
     *
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const ElementDiscretization& elemDisc,
                              const ElementVariables& elemVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source;

        if constexpr (Detail::hasProblemSourceWithIpDataInterface<Problem, ElementDiscretization, ElementVariables, BaseIpData>())
        {
            source = problem.source(elemDisc, elemVars, ipData(elemDisc, scv.center()));

            // ToDo: point source data with ipData
            // add contribution from possible point sources
            if (!problem.pointSourceMap().empty())
                source += problem.scvPointSources(element, elemDisc, elemVars, scv);
        }
        else
            source = ParentType::computeSource(problem, element, elemDisc, elemVars, scv);


        // add rho*g (note that gravity might be zero in case it's disabled in the problem)
        const auto& data = ipData(elemDisc, scv);
        source +=  problem.density(element, elemDisc, data) * problem.gravity();

        // Axisymmetric problems in 2D feature an extra source term arising from the transformation to cylindrical coordinates.
        // See Ferziger/Peric: Computational methods for Fluid Dynamics (2020)
        // https://doi.org/10.1007/978-3-319-99693-6
        // Chapter 9.9 and Eq. (9.81) and comment on finite volume methods
        if constexpr (dim == 2 && isRotationalExtrusion<Extrusion>)
        {
            // the radius with respect to the rotation axis
            const auto& gridDiscretization = Dumux::gridDiscretization(elemDisc);
            const auto r = scv.center()[Extrusion::radialAxis] - gridDiscretization.bBoxMin()[Extrusion::radialAxis];

            // The velocity term is new with respect to Cartesian coordinates and handled below as a source term
            // It only enters the balance of the momentum balance in radial direction
            source[Extrusion::radialAxis] += -2.0*problem.effectiveViscosity(element, elemDisc, data)
                * elemVars[scv].velocity(Extrusion::radialAxis) / (r*r);

            // Pressure term (needed because we incorporate pressure in terms of a surface integral).
            // grad(p) becomes div(pI) + (p/r)*n_r in cylindrical coordinates. The second term
            // is new with respect to Cartesian coordinates and handled below as a source term.
            source[Extrusion::radialAxis] += problem.pressure(element, elemDisc, data)/r;
        }

        return source;
    }

    /*!
     * \brief Calculate the source integral
     *
     * \param elemDisc The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub control volume
     *
     */
    NumEqVector sourceIntegral(const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const SubControlVolume& scv) const
    {
        static_assert(!(dim == 2 && isRotationalExtrusion<Extrusion>), "Rotational extrusion source terms are not implemented for integral interface.");

        const auto& problem = this->asImp().problem();

        NumEqVector source(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv))
        {
            source += qpData.weight() * (problem.source(elemDisc, elemVars, qpData.ipData())
                                        + problem.density(elemDisc.element(), elemDisc, qpData.ipData()) * problem.gravity());
        }

        source *= elemVars[scv].extrusionFactor();

        // add contribution from possible point sources
        const auto& pointSources = problem.pointSources();
        if (!pointSources.empty())
            for (const auto& context : pointSources.contexts(elemDisc, scv))
            {
                auto psValues = pointSources.eval(elemDisc, elemVars, context);
                source += psValues;
            }

        return source;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a sub control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param elemDisc The finite volume geometry context
     * \param elemVars The variables for all local dofs of the element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    template<class ElementFluxVariablesCache>
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const ElementDiscretization& elemDisc,
                            const ElementVariables& elemVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        using FluxContext = NavierStokesMomentumFluxContext<Problem, ElementDiscretization, ElementVariables, ElementFluxVariablesCache>;
        FluxContext context(problem, elemDisc, elemVars, elemFluxVarsCache, scvf);
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
     * \param elemDisc The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scvf The sub control volume face
     *
     */
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        const auto& element = elemDisc.element();

        NumEqVector flux(0.0);
        GlobalPosition velocityIntegral(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
        {
            const auto& faceIpData = qpData.ipData();
            const auto& ipCache = cache(elemVars, faceIpData);
            if (problem.enableInertiaTerms())
                velocityIntegral.axpy(qpData.weight(), FluxTerms::velocity(elemDisc, elemVars, ipCache));

            const auto& normal = faceIpData.unitOuterNormal();
            FluxTerms::diffusiveFluxTerm(problem, element, elemDisc, elemVars, faceIpData, ipCache).usmv(qpData.weight(), normal, flux);
            FluxTerms::pressureFluxTerm(problem, element, elemDisc, faceIpData).usmv(qpData.weight(), normal, flux);
        }

        // the transported momentum is upwinded on the integrated volume flux over the face
        if (problem.enableInertiaTerms())
        {
            static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

            const auto& insideVars = elemVars[elemDisc.scv(scvf.insideScvIdx())];
            const auto& outsideVars = elemVars[elemDisc.scv(scvf.outsideScvIdx())];
            const bool insideIsUpwind = velocityIntegral*scvf.unitOuterNormal() > 0;
            const auto& upwindVelocity = insideIsUpwind ? insideVars.velocity() : outsideVars.velocity();
            const auto& downwindVelocity = insideIsUpwind ? outsideVars.velocity() : insideVars.velocity();

            // both are face quantities here, the density is evaluated at the face center
            const auto faceIpData = ipData(elemDisc, scvf);
            const Scalar density = problem.density(element, elemDisc, faceIpData);
            const auto momentum = [&](const auto&, const auto&, const auto&)
            { return density*(upwindWeight*upwindVelocity + (1.0 - upwindWeight)*downwindVelocity); };
            const auto transportingVelocity = [&](const auto&, const auto&, const auto&) { return velocityIntegral; };
            FluxTerms::advectiveFluxTerm(elemDisc, elemVars, faceIpData, momentum, transportingVelocity)
                .umv(scvf.unitOuterNormal(), flux);
        }

        flux *= elemVars[elemDisc.scv(scvf.insideScvIdx())].extrusionFactor();

        return flux;
    }

    void addToElementStorage(ElementResidualVector& storage,
                             const Problem& problem,
                             const Element& element,
                             const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             bool isPreviousTimeLevel) const
    {
        FeResidual::addStorageTerms(storage, problem, elemDisc, elemVars, isPreviousTimeLevel);
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        FeResidual::addFluxAndSourceTerms(
            residual, problem, elemDisc, elemVars
        );
    }

};

} // end namespace Dumux

#endif
