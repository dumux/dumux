// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_FLUXVARIABLES_HH

#include <array>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>

// #include "staggeredupwindhelper.hh"
#include "velocitygradients.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag>
class NavierStokesMomentumFluxVariables
{
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using VelocityGradients = StaggeredVelocityGradients;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Extrusion = Extrusion_t<GridGeometry>;

    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

public:

    NavierStokesMomentumFluxVariables(const Problem& problem,
                                      const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const SubControlVolumeFace& scvFace,
                                      const ElementVolumeVariables& elemVolVars,
                                      const ElementFluxVariablesCache& elemFluxVarsCache,
                                      const ElementBoundaryTypes& elemBcTypes)
    : problemPtr_(&problem)
    , elementPtr_(&element)
    , fvGeometryPtr_(&fvGeometry)
    , scvFacePtr_(&scvFace)
    , elemVolVarsPtr_(&elemVolVars)
    , elemFluxVarsCachePtr_(&elemFluxVarsCache)
    , elemBcTypesPtr_(&elemBcTypes)
    {
        static_assert(
            std::decay_t<decltype(problem.dirichlet(element, scvFace))>::size()
                == static_cast<std::size_t>(GridView::dimension),
            "Expects problem.dirichlet to return an array with as many entries as dimensions."
        );
    }

    const Problem& problem() const
    { return *problemPtr_; }

    const Element& element() const
    { return *elementPtr_; }

    const SubControlVolumeFace& scvFace() const
    { return *scvFacePtr_; }

    const FVElementGeometry& fvGeometry() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars() const
    { return *elemVolVarsPtr_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return *elemFluxVarsCachePtr_; }

    const ElementBoundaryTypes& elemBcTypes() const
    { return *elemBcTypesPtr_; }

    /*!
     * \brief Returns the diffusive momentum flux due to viscous forces
     */
    NumEqVector advectiveMomentumFlux() const
    {
        if (!this->problem().enableInertiaTerms())
            return NumEqVector(0.0);

        if (this->scvFace().isFrontal())
            return frontalAdvectiveMomentumFlux();
        else
            return lateralAdvectiveMomentumFlux();
    }

    /*!
     * \brief Returns the diffusive momentum flux due to viscous forces
     */
    NumEqVector diffusiveMomentumFlux() const
    {
        if (this->scvFace().isFrontal())
            return frontalDiffusiveMomentumFlux();
        else
            return lateralDiffusiveMomentumFlux();
    }

    /*!
     * \brief Returns the frontal part of the momentum flux.
     *        This treats the flux over the staggered face at the center of an element,
     *        parallel to the current scvf where the velocity dof of interest lives.
     *
     * \verbatim
     *                    scvf
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scv
     *              |       #      |                 # staggered face over which fluxes are calculated
     *      vel.Opp |~~>    O~~~>  x~~~~> vel.Self
     *              |       #      |                 x dof position
     *              |       #      |
     *              --------========                 -- element
     *                   scvf
     *                                               O integration point
     * \endverbatim
     */
    NumEqVector frontalDiffusiveMomentumFlux() const
    {
        const auto& scvf = this->scvFace();
        assert(scvf.isFrontal());

        NumEqVector result(0.0);
        const auto& fvGeometry = this->fvGeometry();
        const auto& elemVolVars = this->elemVolVars();

        if (scvf.boundary())
            return result;

        // get the velocity gradient at the normal face's integration point
        const auto velGradII = VelocityGradients::velocityGradII(fvGeometry, scvf, elemVolVars);

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(this->problem().paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
        static const Scalar factor = enableUnsymmetrizedVelocityGradient ? 1.0 : 2.0;

        const auto mu = this->problem().effectiveViscosity(this->element(), this->fvGeometry(), this->scvFace());
        result -= factor * mu * velGradII * Extrusion::area(fvGeometry, scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor() * scvf.directionSign();

        static const bool enableDilatationTerm = getParamFromGroup<bool>(this->problem().paramGroup(), "FreeFlow.EnableDilatationTerm", false);
        if (enableDilatationTerm)
        {
            Scalar divergence = velGradII;
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto otherFrontalScvf = *(scvfs(fvGeometry, scv).begin());
                assert(otherFrontalScvf.isFrontal() && !otherFrontalScvf.boundary());
                if (otherFrontalScvf.index() != scvf.index())
                    divergence += VelocityGradients::velocityGradII(fvGeometry, otherFrontalScvf, elemVolVars);
            }

            result += 2.0/3.0 * mu * divergence * scvf.directionSign() * Extrusion::area(fvGeometry, scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        }

        return result;
    }

    /*!
     * \brief Returns the diffusive momentum flux over the staggered face perpendicular to the scvf
     *        where the velocity dof of interest lives (coinciding with the element's scvfs).
     *
     * \verbatim
     *              ----------------
     *              |              |vel.
     *              |    in.norm.  |Parallel
     *              |       vel.   |~~~~>
     *              |       ^      |        ^ out.norm.vel.
     *              |       |      |        |
     *       scvf   ---------#######:::::::::       || and # staggered half-control-volume (own element)
     *              |      ||      | curr. ::
     *              |      ||      | scvf  ::       :: staggered half-control-volume (neighbor element)
     *              |      ||      x~~~~>  ::
     *              |      ||      | vel.  ::       # lateral staggered faces over which fluxes are calculated
     *        scvf  |      ||      | Self  ::
     *              ---------#######:::::::::       x dof position
     *                 scvf
     *                                              -- elements
     * \endverbatim
     */
    NumEqVector lateralDiffusiveMomentumFlux() const
    {
        const auto& scvf = this->scvFace();
        assert(scvf.isLateral());

        NumEqVector result(0.0);
        const auto& fvGeometry = this->fvGeometry();
        const auto& elemVolVars = this->elemVolVars();
        const auto& problem = this->problem();
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        const auto mu = this->problem().effectiveViscosity(this->element(), this->fvGeometry(), this->scvFace());

        // get the velocity gradient at the lateral face's integration point
        const auto gradV = VelocityGradients::velocityGradient(fvGeometry, scvf, elemVolVars, this->elemBcTypes(), false);

        // Consider the shear stress caused by the gradient of the velocities parallel to our face of interest.
        GlobalPosition gradVn(0.0);
        gradV.mv(scvf.unitOuterNormal(), gradVn);
        const Scalar velocityGrad_ij = gradVn[scv.dofAxis()];
        result -= mu * velocityGrad_ij;

        // Consider the shear stress caused by the gradient of the velocities normal to our face of interest.
        if (!enableUnsymmetrizedVelocityGradient)
        {
            GlobalPosition gradVTransposedN(0.0);
            gradV.mtv(scvf.unitOuterNormal(), gradVTransposedN);
            const Scalar velocityGrad_ji = gradVTransposedN[scv.dofAxis()];
            result -= mu * velocityGrad_ji;
        }

        // Account for the area of the staggered lateral face.
        return result * Extrusion::area(fvGeometry, scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
    }

    /*!
     * \brief Returns the frontal pressure contribution.
     *
     * \verbatim
     *
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scv
     *              |       #      |                 # frontal staggered face for which pressure contribution is calculated
     *              |       P      x
     *              |       #      |                 x dof position
     *              |       #      |
     *              --------========                 -- element
     *
     *                                               P integration point, pressure DOF lives here
     * \endverbatim
     */
    NumEqVector pressureContribution() const
    {
        NumEqVector result(0.0);
        const auto& scvf = this->scvFace();
        if (scvf.isLateral() || scvf.boundary())
            return result;

        // The pressure force needs to take the extruded scvf area into account.
        const auto pressure = this->problem().pressure(this->element(), this->fvGeometry(), scvf);
        result = pressure*Extrusion::area(this->fvGeometry(), scvf)*this->elemVolVars()[scvf.insideScvIdx()].extrusionFactor();

        // The pressure contribution calculated above might have a much larger numerical value compared to the viscous or inertial forces.
        // This may lead to numerical inaccuracies due to loss of significance (cancellantion) for the final residual value.
        // In the end, we are only interested in a pressure difference between the two relevant faces so we can
        // subtract a reference value from the actual pressure contribution. Assuming an axisparallel cartesian grid,
        // scvf.area() will have the same value at both opposing faces such that the reference pressure contribution
        // cancels out in the final residual which combines the pressure contribution of two adjacent elements
        // We explicitly do extrude the area here because that might yield different results in both elements.
        // The multiplication by scvf.area() aims at having a reference value of the same order of magnitude as the actual pressure contribution.
        const auto referencePressure = this->problem().referencePressure(this->element(), this->fvGeometry(), scvf);
        result -= referencePressure*scvf.area();

        // Account for the orientation of the face.
        result *= scvf.directionSign();
        return result;
    }

    /*!
     * \brief Returns the frontal part of the momentum flux.
     *        This treats the flux over the staggered face at the center of an element,
     *        parallel to the current scvf where the velocity dof of interest lives.
     *
     * \verbatim
     *                    scvf
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scv
     *              |       #      |                 # staggered face over which fluxes are calculated
     *      vel.Opp |~~>    O~~~>  x~~~~> vel.Self
     *              |       #      |                 x dof position
     *              |       #      |
     *              --------========                 -- element
     *                   scvf
     *                                               O integration point
     * \endverbatim
     */
    NumEqVector frontalAdvectiveMomentumFlux() const
    {
        NumEqVector flux(0.0);
        const auto& scvf = this->scvFace();
        assert(scvf.isFrontal());

        const auto& problem = this->problem();
        const auto& elemVolVars = this->elemVolVars();
        const auto velocitySelf = elemVolVars[scvf.insideScvIdx()].velocity();
        const auto velocityOpposite = elemVolVars[scvf.outsideScvIdx()].velocity();

        // Get the average velocity at the center of the element (i.e. the location of the staggered face).
        const Scalar transportingVelocity = (velocitySelf + velocityOpposite) * 0.5;
        const Scalar density = this->problem().density(this->element(), this->fvGeometry(), scvf);
        const bool selfIsUpstream = scvf.directionSign() == sign(transportingVelocity);

        // TODO use higher order helper
        static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");
        const Scalar transportedMomentum = selfIsUpstream ? (upwindWeight * velocitySelf + (1.0 - upwindWeight) * velocityOpposite) * density
                                                          : (upwindWeight * velocityOpposite + (1.0 - upwindWeight) * velocitySelf) * density;

        return  transportingVelocity * transportedMomentum * scvf.directionSign() * Extrusion::area(this->fvGeometry(), scvf) * extrusionFactor_(elemVolVars, scvf);
    }

    /*!
     * \brief Returns the advective momentum flux over the staggered face perpendicular to the scvf
     *        where the velocity dof of interest lives (coinciding with the element's scvfs).
     *
     * \verbatim
     *              ----------------
     *              |     inner    |
     *              |    transp.   |
     *              |      vel.    |~~~~> outer vel.
     *              |       ^      |
     *              |       |      |
     *              ---------######O                 || and # staggered half-control-volume
     *              |      ||      | scv
     *              |      ||      |                 # lateral staggered faces over which fluxes are calculated
     *              |      ||      x~~~~> inner vel.
     *              |      ||      |                 x dof position
     *              |      ||      |
     *              ---------#######                -- elements
     *
     *                                               O integration point
     * \endverbatim
     */
    NumEqVector lateralAdvectiveMomentumFlux() const
    {
        NumEqVector flux(0.0);
        const auto& scvf = this->scvFace();
        assert(scvf.isLateral());

        const auto& problem = this->problem();
        const auto& elemVolVars = this->elemVolVars();
        const auto fvGeometry = this->fvGeometry();

        // get the transporting velocity which is perpendicular to our own (inner) velocity
        const Scalar transportingVelocity = [&]()
        {
            const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
            const Scalar innerTransportingVelocity = elemVolVars[orthogonalScvf.insideScvIdx()].velocity();

            static const bool useOldScheme = getParam<bool>("FreeFlow.UseOldTransportingVelocity", true); // TODO how to deprecate?

            if (useOldScheme)
            {
                if (scvf.boundary() && fvGeometry.scv(scvf.insideScvIdx()).boundary())
                {
                    if (this->elemBcTypes()[scvf.localIndex()].isDirichlet(scvf.normalAxis()))
                        return problem.dirichlet(this->element(), scvf)[scvf.normalAxis()];
                }
                else
                    return innerTransportingVelocity;
            }

            // use the Dirichlet velocity as transporting velocity if the lateral face is on a Dirichlet boundary
            if (scvf.boundary())
            {
                if (this->elemBcTypes()[scvf.localIndex()].isDirichlet(scvf.normalAxis()))
                    return 0.5*(problem.dirichlet(this->element(), scvf)[scvf.normalAxis()] + innerTransportingVelocity);
            }

            if (orthogonalScvf.boundary())
            {
                if (this->elemBcTypes()[orthogonalScvf.localIndex()].isDirichlet(scvf.normalAxis()))
                    return 0.5*(problem.dirichlet(this->element(), scvf)[scvf.normalAxis()] + innerTransportingVelocity);
                else
                    return innerTransportingVelocity; // fallback value, should actually never be called
            }

            // average the transporting velocity by weighting with the scv volumes
            const auto insideVolume = fvGeometry.scv(orthogonalScvf.insideScvIdx()).volume();
            const auto outsideVolume = fvGeometry.scv(orthogonalScvf.outsideScvIdx()).volume();
            const auto outerTransportingVelocity = elemVolVars[orthogonalScvf.outsideScvIdx()].velocity();
            return (insideVolume*innerTransportingVelocity + outsideVolume*outerTransportingVelocity) / (insideVolume + outsideVolume);
        }();

        const Scalar transportedMomentum = [&]()
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

            auto getDirichletMomentumFlux = [&]()
            {
                return problem.dirichlet(this->element(), scvf)[insideScv.dofAxis()] * this->problem().density(this->element(), insideScv);
            };

            // use the Dirichlet velocity as for transported momentum if the lateral face is on a Dirichlet boundary
            if (scvf.boundary())
            {
                if (!this->elemBcTypes()[scvf.localIndex()].isDirichlet(insideScv.dofAxis()))
                    DUNE_THROW(Dune::InvalidStateException, "Neither Dirichlet nor Neumann BC set at " << scvf.ipGlobal());

                return getDirichletMomentumFlux();
            }
            else
            {
                if (fvGeometry.scvfIntegrationPointInConcaveCorner(scvf))
                {
                    // TODO: we could put this into an assert, as the construction of outsideScvfWithSameIntegrationPoint is quite expensive
                    const auto& outsideScvfWithSameIntegrationPoint = fvGeometry.outsideScvfWithSameIntegrationPoint(scvf);
                    if (!this->problem().boundaryTypes(this->element(), outsideScvfWithSameIntegrationPoint).isDirichlet(insideScv.dofAxis()))
                        DUNE_THROW(Dune::InvalidStateException, "Neither Dirichlet nor Neumann BC set at " << outsideScvfWithSameIntegrationPoint.ipGlobal());

                    return getDirichletMomentumFlux();
                }
            }

            const bool selfIsUpstream = scvf.directionSign() == sign(transportingVelocity);

            const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
            const auto outerVelocity = elemVolVars[scvf.outsideScvIdx()].velocity();

            const auto rho = this->problem().insideAndOutsideDensity(this->element(), fvGeometry, scvf);

            const auto insideMomentum = innerVelocity * rho.first;
            const auto outsideMomentum = outerVelocity * rho.second;

            // TODO use higher order helper
            static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

            return selfIsUpstream ? (upwindWeight * insideMomentum + (1.0 - upwindWeight) * outsideMomentum)
                                  : (upwindWeight * outsideMomentum + (1.0 - upwindWeight) * insideMomentum);
        }();

        return  transportingVelocity * transportedMomentum * scvf.directionSign() * Extrusion::area(fvGeometry, scvf) * extrusionFactor_(elemVolVars, scvf);
    }

private:

    template<class ElementVolumeVariables, class SubControlVolumeFace>
    Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf) const
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
    }


    const Problem* problemPtr_;                             //!< Pointer to the problem
    const Element* elementPtr_;                             //!< Pointer to the element at hand
    const FVElementGeometry* fvGeometryPtr_;                //!< Pointer to the current FVElementGeometry
    const SubControlVolumeFace* scvFacePtr_;                //!< Pointer to the sub control volume face for which the flux variables are created
    const ElementVolumeVariables* elemVolVarsPtr_;          //!< Pointer to the current element volume variables
    const ElementFluxVariablesCache* elemFluxVarsCachePtr_; //!< Pointer to the current element flux variables cache
    const ElementBoundaryTypes* elemBcTypesPtr_; //!< Pointer to element boundary types
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
