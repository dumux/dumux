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
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_FLUXVARIABLES_HH

#include <array>
#include <optional>

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

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

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
    {}

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
        const auto gradV = VelocityGradients::velocityGradient(fvGeometry, scvf, elemVolVars);

        GlobalPosition gradVn(0.0);
        gradV.mv(scvf.unitOuterNormal(), gradVn);
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const Scalar velocityGrad_ii = gradVn[scv.directionIndex()];

        static const bool enableUnsymmetrizedVelocityGradient
        = getParamFromGroup<bool>(this->problem().paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
        const Scalar factor = enableUnsymmetrizedVelocityGradient ? 1.0 : 2.0;

        const auto mu = this->problem().effectiveViscosity(this->element(), this->fvGeometry(), this->scvFace());
        result -= factor * mu * velocityGrad_ii * Extrusion::area(scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();

        static const bool enableDilatationTerm = getParamFromGroup<bool>(this->problem().paramGroup(), "FreeFlow.EnableDilatationTerm", false);
        if (enableDilatationTerm)
        {
            Scalar divergence = 0.0;
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto frontalScvf  = *(scvfs(fvGeometry, scv).begin());
                assert(frontalScvf.isFrontal() && !frontalScvf.boundary());
                divergence += VelocityGradients::velocityGradII(fvGeometry, frontalScvf, elemVolVars);
            }
            // std::cout << "divergence at " << scvf.center() << " is " << divergence << std::endl;
            // std::cout << std::setprecision(15) << "old term " << factor * mu * velocityGrad_ii * scvf.area() * elemVolVars[scvf.insideScvIdx()].extrusionFactor() << ", div term " << 2.0/3.0 * mu * divergence * scvf.directionSign() * scvf.area() * elemVolVars[scvf.insideScvIdx()].extrusionFactor() << std::endl;
            result += 2.0/3.0 * mu * divergence * scvf.directionSign() * Extrusion::area(scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
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
        const auto gradV = VelocityGradients::velocityGradient(fvGeometry, scvf, elemVolVars);

        // Consider the shear stress caused by the gradient of the velocities parallel to our face of interest.
        GlobalPosition gradVn(0.0);
        gradV.mv(scvf.unitOuterNormal(), gradVn);
        const Scalar velocityGrad_ij = gradVn[scv.directionIndex()];
        result -= mu * velocityGrad_ij;

        // Consider the shear stress caused by the gradient of the velocities normal to our face of interest.
        if (!enableUnsymmetrizedVelocityGradient)
        {
            GlobalPosition gradVTransposedN(0.0);
            gradV.mtv(scvf.unitOuterNormal(), gradVTransposedN);
            const Scalar velocityGrad_ji = gradVTransposedN[scv.directionIndex()];
            result -= mu * velocityGrad_ji;
        }

        // Account for the area of the staggered lateral face.
        return result * Extrusion::area(scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
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
        result = pressure*Extrusion::area(scvf)*this->elemVolVars()[scvf.insideScvIdx()].extrusionFactor();

        // The pressure contribution calculated above might have a much larger numerical value compared to the viscous or inertial forces.
        // This may lead to numerical inaccuracies due to loss of significance (cancellantion) for the final residual value.
        // In the end, we are only interested in a pressure difference between the two relevant faces so we can
        // substract a reference value from the actual pressure contribution. Assuming an axisparallel cartesian grid,
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

        return  transportingVelocity * transportedMomentum * scvf.directionSign() * Extrusion::area(scvf) * extrusionFactor_(elemVolVars, scvf);
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
            static const bool useOldScheme = getParam<bool>("FreeFlow.UseOldTransportingVelocity", true); // TODO how to deprecate?

            // use the Dirichlet velocity as transporting velocity if the lateral face is on a Dirichlet boundary
            if (!useOldScheme && scvf.boundary())
            {
                if (this->elemBcTypes()[scvf.localIndex()].isDirichlet(scvf.directionIndex()))
                    return problem.dirichlet(this->element(), scvf)[scvf.directionIndex()];
            }

            const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
            const Scalar innerTransportingVelocity = elemVolVars[orthogonalScvf.insideScvIdx()].velocity();

            if (useOldScheme)
                return innerTransportingVelocity;

            if (orthogonalScvf.boundary())
            {
                if (this->elemBcTypes()[orthogonalScvf.localIndex()].isDirichlet(scvf.directionIndex()))
                    return problem.dirichlet(this->element(), scvf)[scvf.directionIndex()];
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
            // use the Dirichlet velocity as for transported momentum if the lateral face is on a Dirichlet boundary
            if (scvf.boundary())
            {
                if (const auto& scv = fvGeometry.scv(scvf.insideScvIdx()); this->elemBcTypes()[scvf.localIndex()].isDirichlet(scv.directionIndex()))
                    return problem.dirichlet(this->element(), scvf)[scv.directionIndex()] * this->problem().density(this->element(), scv);
            }

            const bool selfIsUpstream = scvf.directionSign() == sign(transportingVelocity);

            const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
            const auto outerVelocity = elemVolVars[scvf.outsideScvIdx()].velocity();

            const auto rho = this->problem().getInsideAndOutsideDensity(this->element(), fvGeometry, scvf);

            const auto insideMomentum = innerVelocity * rho.first;
            const auto outsideMomentum = outerVelocity * rho.second;

            // TODO use higher order helper
            static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

            return selfIsUpstream ? (upwindWeight * insideMomentum + (1.0 - upwindWeight) * outsideMomentum)
                                  : (upwindWeight * outsideMomentum + (1.0 - upwindWeight) * insideMomentum);
        }();

        return  transportingVelocity * transportedMomentum * scvf.directionSign() * Extrusion::area(scvf) * extrusionFactor_(elemVolVars, scvf);
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






//     /*!
//      * \brief Returns the momentum flux over the staggered faces
//      *        perpendicular to the scvf where the velocity dof of interest
//      *         lives (coinciding with the element's scvfs).
//      *
//      * \verbatim
//      *                scvf
//      *              ---------#######                 || and # staggered half-control-volume
//      *              |      ||      | current scvf
//      *              |      ||      |                 # normal staggered sub faces over which fluxes are calculated
//      *              |      ||      x~~~~> vel.Self
//      *              |      ||      |                 x dof position
//      *        scvf  |      ||      |
//      *              --------########                -- element
//      *                 scvf
//      * \endverbatim
//      */
//     FacePrimaryVariables computeLateralMomentumFlux(const Problem& problem,
//                                                     const Element& element,
//                                                     const SubControlVolumeFace& scvf,
//                                                     const FVElementGeometry& fvGeometry,
//                                                     const FaceFVElementGeometry& faceFVGeometry,
//                                                     const StaggeredScv& staggeredScv,
//                                                     const ElementVolumeVariables& elemVolVars,
//                                                     const ElementFaceVariables& elemFaceVars,
//                                                     const GridFluxVariablesCache& gridFluxVarsCache)
//     {
//         FacePrimaryVariables lateralFlux(0.0);
//         const auto& faceVars = elemFaceVars[staggeredScv];

//         // If the current scvf is on a boundary, check if there is a Neumann BC for the stress in tangential direction.
//         // Create a boundaryTypes object (will be empty if not at a boundary).
//         std::optional<BoundaryTypes> currentScvfBoundaryTypes;
//         if (staggeredScv.boundary())
//         {
//             // TODO
//             for (const auto& staggeredScvf : scvfs(faceFVGeometry))
//             {
//                 if (staggeredScvf.boundary() && staggeredScvf.isFrontal())
//                     currentScvfBoundaryTypes.emplace(problem.boundaryTypes(element, staggeredScvf));
//             }
//         }

//         // Account for all sub faces.
//         int localSubFaceIdx = 0;
//         for (const auto& staggeredScvf : scvfs(faceFVGeometry))
//         {
//             if (!staggeredScvf.isLateral() || staggeredScvf.insideScvIdx() != staggeredScv.dofIndex())
//                 continue;

//             // Create a boundaryTypes object (will be empty if not at a boundary).
//             std::optional<BoundaryTypes> lateralFaceBoundaryTypes;

//             // Check if there is face/element parallel to our face of interest where the dof lives on. If there is no parallel neighbor,
//             // we are on a boundary where we have to check for boundary conditions.
//             if (staggeredScvf.boundary())
//             {
//                 // Retrieve the boundary types that correspond to the center of the lateral scvf. As a convention, we always query
//                 // the type of BCs at the center of the element's "actual" lateral scvf (not the face of the staggered control volume).
//                 // The value of the BC will be evaluated at the center of the staggered face.
//                 //     --------###T##V                 || frontal face of staggered half-control-volume
//                 //     |      ||      | current scvf    #  lateral staggered face of interest (may lie on a boundary)
//                 //     |      ||      |                 x  dof position
//                 //     |      ||      x~~~~> vel.Self   -- element boundaries
//                 //     |      ||      |                 T  position at which the type of boundary conditions will be evaluated
//                 //     |      ||      |                    (center of lateral scvf)
//                 //     ----------------                 V  position at which the value of the boundary conditions will be evaluated
//                 //                                         (center of the staggered lateral face)
//                 lateralFaceBoundaryTypes.emplace(problem.boundaryTypes(element, staggeredScvf));
//             }

//             // If the current scvf is on a bounary and if there is a Neumann or Beavers-Joseph-(Saffmann) BC for the stress in tangential direction,
//             // assign this value for the lateral momentum flux. No further calculations are required. We assume that all lateral faces
//             // have the same type of BC (Neumann or Beavers-Joseph-(Saffmann)), but we sample the value at their actual positions.
//             if (currentScvfBoundaryTypes)
//             {
//                 // Handle Neumann BCs.
//                 if (currentScvfBoundaryTypes->isNeumann(Indices::velocity(staggeredScvf.directionIndex())))
//                 {
//                     // TODO pass elemFluxVarsCache
//                     lateralFlux += problem.neumann(element, faceFVGeometry, elemFaceVars, staggeredScvf)[Indices::velocity(staggeredScvf.directionIndex())] * staggeredScvf.area() * staggeredScvf.directionSign();
//                     ++localSubFaceIdx;
//                     continue;
//                 }
//             }

//             // Check if the lateral face (perpendicular to our current scvf) lies on a boundary. If yes, boundary conditions might need to be treated
//             // and further calculations can be skipped.
//             if (staggeredScvf.boundary())
//             {
//                 // Check if we have a symmetry boundary condition. If yes, the tangential part of the momentum flux can be neglected
//                 // and we may skip any further calculations for the given sub face.
//                 if (lateralFaceBoundaryTypes->isSymmetry())
//                 {
//                     ++localSubFaceIdx;
//                     continue;
//                 }

//                 // Handle Neumann boundary conditions. No further calculations are then required for the given sub face.
//                 if (lateralFaceBoundaryTypes->isNeumann(Indices::velocity(staggeredScv.directionIndex())))
//                 {
//                     lateralFlux +=  problem.neumann(element, faceFVGeometry, elemFaceVars, staggeredScvf)[Indices::velocity(staggeredScv.directionIndex())] * staggeredScvf.area() * staggeredScvf.directionSign();
//                     ++localSubFaceIdx;
//                     continue;
//                 }

//                 // Handle wall-function fluxes (only required for RANS models)
//                 if (incorporateWallFunction_(lateralFlux, problem, element, fvGeometry, scvf, elemVolVars, elemFaceVars, localSubFaceIdx))
//                 {
//                     ++localSubFaceIdx;
//                     continue;
//                 }
//             }

//             // Check the consistency of the boundary conditions, exactly one of the following must be set
//             if (lateralFaceBoundaryTypes)
//             {
//                 std::bitset<3> admittableBcTypes;
//                 admittableBcTypes.set(0, lateralFaceBoundaryTypes->isDirichlet(Indices::pressureIdx));
//                 admittableBcTypes.set(1, lateralFaceBoundaryTypes->isDirichlet(Indices::velocity(staggeredScv.directionIndex())));
//                 admittableBcTypes.set(2, lateralFaceBoundaryTypes->isBeaversJoseph(Indices::velocity(staggeredScv.directionIndex())));
//                 if (admittableBcTypes.count() != 1)
//                 {
//                     DUNE_THROW(Dune::InvalidStateException, "Invalid boundary conditions for lateral scvf "
//                     "for the momentum equations at global position " << lateralStaggeredFaceCenter_(scvf, localSubFaceIdx)
//                     << ", current scvf global position " << scvf.center());
//                 }
//             }

//             // If none of the above boundary conditions apply for the given sub face, proceed to calculate the tangential momentum flux.
//             // TODO advective fluxes
//             // if (problem.enableInertiaTerms())
//             //     lateralFlux += computeAdvectivePartOfLateralMomentumFlux_(problem, fvGeometry, element,
//             //                                                               scvf, faceFVGeometry, staggeredScvf,
//             //                                                               elemVolVars, faceVars,
//             //                                                               gridFluxVarsCache,
//             //                                                               currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
//             //                                                               localSubFaceIdx);

//             lateralFlux += computeDiffusivePartOfLateralMomentumFlux_(problem, fvGeometry, element,
//                                                                       scvf, faceFVGeometry, staggeredScvf,
//                                                                       elemVolVars, faceVars,
//                                                                       currentScvfBoundaryTypes, lateralFaceBoundaryTypes,
//                                                                       localSubFaceIdx);

//             ++localSubFaceIdx;
//         }
//         return lateralFlux;
//     }

//     /*!
//      * \brief Returns the momentum flux over an inflow or outflow boundary face.
//      *
//      * \verbatim
//      *                    scvf      //
//      *              ---------=======//               == and # staggered half-control-volume
//      *              |      ||      #// current scvf
//      *              |      ||      #//               # staggered boundary face over which fluxes are calculated
//      *              |      ||      x~~~~> vel.Self
//      *              |      ||      #//               x dof position
//      *        scvf  |      ||      #//
//      *              --------========//               -- element
//      *                   scvf       //
//      *                                              // boundary
//      * \endverbatim
//      */
//     FacePrimaryVariables inflowOutflowBoundaryFlux(const Problem& problem,
//                                                    const Element& element,
//                                                    const FaceFVElementGeometry& faceFVGeometry,
//                                                    const StaggeredScvf& scvf,
//                                                    const ElementFaceVariables& elemFaceVars) const
//     {
//         FacePrimaryVariables inOrOutflow(0.0);
//         const auto& scv = faceFVGeometry.scv(scvf.insideScvIdx());
//         const auto& context = problem.couplingManager().faceCouplingContext(element, scv);

//         const auto eIdx = problem.gridGeometry().elementMapper().index(element);
//         const auto& insideVolVars = context.elemVolVars[eIdx];

//         // Advective momentum flux.
//         if (problem.enableInertiaTerms())
//         {
//             // TODO
//             // const Scalar velocitySelf = elemFaceVars[scv].velocitySelf();
//             // const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
//             // const auto& upVolVars = (scvf.directionSign() == sign(velocitySelf)) ?
//             //                         insideVolVars : outsideVolVars;
//             //
//             // inOrOutflow += velocitySelf * velocitySelf * upVolVars.density();
//         }

//         // Apply a pressure at the boundary.
//         const Scalar boundaryPressure = normalizePressure
//                                         ? (problem.dirichlet(element, scvf)[Indices::pressureIdx] -
//                                            problem.initial(scvf)[Indices::pressureIdx])
//                                         : problem.dirichlet(element, scvf)[Indices::pressureIdx];
//         inOrOutflow += boundaryPressure;

//         // Account for the orientation of the face at the boundary,
//         return inOrOutflow * scvf.directionSign() * scvf.area() * insideVolVars.extrusionFactor();
//     }

// private:

//     /*!
//      * \brief Returns the advective momentum flux over the staggered face perpendicular to the scvf
//      *        where the velocity dof of interest lives (coinciding with the element's scvfs).
//      *
//      * \verbatim
//      *              ----------------
//      *              |              |
//      *              |    transp.   |
//      *              |      vel.    |~~~~> vel.Parallel
//      *              |       ^      |
//      *              |       |      |
//      *       scvf   ---------#######                 || and # staggered half-control-volume
//      *              |      ||      | current scvf
//      *              |      ||      |                 # normal staggered faces over which fluxes are calculated
//      *              |      ||      x~~~~> vel.Self
//      *              |      ||      |                 x dof position
//      *        scvf  |      ||      |
//      *              ---------#######                -- elements
//      *                 scvf
//      * \endverbatim
//      */
//     template<class StaggeredFVElementGeometry>
//     FacePrimaryVariables computeAdvectivePartOfLateralMomentumFlux_(const Problem& problem,
//                                                                     const FVElementGeometry& fvGeometry,
//                                                                     const Element& element,
//                                                                     const SubControlVolumeFace& scvf,
//                                                                     const StaggeredFVElementGeometry staggeredFVGeometry,
//                                                                     const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredScvf,
//                                                                     const ElementVolumeVariables& elemVolVars,
//                                                                     const ElementFaceVariables& elemFaceVars,
//                                                                     const GridFluxVariablesCache& gridFluxVarsCache,
//                                                                     const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
//                                                                     const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
//                                                                     const int localSubFaceIdx)
//     {
//         const auto eIdx = scvf.insideScvIdx();
//         const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

//         // Get the transporting velocity, located at the scvf perpendicular to the current scvf where the dof
//         // of interest is located.
//         const Scalar transportingVelocity = [&]()
//         {
//             const auto& faceVars = elemFaceVars[scvf];
//             if (!scvf.boundary())
//                 return faceVars.velocityLateralInside(localSubFaceIdx);
//             else
//             {
//                 // Create a boundaryTypes object. Get the boundary conditions. We sample the type of BC at the center of the current scvf.
//                 const BoundaryTypes bcTypes = problem.boundaryTypes(element, scvf);

//                 if (bcTypes.isDirichlet(Indices::velocity(lateralFace.directionIndex())))
//                 {
//                     // Construct a temporary scvf which corresponds to the staggered sub face, featuring the location
//                     // the staggered faces's center.
//                     const auto& lateralBoundaryFacePos = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
//                     return problem.dirichlet(element, scvf.makeBoundaryFace(lateralBoundaryFacePos))[Indices::velocity(lateralFace.directionIndex())];
//                 }
//                 else if (bcTypes.isBeaversJoseph(Indices::velocity(lateralFace.directionIndex())))
//                 {
//                     return VelocityGradients::beaversJosephVelocityAtCurrentScvf(problem, element, staggeredFVGeometry, staggeredScvf, faceVars,
//                                                                                  currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
//                 }
//                 else
//                     return faceVars.velocityLateralInside(localSubFaceIdx);
//             }
//         }();

//         const bool selfIsUpstream = lateralFace.directionSign() == sign(transportingVelocity);
//         StaggeredUpwindHelper<TypeTag, upwindSchemeOrder> upwindHelper(element, fvGeometry, scvf, elemFaceVars, elemVolVars, gridFluxVarsCache.staggeredUpwindMethods());
//         return upwindHelper.computeUpwindLateralMomentum(selfIsUpstream, lateralFace, localSubFaceIdx, currentScvfBoundaryTypes, lateralFaceBoundaryTypes)
//                * transportingVelocity * lateralFace.directionSign() * lateralFace.area() * 0.5 * extrusionFactor_(elemVolVars, lateralFace);
//     }

//     /*!
//      * \brief Returns the diffusive momentum flux over the staggered face perpendicular to the scvf
//      *        where the velocity dof of interest lives (coinciding with the element's scvfs).
//      *
//      * \verbatim
//      *              ----------------
//      *              |              |vel.
//      *              |    in.norm.  |Parallel
//      *              |       vel.   |~~~~>
//      *              |       ^      |        ^ out.norm.vel.
//      *              |       |      |        |
//      *       scvf   ---------#######:::::::::       || and # staggered half-control-volume (own element)
//      *              |      ||      | curr. ::
//      *              |      ||      | scvf  ::       :: staggered half-control-volume (neighbor element)
//      *              |      ||      x~~~~>  ::
//      *              |      ||      | vel.  ::       # lateral staggered faces over which fluxes are calculated
//      *        scvf  |      ||      | Self  ::
//      *              ---------#######:::::::::       x dof position
//      *                 scvf
//      *                                              -- elements
//      * \endverbatim
//      */
//     template<class StaggeredFVElementGeometry>
//     FacePrimaryVariables computeDiffusivePartOfLateralMomentumFlux_(const Problem& problem,
//                                                                     const FVElementGeometry& fvGeometry,
//                                                                     const Element& element,
//                                                                     const SubControlVolumeFace& scvf,
//                                                                     const StaggeredFVElementGeometry staggeredFVGeometry,
//                                                                     const typename StaggeredFVElementGeometry::StaggeredSubControlVolumeFace& staggeredScvf,
//                                                                     const ElementVolumeVariables& elemVolVars,
//                                                                     const FaceVariables& faceVars,
//                                                                     const std::optional<BoundaryTypes>& currentScvfBoundaryTypes,
//                                                                     const std::optional<BoundaryTypes>& lateralFaceBoundaryTypes,
//                                                                     const int localSubFaceIdx)
//     {
//         assert(staggeredScvf.isLateral());

//         const auto eIdx = scvf.insideScvIdx();
//         const auto& lateralFace = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

//         FacePrimaryVariables lateralDiffusiveFlux(0.0);

//         static const bool enableUnsymmetrizedVelocityGradient
//             = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

//         // Get the volume variables of the own and the neighboring element. The neighboring
//         // element is adjacent to the staggered face normal to the current scvf
//         // where the dof of interest is located.
//         // TODO use coupling context to get interpolated mu
//         const auto& insideVolVars = elemVolVars[lateralFace.insideScvIdx()];
//         const auto& outsideVolVars = elemVolVars[lateralFace.outsideScvIdx()];

//         // Get the averaged viscosity at the staggered face normal to the current scvf.
//         const Scalar muAvg = lateralFace.boundary()
//                              ? insideVolVars.effectiveViscosity()
//                              : (insideVolVars.effectiveViscosity() + outsideVolVars.effectiveViscosity()) * 0.5;

//         // Consider the shear stress caused by the gradient of the velocities normal to our face of interest.
//         if (!enableUnsymmetrizedVelocityGradient)
//         {
//             if (!scvf.boundary() ||
//                 currentScvfBoundaryTypes->isDirichlet(Indices::velocity(staggeredScvf.directionIndex())) ||
//                 currentScvfBoundaryTypes->isBeaversJoseph(Indices::velocity(staggeredScvf.directionIndex())))
//             {
//                 const Scalar velocityGrad_ji = VelocityGradients::velocityGradJI(problem, element, staggeredFVGeometry, staggeredScvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
//                 // Account for the orientation of the staggered normal face's outer normal vector.
//                 lateralDiffusiveFlux -= muAvg * velocityGrad_ji * staggeredScvf.directionSign();
//             }
//         }

//         // Consider the shear stress caused by the gradient of the velocities parallel to our face of interest.
//         const Scalar velocityGrad_ij = VelocityGradients::velocityGradIJ(problem, element, staggeredFVGeometry, staggeredScvf, faceVars, currentScvfBoundaryTypes, lateralFaceBoundaryTypes);
//         lateralDiffusiveFlux -= muAvg * velocityGrad_ij * staggeredScvf.directionSign();

//         // Account for the area of the staggered lateral face (0.5 of the coinciding scfv).
//         return lateralDiffusiveFlux * lateralFace.area() * 0.5 * extrusionFactor_(elemVolVars, lateralFace);
//     }

//     /*!
//      * \brief Get the location of the lateral staggered face's center.
//      *        Only needed for boundary conditions if the current scvf or the lateral one is on a bounary.
//      *
//      * \verbatim
//      *      --------#######o                 || frontal face of staggered half-control-volume
//      *      |      ||      | current scvf    #  lateral staggered face of interest (may lie on a boundary)
//      *      |      ||      |                 x  dof position
//      *      |      ||      x~~~~> vel.Self   -- element boundaries, current scvf may lie on a boundary
//      *      |      ||      |                 o  position at which the boundary conditions will be evaluated
//      *      |      ||      |                    (lateralStaggeredFaceCenter)
//      *      ----------------
//      * \endverbatim
//      */
//     const GlobalPosition& lateralStaggeredFaceCenter_(const SubControlVolumeFace& scvf, const int localSubFaceIdx) const
//     {
//         return scvf.pairData(localSubFaceIdx).lateralStaggeredFaceCenter;
//     };

//     //! helper function to get the averaged extrusion factor for a face
//     static Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf)
//     {
//         const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
//         const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
//         return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
//     }

//     //! do nothing if no turbulence model is used
//     template<class ...Args, bool turbulenceModel = ModelTraits::usesTurbulenceModel(), std::enable_if_t<!turbulenceModel, int> = 0>
//     bool incorporateWallFunction_(Args&&... args) const
//     { return false; }

//     //! if a turbulence model is used, ask the problem is a wall function shall be employed and get the flux accordingly
//     template<bool turbulenceModel = ModelTraits::usesTurbulenceModel(), std::enable_if_t<turbulenceModel, int> = 0>
//     bool incorporateWallFunction_(FacePrimaryVariables& lateralFlux,
//                                   const Problem& problem,
//                                   const Element& element,
//                                   const FVElementGeometry& fvGeometry,
//                                   const SubControlVolumeFace& scvf,
//                                   const ElementVolumeVariables& elemVolVars,
//                                   const ElementFaceVariables& elemFaceVars,
//                                   const std::size_t localSubFaceIdx) const
//     {
//         const auto eIdx = scvf.insideScvIdx();
//         const auto& lateralScvf = fvGeometry.scvf(eIdx, scvf.pairData(localSubFaceIdx).localLateralFaceIdx);

//         if (problem.useWallFunction(element, lateralScvf, Indices::velocity(scvf.directionIndex())))
//         {
//             const auto& lateralStaggeredFaceCenter = lateralStaggeredFaceCenter_(scvf, localSubFaceIdx);
//             const auto lateralBoundaryFace = lateralScvf.makeBoundaryFace(lateralStaggeredFaceCenter);
//             lateralFlux += problem.wallFunction(element, fvGeometry, elemVolVars, elemFaceVars, scvf, lateralBoundaryFace)[Indices::velocity(scvf.directionIndex())]
//                                                * extrusionFactor_(elemVolVars, lateralScvf) * lateralScvf.area() * 0.5;
//             return true;
//         }
//         else
//             return false;
//     }
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
