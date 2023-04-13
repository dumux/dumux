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
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_LOCAL_RESIDUAL_HH

#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/fclocalresidual.hh>
#include <dune/common/hybridutilities.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using the staggered discretization
 */
template<class TypeTag>
class NavierStokesMomentumResidual
: public FaceCenteredLocalResidual<TypeTag>
{
    using ParentType = FaceCenteredLocalResidual<TypeTag>;
    friend class FaceCenteredLocalResidual<TypeTag>;

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
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr auto dim = GridView::dimension;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param scv The sub-control volume over which we integrate the storage term
     * \param volVars The volume variables associated with the scv
     * \param isPreviousStorage Bool transferring the information if the storage term is computed at the current or previous time step
     * \note has to be implemented by the model specific residual class
     *
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars,
                               const bool isPreviousStorage = false) const
    {
        const auto& element = problem.gridGeometry().element(scv.elementIndex());
        return problem.density(element, scv, isPreviousStorage) * volVars.velocity();
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
     * \note This is the default implementation for all models as sources are computed
     *       in the user interface of the problem
     *
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source = ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);
        source += problem.gravity()[scv.dofAxis()] * problem.density(element, scv);

        // Axisymmetric problems in 2D feature an extra source terms arising from the transformation to cylindrical coordinates.
        // See Ferziger/Peric: Computational methods for fluid dynamics chapter 8.
        // https://doi.org/10.1007/978-3-540-68228-8 (page 301)
        if constexpr (dim == 2 && isRotationalExtrusion<Extrusion>)
        {
            if (scv.dofAxis() == Extrusion::radialAxis)
            {
                const auto r = scv.center()[scv.dofAxis()] - fvGeometry.gridGeometry().bBoxMin()[scv.dofAxis()];
                const auto& scvf = (*scvfs(fvGeometry, scv).begin()); // the frontal scvf belonging to the scv

                // Velocity term
                source -= -2.0*problem.effectiveViscosity(element, fvGeometry, scvf) * elemVolVars[scv].velocity() / (r*r);

                // Pressure term (needed because we incorporate pressure in terms of a surface integral).
                 // grad(p) becomes div(pI) + (p/r)*n_r in cylindrical coordinates. The second term
                 // is new with respect to Cartesian coordinates and handled below as a source term.
                source += problem.pressure(element, fvGeometry, scvf)/r;
            }
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
     * \param elemBcTypes The element boundary condition types
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache,
                            const ElementBoundaryTypes& elemBcTypes) const
    {
        FluxVariables fluxVars(problem, element, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, elemBcTypes);

        NumEqVector flux(0.0);
        flux += fluxVars.advectiveMomentumFlux();
        flux += fluxVars.diffusiveMomentumFlux();
        flux += fluxVars.pressureContribution();
        return flux;
    }

    //! Evaluate flux residuals for one sub control volume face when the element features at least one Dirichlet boundary condition
    //! Skip the flux calculation if the DOF of the associated sub control volume features an appropriate Dirichlet condition.
    NumEqVector maybeHandleDirichletBoundary(const Problem& problem,
                                             const Element& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const ElementBoundaryTypes& elemBcTypes,
                                             const ElementFluxVariablesCache& elemFluxVarsCache,
                                             const SubControlVolumeFace& scvf) const
    {
        if (const auto& scv = fvGeometry.scv(scvf.insideScvIdx()); scv.boundary())
        {
            const auto& frontalScvfOnBoundary = fvGeometry.frontalScvfOnBoundary(scv);
            const auto& bcTypes = elemBcTypes[frontalScvfOnBoundary.localIndex()];
            if (bcTypes.isDirichlet(scv.dofAxis()))
                return NumEqVector(0.0); // skip calculation as Dirichlet BC will be incorporated explicitly by assembler.
        }

        // Default: The scv does not lie on a boundary with a Dirichlet condition for the velocity in normal direction.
        // Check for a Neumann condition or  calculate the flux as normal.
        if (elemBcTypes.hasNeumann())
            return this->asImp().maybeHandleNeumannBoundary(problem, element, fvGeometry, elemVolVars, elemBcTypes, elemFluxVarsCache, scvf);
        else
            return this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache, elemBcTypes);
    }

    //! Evaluate flux residuals for one sub control volume face when the element features at least one Neumann boundary condition
    //! This requires special care.
    NumEqVector maybeHandleNeumannBoundary(const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVolumeVariables& elemVolVars,
                                           const ElementBoundaryTypes& elemBcTypes,
                                           const ElementFluxVariablesCache& elemFluxVarsCache,
                                           const SubControlVolumeFace& scvf) const
    {
        static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::fcstaggered); // TODO overload this method for different discretizations
        static_assert(
            std::decay_t<decltype(
                problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf)
            )>::size() == ModelTraits::dim(),
            "The momentum model expects problem.neumann to return a vector of size dim. "
            "When in doubt you should be able to use 'using NumEqVector = typename ParentType::NumEqVector;'."
        );
        assert(elemBcTypes.hasNeumann());

        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        if (scv.boundary())
        {
            // This is the most simple case for a frontal scvf lying directly on the boundary.
            // Just evaluate the Neumann flux for the normal velocity component.
            //
            //     ---------------- *                || frontal face of staggered half-control-volume
            //     |      ||      # *
            //     |      ||      # *                D  dof position (coincides with integration point where
            //     |      || scv  D *                   the Neumann flux is evaluated, here for x component)
            //     |      ||      # *
            //     |      ||      # *                -- element boundaries
            //     ---------------- *
            //  y                                     * domain boundary
            //  ^
            //  |                                     # current scvf over which Neumann flux is evaluated
            //  ----> x
            //
            if (scvf.boundary() && scvf.isFrontal())
            {
                const auto& bcTypes = elemBcTypes[scvf.localIndex()];
                if (bcTypes.hasNeumann() && bcTypes.isNeumann(scv.dofAxis()))
                {
                    const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                    return neumannFluxes[scv.dofAxis()] * Extrusion::area(fvGeometry, scvf) * elemVolVars[scv].extrusionFactor();
                }
            }
            else if (scvf.isLateral())
            {
                // If a sub control volume lies on a boundary, Neumann fluxes also need to be considered for its lateral faces,
                // even though those do not necessarily lie on a boundary (only their edges / integration points are touching the boundary).
                // We cannot simply calculate momentum fluxes across the lateral boundary in this situation because we might not
                // have enough information at the boundary to do so.
                // We first need to find the scv's frontal face on the boundary to check if there actually is a Neumann BC.
                // Afterwards, we use the actual (lateral) scvf to retrieve the Neumann flux at its integration point intersecting with the boundary.
                //
                //
                //     ---------######O *                || frontal face of staggered half-control-volume
                //     |      ||      : *
                //     |      ||      : *                D  dof position
                //     |      || scv  D *
                //     |      ||      : *               -- element boundaries
                //     |      ||      : *
                //     ---------------- *                * domain boundary
                //
                //  y                                    # current scvf over which Neumann flux is evaluated
                //  ^
                //  |                                    : frontal scvf on boundary
                //  ----> x
                //                                       O integration point at which Neumann flux is evaluated (here for y component)
                //
                const auto& frontalScvfOnBoundary = fvGeometry.frontalScvfOnBoundary(scv);
                assert(frontalScvfOnBoundary.isFrontal() && frontalScvfOnBoundary.boundary());
                const auto& bcTypes = elemBcTypes[frontalScvfOnBoundary.localIndex()];
                if (bcTypes.hasNeumann() && bcTypes.isNeumann(scvf.normalAxis()))
                {
                    const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                    return neumannFluxes[scv.dofAxis()] * Extrusion::area(fvGeometry, scvf) * elemVolVars[scv].extrusionFactor();
                }
            }
        }
        else if (scvf.boundary())
        {
            // Here, the lateral face does lie on a boundary. Retrieve the Neumann flux component for
            // the corresponding scvs' DOF orientation.
            //
            //
            //     *****************
            //     ---------######O                 || frontal face of staggered half-control-volume
            //     |      ||      :
            //     |      ||      :                 D  dof position
            //     |      || scv  D
            //     |      ||      :                -- element boundaries
            //     |      ||      :
            //     ----------------                 * domain boundary
            //
            //  y                                    # current scvf over which Neumann flux is evaluated
            //  ^
            //  |                                    : frontal scvf on boundary
            //  ----> x
            //                                       O integration point at which Neumann flux is evaluated (here for x component)
            //
            assert(scvf.isLateral());
            const auto& bcTypes = elemBcTypes[scvf.localIndex()];
            if (bcTypes.hasNeumann() && bcTypes.isNeumann(scv.dofAxis()))
            {
                const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                return neumannFluxes[scv.dofAxis()] * Extrusion::area(fvGeometry, scvf) * elemVolVars[scv].extrusionFactor();
            }
        }

        // Default: Neither the scvf itself nor a lateral one considers a Neumann flux. We just calculate the flux as normal.
        return this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache, elemBcTypes);
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
