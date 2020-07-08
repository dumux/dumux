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
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
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
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;


    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param scv The sub-control volume over which we integrate the storage term
     * \param volVars The volume variables associated with the scv
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
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        NumEqVector flux(0.0);
        flux += fluxVars.advectiveMomentumFlux();
        flux += fluxVars.diffusiveMomentumFlux();
        flux += fluxVars.pressureContribution();
        return flux;
    }

    //! Evaluate flux residuals for one sub control volume face when the element features at least one Neumann boundary condition
    //! This requires special care.
    NumEqVector computeFluxWithNeumannBoundaries(const Problem& problem,
                                                 const Element& element,
                                                 const FVElementGeometry& fvGeometry,
                                                 const ElementVolumeVariables& elemVolVars,
                                                 const ElementBoundaryTypes& elemBcTypes,
                                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                                 const SubControlVolumeFace& scvf) const
    {
        static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcstaggered); // TODO overload this method for different discretizations
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
                if (bcTypes.hasNeumann() && bcTypes.isNeumann(scvf.directionIndex()))
                {
                    const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                    return neumannFluxes[scvf.directionIndex()] * scvf.area() * elemVolVars[scv].extrusionFactor();
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
                for (const auto& otherScvf : scvfs(fvGeometry))
                {
                    if (otherScvf.isFrontal() && otherScvf.boundary())
                    {
                        const auto& bcTypes = elemBcTypes[otherScvf.localIndex()];
                        if (bcTypes.hasNeumann() && bcTypes.isNeumann(scvf.directionIndex()))
                        {
                            const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                            return neumannFluxes[scvf.directionIndex()] * scvf.area() * elemVolVars[scv].extrusionFactor();
                        }
                    }
                }
            }
        }
        else if (scvf.boundary())
        {
            // Here, the lateral face does lie on a boundary. Retrive the Neumann flux component for
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
            if (bcTypes.hasNeumann() && bcTypes.isNeumann(scv.directionIndex()))
            {
                const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                return neumannFluxes[scv.directionIndex()] * scvf.area() * elemVolVars[scv].extrusionFactor();
            }
        }

        // Default: Neither the scvf itself nor a lateral one considers a Neumann flux. We just calculate the flux as normal.
        return this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif   // DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH
