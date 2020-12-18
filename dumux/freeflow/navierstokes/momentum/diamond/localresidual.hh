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
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_DIAMOND_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_DIAMOND_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/fcdiamondlocalresidual.hh>
#include <dune/common/hybridutilities.hh>

namespace Dumux {

namespace Impl {
template<class T>
static constexpr bool isRotationalExtrusion = false;

template<int radialAxis>
static constexpr bool isRotationalExtrusion<RotationalExtrusion<radialAxis>> = true;
} // end namespace Impl

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using the staggered discretization
 */
template<class TypeTag>
class NavierStokesMomentumDiamondResidual
: public FaceCenteredDiamondLocalResidual<TypeTag>
{
    using ParentType = FaceCenteredDiamondLocalResidual<TypeTag>;
    friend class FaceCenteredDiamondLocalResidual<TypeTag>;

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
        source += problem.gravity() * problem.density(element, scv);

        // TODO Axisymmetric problems

        // // Axisymmetric problems in 2D feature an extra source terms arising from the transformation to cylindrical coordinates.
        // // See Ferziger/Peric: Computational methods for fluid dynamics chapter 8.
        // // https://doi.org/10.1007/978-3-540-68228-8 (page 301)
        // if constexpr (dim == 2 && Impl::isRotationalExtrusion<Extrusion>)
        // {
        //     if (scv.directionIndex() == Extrusion::radialAxis)
        //     {
        //         const auto r = scv.center()[scv.directionIndex()] - fvGeometry.gridGeometry().bBoxMin()[scv.directionIndex()];
        //         const auto& scvf = (*scvfs(fvGeometry, scv).begin()); // the frontal scvf belonging to the scv

        //         // Velocity term
        //         source -= -2.0*problem.effectiveViscosity(element, fvGeometry, scvf) * elemVolVars[scv].velocity() / (r*r);

        //         // Pressure term (needed because we incorporate pressure in terms of a surface integral).
        //          // grad(p) becomes div(pI) + (p/r)*n_r in cylindrical coordinates. The second term
        //          // is new with respect to Cartesian coordinates and handled below as a source term.
        //         source += problem.pressure(element, fvGeometry, scvf)/r;
        //     }
        // }

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
        assert(elemBcTypes.hasNeumann());

        NumEqVector flux = this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache, elemBcTypes);

        if (scvf.boundary())
        {
            const auto& bcTypes = elemBcTypes[scvf.localIndex()];
            const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
            for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
            {
                if (bcTypes.isNeumann(eqIdx))
                    flux[eqIdx] = neumannFluxes[eqIdx];
                else if(bcTypes.isDirichlet(eqIdx))
                    flux[eqIdx] = 0.0;
            }
        }

        return flux;
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
