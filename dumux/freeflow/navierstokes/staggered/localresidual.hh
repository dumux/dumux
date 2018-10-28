// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/assembly/staggeredlocalresidual.hh>
#include <dune/common/hybridutilities.hh>
#include <dumux/freeflow/nonisothermal/localresidual.hh>

namespace Dumux
{

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class NavierStokesResidualImpl;

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using the staggered discretization
 */
template<class TypeTag>
class NavierStokesResidualImpl<TypeTag, DiscretizationMethod::staggered>
: public StaggeredLocalResidual<TypeTag>
{
    using ParentType = StaggeredLocalResidual<TypeTag>;
    friend class StaggeredLocalResidual<TypeTag>;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    using CellCenterResidual = CellCenterPrimaryVariables;
    using FaceResidual = FacePrimaryVariables;

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();
    static constexpr bool isCompositional = ModelTraits::numComponents() > 1;

public:
    using EnergyLocalResidual = FreeFlowEnergyLocalResidual<FVGridGeometry, FluxVariables, ModelTraits::enableEnergyBalance(), (ModelTraits::numComponents() > 1)>;

    // account for the offset of the cell center privars within the PrimaryVariables container
    static constexpr auto cellCenterOffset = ModelTraits::numEq() - CellCenterPrimaryVariables::dimension;
    static_assert(cellCenterOffset == ModelTraits::dim(), "cellCenterOffset must equal dim for staggered NavierStokes");

    //! Use the parent type's constructor
    using ParentType::ParentType;

    //! Evaluate fluxes entering or leaving the cell center control volume.
    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const ElementFaceVariables& elemFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        CellCenterPrimaryVariables flux = fluxVars.computeMassFlux(problem, element, fvGeometry, elemVolVars,
                                                                   elemFaceVars, scvf, elemFluxVarsCache[scvf]);

        EnergyLocalResidual::heatFlux(flux, problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf);

        return flux;
    }

    //! Evaluate the source term for the cell center control volume.
    CellCenterPrimaryVariables computeSourceForCellCenter(const Problem& problem,
                                                          const Element &element,
                                                          const FVElementGeometry& fvGeometry,
                                                          const ElementVolumeVariables& elemVolVars,
                                                          const ElementFaceVariables& elemFaceVars,
                                                          const SubControlVolume &scv) const
    {
        CellCenterPrimaryVariables result(0.0);

        // get the values from the problem
        const auto sourceValues = problem.source(element, fvGeometry, elemVolVars, elemFaceVars, scv);

        // copy the respective cell center related values to the result
        for (int i = 0; i < result.size(); ++i)
            result[i] = sourceValues[i + cellCenterOffset];

        return result;
    }


    //! Evaluate the storage term for the cell center control volume.
    CellCenterPrimaryVariables computeStorageForCellCenter(const Problem& problem,
                                                           const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage;
        storage[Indices::conti0EqIdx - ModelTraits::dim()] = volVars.density();

        EnergyLocalResidual::fluidPhaseStorage(storage, volVars);

        return storage;
    }

    //! Evaluate the storage term for the face control volume.
    FacePrimaryVariables computeStorageForFace(const Problem& problem,
                                               const SubControlVolumeFace& scvf,
                                               const VolumeVariables& volVars,
                                               const ElementFaceVariables& elemFaceVars) const
    {
        FacePrimaryVariables storage(0.0);
        const Scalar velocity = elemFaceVars[scvf].velocitySelf();
        storage[0] = volVars.density() * velocity;
        return storage;
    }

    //! Evaluate the source term for the face control volume.
    FacePrimaryVariables computeSourceForFace(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const ElementFaceVariables& elemFaceVars) const
    {
        FacePrimaryVariables source(0.0);
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        source += problem.gravity()[scvf.directionIndex()] * insideVolVars.density();
        source += problem.source(element, fvGeometry, elemVolVars, elemFaceVars, scvf)[Indices::velocity(scvf.directionIndex())];

        return source;
    }

    //! Evaluate the momentum flux for the face control volume.
    FacePrimaryVariables computeFluxForFace(const Problem& problem,
                                            const Element& element,
                                            const SubControlVolumeFace& scvf,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFaceVariables& elemFaceVars,
                                            const ElementFaceVariables& defElemFaceVars,
                                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        return fluxVars.computeMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, defElemFaceVars);
    }

    /*!
     * \brief Sets a fixed Dirichlet value for a cell (such as pressure) at the boundary.
     *        This is a provisional alternative to setting the Dirichlet value on the boundary directly.
     */
    template<class BoundaryTypes>
    void setFixedCell(CellCenterResidual& residual,
                      const Problem& problem,
                      const Element& element,
                      const SubControlVolume& insideScv,
                      const ElementVolumeVariables& elemVolVars,
                      const BoundaryTypes& bcTypes) const
    {
        // set a fixed pressure for cells adjacent to a wall
        if(bcTypes.isDirichletCell(Indices::pressureIdx))
        {
            const auto& insideVolVars = elemVolVars[insideScv];
            residual[Indices::pressureIdx - cellCenterOffset] = insideVolVars.pressure() - problem.dirichlet(element, insideScv)[Indices::pressureIdx];
        }
    }

protected:

     /*!
     * \brief Evaluate boundary conditions for a cell center dof
     */
    template<class ElementBoundaryTypes>
    void evalBoundaryForCellCenter_(CellCenterResidual& residual,
                                    const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const ElementFaceVariables& elemFaceVars,
                                    const ElementBoundaryTypes& elemBcTypes,
                                    const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
            {
                const auto bcTypes = problem.boundaryTypes(element, scvf);
                const auto extrusionFactor = elemVolVars[scvf.insideScvIdx()].extrusionFactor();

                // treat Dirichlet and outflow BCs
                FluxVariables fluxVars;
                auto boundaryFlux = fluxVars.computeMassFlux(problem, element, fvGeometry, elemVolVars,
                                                             elemFaceVars, scvf, elemFluxVarsCache[scvf]);

                EnergyLocalResidual::heatFlux(boundaryFlux, problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf);

                // treat Neumann BCs, i.e. overwrite certain fluxes by user-specified values
                static constexpr auto numEqCellCenter = CellCenterResidual::dimension;
                if(bcTypes.hasNeumann())
                {
                    for(int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
                    {
                        if(bcTypes.isNeumann(eqIdx + cellCenterOffset))
                        {
                            boundaryFlux[eqIdx] = problem.neumann(element, fvGeometry, elemVolVars, elemFaceVars, scvf)[eqIdx + cellCenterOffset]
                                                                  * extrusionFactor * scvf.area();
                        }
                    }
                }
                for(int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
                {
                    // use a wall function
                    if(problem.useWallFunction(element, scvf, eqIdx + cellCenterOffset))
                    {
                        boundaryFlux[eqIdx] = problem.wallFunction(element, fvGeometry, elemVolVars, elemFaceVars, scvf)[eqIdx]
                                                                   * extrusionFactor * scvf.area();
                    }
                }

                residual += boundaryFlux;

                // if specified, set a fixed value at the center of a cell at the boundary
                const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
                asImp_().setFixedCell(residual, problem, element, scv, elemVolVars, bcTypes);
            }
        }
    }

     /*!
     * \brief Evaluate boundary conditions for a face dof
     */
    template<class ElementBoundaryTypes>
    void evalBoundaryForFace_(FaceResidual& residual,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf,
                              const ElementVolumeVariables& elemVolVars,
                              const ElementFaceVariables& elemFaceVars,
                              const ElementFaceVariables& defElemFaceVars,
                              const ElementBoundaryTypes& elemBcTypes,
                              const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        if (scvf.boundary())
        {
            // handle the actual boundary conditions:
            const auto bcTypes = problem.boundaryTypes(element, scvf);

            if(bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())))
            {
                // set a fixed value for the velocity for Dirichlet boundary conditions
                const Scalar velocity = elemFaceVars[scvf].velocitySelf();
                const Scalar dirichletValue = problem.dirichlet(element, scvf)[Indices::velocity(scvf.directionIndex())];
                residual = velocity - dirichletValue;
            }
            else if(bcTypes.isNeumann(Indices::velocity(scvf.directionIndex())))
            {
                // set a given Neumann flux for the face on the boundary itself
                const auto extrusionFactor = elemVolVars[scvf.insideScvIdx()].extrusionFactor();
                residual += problem.neumann(element, fvGeometry, elemVolVars, elemFaceVars, scvf)[Indices::velocity(scvf.directionIndex())]
                                           * extrusionFactor * scvf.area();

                // treat the remaining (frontal and lateral) faces of the staggered control volume
                residual += computeFluxForFace(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, defElemFaceVars, elemFluxVarsCache);
            }
            else if(bcTypes.isSymmetry())
            {
                // For symmetry boundary conditions, there is no flow accross the boundary and
                // we therefore treat it like a Dirichlet boundary conditions with zero velocity
                const Scalar velocity = elemFaceVars[scvf].velocitySelf();
                const Scalar fixedValue = 0.0;
                residual = velocity - fixedValue;
            }
            else if(bcTypes.isDirichlet(Indices::pressureIdx))
            {
                // If none of the above conditions apply, we are at an "fixed pressure" boundary for which the velocity needs to be assembled
                residual += computeFluxForFace(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, defElemFaceVars, elemFluxVarsCache);
            }
            else
                DUNE_THROW(Dune::InvalidStateException, "Something went wrong with the boundary conditions for the momentum equations.");
        }
    }

private:

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif   // DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH
