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

public:

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
        CellCenterPrimaryVariables flux = fluxVars.computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars,
                                                 elemFaceVars, scvf, elemFluxVarsCache[scvf]);

        computeFluxForCellCenterNonIsothermal_(std::integral_constant<bool, ModelTraits::enableEnergyBalance()>(),
                                               problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache, flux);

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

        computeStorageForCellCenterNonIsothermal_(std::integral_constant<bool, ModelTraits::enableEnergyBalance() >(),
                                                  problem, scv, volVars, storage);

        return storage;
    }

    //! Evaluate the storage term for the face control volume.
    FacePrimaryVariables computeStorageForFace(const Problem& problem,
                                               const SubControlVolumeFace& scvf,
                                               const VolumeVariables& volVars,
                                               const ElementFaceVariables& elementFaceVars) const
    {
        FacePrimaryVariables storage(0.0);
        const Scalar velocity = elementFaceVars[scvf].velocitySelf();
        storage[0] = volVars.density() * velocity;
        return storage;
    }

    //! Evaluate the source term for the face control volume.
    FacePrimaryVariables computeSourceForFace(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const ElementFaceVariables& elementFaceVars) const
    {
        FacePrimaryVariables source(0.0);
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        source += problem.gravity()[scvf.directionIndex()] * insideVolVars.density();
        source += problem.source(element, fvGeometry, elemVolVars, elementFaceVars, scvf)[Indices::velocity(scvf.directionIndex())];

        return source;
    }

    //! Evaluate the momentum flux for the face control volume.
    FacePrimaryVariables computeFluxForFace(const Problem& problem,
                                            const Element& element,
                                            const SubControlVolumeFace& scvf,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFaceVariables& elementFaceVars,
                                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        return fluxVars.computeMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elementFaceVars);
    }

    /*!
     * \brief Sets a fixed Dirichlet value for a cell (such as pressure) at the boundary.
     *        This is a provisional alternative to setting the Dirichlet value on the boundary directly.
     */
    template<class BoundaryTypes>
    void setFixedCell(CellCenterResidual& residual,
                      const Problem& problem,
                      const SubControlVolume& insideScv,
                      const ElementVolumeVariables& elemVolVars,
                      const BoundaryTypes& bcTypes) const
    {
        // set a fixed pressure for cells adjacent to a wall
        if(bcTypes.isDirichletCell(Indices::conti0EqIdx))
        {
            const auto& insideVolVars = elemVolVars[insideScv];
            residual[Indices::conti0EqIdx - cellCenterOffset] = insideVolVars.pressure() - problem.dirichletAtPos(insideScv.center())[Indices::pressureIdx];
        }
    }

protected:

    //  /*!
    //  * \brief Evaluate boundary conditions
    //  */
    // template<class ElementBoundaryTypes>
    // void evalBoundary_(const Element& element,
    //                    const FVElementGeometry& fvGeometry,
    //                    const ElementVolumeVariables& elemVolVars,
    //                    const ElementFaceVariables& elemFaceVars,
    //                    const ElementBoundaryTypes& elemBcTypes,
    //                    const ElementFluxVariablesCache& elemFluxVarsCache)
    // {
    //     evalBoundaryForCellCenter_(element, fvGeometry, elemVolVars, elemFaceVars, elemBcTypes, elemFluxVarsCache);
    //     for (auto&& scvf : scvfs(fvGeometry))
    //     {
    //         evalBoundaryForFace_(element, fvGeometry, scvf, elemVolVars, elemFaceVars, elemBcTypes, elemFluxVarsCache);
    //     }
    // }

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
                auto boundaryFlux = computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);
                const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

                // handle the actual boundary conditions:
                const auto bcTypes = problem.boundaryTypes(element, scvf);

                if(bcTypes.hasNeumann())
                {
                    static constexpr auto numEqCellCenter = CellCenterResidual::dimension;
                    // handle Neumann BCs, i.e. overwrite certain fluxes by user-specified values
                    for(int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
                    {
                        if(bcTypes.isNeumann(eqIdx + cellCenterOffset))
                        {
                            const auto extrusionFactor = elemVolVars[scvf.insideScvIdx()].extrusionFactor();
                            boundaryFlux[eqIdx] = problem.neumann(element, fvGeometry, elemVolVars, elemFaceVars, scvf)[eqIdx + cellCenterOffset]
                                                  * extrusionFactor * scvf.area();
                        }
                    }
                }

                residual += boundaryFlux;

                asImp_().setFixedCell(residual, problem, scv, elemVolVars, bcTypes);
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
                              const ElementFaceVariables& elementFaceVars,
                              const ElementBoundaryTypes& elemBcTypes,
                              const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        if (scvf.boundary())
        {
            // handle the actual boundary conditions:
            const auto bcTypes = problem.boundaryTypes(element, scvf);

            // set a fixed value for the velocity for Dirichlet boundary conditions
            if(bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())))
            {
                const Scalar velocity = elementFaceVars[scvf].velocitySelf();
                const Scalar dirichletValue = problem.dirichlet(element, scvf)[Indices::velocity(scvf.directionIndex())];
                residual = velocity - dirichletValue;
            }

            // handle Neumann boundary conditions
            if(bcTypes.isNeumann(Indices::velocity(scvf.directionIndex())))
            {
                // set the given Neumann flux for the face on the boundary itself
                const auto extrusionFactor = elemVolVars[scvf.insideScvIdx()].extrusionFactor();
                residual = problem.neumann(element, fvGeometry, elemVolVars, elementFaceVars, scvf)[Indices::velocity(scvf.directionIndex())]
                                           * extrusionFactor * scvf.area();

                // treat the remaining (normal) faces of the staggered control volume
                FluxVariables fluxVars;
                residual += fluxVars.computeNormalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elementFaceVars);
            }

            // For symmetry boundary conditions, there is no flow accross the boundary and
            // we therefore treat it like a Dirichlet boundary conditions with zero velocity
            if(bcTypes.isSymmetry())
            {
                const Scalar velocity = elementFaceVars[scvf].velocitySelf();
                const Scalar fixedValue = 0.0;
                residual = velocity - fixedValue;
            }

            // outflow condition for the momentum balance equation
            if(bcTypes.isOutflow(Indices::velocity(scvf.directionIndex())))
            {
                if(bcTypes.isDirichlet(Indices::conti0EqIdx))
                    residual += computeFluxForFace(problem, element, scvf, fvGeometry, elemVolVars, elementFaceVars, elemFluxVarsCache);
                else
                    DUNE_THROW(Dune::InvalidStateException, "Face at " << scvf.center()  << " has an outflow BC for the momentum balance but no Dirichlet BC for the pressure!");
            }
        }
    }

    //! Evaluate energy fluxes entering or leaving the cell center control volume for non isothermal models
    void computeFluxForCellCenterNonIsothermal_(std::true_type,
                                                const Problem& problem,
                                                const Element &element,
                                                const FVElementGeometry& fvGeometry,
                                                const ElementVolumeVariables& elemVolVars,
                                                const ElementFaceVariables& elemFaceVars,
                                                const SubControlVolumeFace &scvf,
                                                const ElementFluxVariablesCache& elemFluxVarsCache,
                                                CellCenterPrimaryVariables& flux) const
    {
        // if we are on an inflow/outflow boundary, use the volVars of the element itself
        // TODO: catch neumann and outflow in localResidual's evalBoundary_()
        bool isOutflow = false;
        if(scvf.boundary())
        {
            const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
                if(bcTypes.isOutflow(Indices::energyBalanceIdx))
                    isOutflow = true;
        }

        auto upwindTerm = [](const auto& volVars) { return volVars.density() * volVars.enthalpy(); };
        using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

        flux[Indices::energyBalanceIdx - cellCenterOffset] = FluxVariables::advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTerm, isOutflow);
        flux[Indices::energyBalanceIdx - cellCenterOffset] += HeatConductionType::diffusiveFluxForCellCenter(problem, element, fvGeometry, elemVolVars, scvf);
    }

    //! Evaluate energy fluxes entering or leaving the cell center control volume for non isothermal models
    template <typename... Args>
    void computeFluxForCellCenterNonIsothermal_(std::false_type, Args&&... args) const {}

    //! Evaluate energy storage for non isothermal models
    void computeStorageForCellCenterNonIsothermal_(std::true_type,
                                                   const Problem& problem,
                                                   const SubControlVolume& scv,
                                                   const VolumeVariables& volVars,
                                                   CellCenterPrimaryVariables& storage) const
    {
        storage[Indices::energyBalanceIdx - cellCenterOffset] = volVars.density() * volVars.internalEnergy();
    }

    //! Evaluate energy storage for isothermal models
    template <typename... Args>
    void computeStorageForCellCenterNonIsothermal_(std::false_type, Args&&... args) const {}

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
