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
template<class TypeTag, DiscretizationMethods Method>
class NavierStokesResidualImpl;

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using the staggered discretization
 */
template<class TypeTag>
class NavierStokesResidualImpl<TypeTag, DiscretizationMethods::Staggered>
: public StaggeredLocalResidual<TypeTag>
{
    using ParentType = StaggeredLocalResidual<TypeTag>;
    friend class StaggeredLocalResidual<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);

    using CellCenterResidual = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceResidual = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    static constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);

    enum {
        pressureIdx = Indices::pressureIdx,

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    static constexpr bool normalizePressure = GET_PROP_VALUE(TypeTag, NormalizePressure);

public:

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

        computeFluxForCellCenterNonIsothermal_(std::integral_constant<bool, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>(),
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
        const auto sourceValues = problem.sourceAtPos(scv.center());

        // copy the respective cell center related values to the result
        for(int i = 0; i < numEqCellCenter; ++i)
            result[i] = sourceValues[i];

        return result;
    }


    //! Evaluate the storage term for the cell center control volume.
    CellCenterPrimaryVariables computeStorageForCellCenter(const Problem& problem,
                                                           const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage;
        storage[massBalanceIdx] = volVars.density();

        computeStorageForCellCenterNonIsothermal_(std::integral_constant<bool, GET_PROP_VALUE(TypeTag, EnableEnergyBalance) >(),
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
                                              const SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const ElementFaceVariables& elementFaceVars) const
    {
        FacePrimaryVariables source(0.0);
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        source += problem.gravity()[scvf.directionIndex()] * insideVolVars.density();

        source += problem.sourceAtPos(scvf.center())[Indices::velocity(scvf.directionIndex())];

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
        FacePrimaryVariables flux(0.0);
        FluxVariables fluxVars;
        flux += fluxVars.computeNormalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elementFaceVars);
        flux += fluxVars.computeTangetialMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, elementFaceVars);
        flux += computePressureTerm_(problem, element, scvf, fvGeometry, elemVolVars, elementFaceVars);
        return flux;
    }

protected:

     /*!
     * \brief Evaluate boundary conditions
     */
    void evalBoundary_(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const ElementFaceVariables& elemFaceVars,
                       const ElementBoundaryTypes& elemBcTypes,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        evalBoundaryForCellCenter_(element, fvGeometry, elemVolVars, elemFaceVars, elemBcTypes, elemFluxVarsCache);
        for (auto&& scvf : scvfs(fvGeometry))
        {
            evalBoundaryForFace_(element, fvGeometry, scvf, elemVolVars, elemFaceVars, elemBcTypes, elemFluxVarsCache);
        }
    }

     /*!
     * \brief Evaluate boundary conditions for a cell center dof
     */
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

                // handle the actual boundary conditions:
                const auto bcTypes = problem.boundaryTypes(element, scvf);

                if(bcTypes.hasNeumann())
                {
                    // handle Neumann BCs, i.e. overwrite certain fluxes by user-specified values
                    for(int eqIdx = 0; eqIdx < GET_PROP_VALUE(TypeTag, NumEqCellCenter); ++eqIdx)
                        if(bcTypes.isNeumann(eqIdx))
                        {
                            const auto extrusionFactor = 1.0; //TODO: get correct extrusion factor
                            boundaryFlux[eqIdx] = problem.neumann(element, fvGeometry, elemVolVars, scvf)[eqIdx]
                                                   * extrusionFactor * scvf.area();
                        }
                }

                residual += boundaryFlux;

                asImp_().setFixedCell_(residual, problem, fvGeometry.scv(scvf.insideScvIdx()), elemVolVars, bcTypes);
            }
        }
    }

    /*!
     * \brief Sets a fixed Dirichlet value for a cell (such as pressure) at the boundary.
     *        This is a provisional alternative to setting the Dirichlet value on the boundary directly.
     */
    void setFixedCell_(CellCenterResidual& residual,
                       const Problem& problem,
                       const SubControlVolume& insideScv,
                       const ElementVolumeVariables& elemVolVars,
                       const BoundaryTypes& bcTypes) const
    {
        // set a fixed pressure for cells adjacent to a wall
        if(bcTypes.isDirichletCell(massBalanceIdx))
        {
            const auto& insideVolVars = elemVolVars[insideScv];
            residual[pressureIdx] = insideVolVars.pressure() - problem.dirichletAtPos(insideScv.center())[pressureIdx];
        }
    }

     /*!
     * \brief Evaluate boundary conditions for a face dof
     */
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
            if(bcTypes.isDirichlet(momentumBalanceIdx))
            {
                const Scalar velocity = elementFaceVars[scvf].velocitySelf();
                const Scalar dirichletValue = problem.dirichlet(element, scvf)[Indices::velocity(scvf.directionIndex())];
                residual = velocity - dirichletValue;
            }

            // For symmetry boundary conditions, there is no flow accross the boundary and
            // we therefore treat it like a Dirichlet boundary conditions with zero velocity
            if(bcTypes.isSymmetry())
            {
                // const Scalar velocity = faceVars.faceVars(scvf.dofIndex()).velocity();
                const Scalar velocity = elementFaceVars[scvf].velocitySelf();
                const Scalar fixedValue = 0.0;
                residual = velocity - fixedValue;
            }

            // outflow condition for the momentum balance equation
            if(bcTypes.isOutflow(momentumBalanceIdx))
            {
                if(bcTypes.isDirichlet(massBalanceIdx))
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

        flux[Indices::energyBalanceIdx] = FluxVariables::advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTerm, isOutflow);
        flux[Indices::energyBalanceIdx] += HeatConductionType::diffusiveFluxForCellCenter(problem, element, fvGeometry, elemVolVars, scvf);
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
        storage[Indices::energyBalanceIdx] = volVars.density() * volVars.internalEnergy();
    }

    //! Evaluate energy storage for isothermal models
    template <typename... Args>
    void computeStorageForCellCenterNonIsothermal_(std::false_type, Args&&... args) const {}

private:

     /*!
     * \brief Returns the pressure term
     * \param scvf The sub control volume face
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elementFaceVars The face variables
     */
    FacePrimaryVariables computePressureTerm_(const Problem& problem,
                                              const Element& element,
                                              const SubControlVolumeFace& scvf,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const ElementFaceVariables& elementFaceVars) const
    {
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        const Scalar deltaP = normalizePressure ? problem.initialAtPos(scvf.center())[pressureIdx] : 0.0;

        Scalar result = (insideVolVars.pressure() - deltaP) * scvf.area() * -1.0 * scvf.directionSign();

        // treat outflow BCs
        if(scvf.boundary())
        {
            const Scalar pressure = problem.dirichlet(element, scvf)[pressureIdx] - deltaP;
            result += pressure * scvf.area() * scvf.directionSign();
        }
        return result;
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};
}

#endif   // DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH
