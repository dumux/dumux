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
 * \ingroup BoundaryCoupling
 * \brief Coupling manager for two domains with the same dimension.
 *        Intersection computation relies on boundingBoxTree and therefore
          does not require grid-glue.
 */
#ifndef DUMUX_COUPLINGMANAGER_STOKES_DARCY_HH
#define DUMUX_COUPLINGMANAGER_STOKES_DARCY_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/boundingboxtree.hh>
#include "couplingmapper.hh"
#include <dumux/common/exceptions.hh>

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(StokesProblemTypeTag);
NEW_PROP_TAG(DarcyProblemTypeTag);
NEW_PROP_TAG(GridView);
} // namespace Properties

/*!
 * \brief Manages the coupling between Stokes and Darcy elements
 * \ingroup BoundaryCoupling
 */
template<class TypeTag>
class CouplingManagerStokesDarcy
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

    using StokesGrid = typename GET_PROP_TYPE(StokesProblemTypeTag, Grid);
    using DarcyGrid = typename GET_PROP_TYPE(DarcyProblemTypeTag, Grid);

    using StokesProblem = typename GET_PROP_TYPE(StokesProblemTypeTag, Problem);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);

    using DarcySubControlVolume = typename GET_PROP_TYPE(DarcyProblemTypeTag, SubControlVolume);
    using DarcySubControlVolumeFace = typename GET_PROP_TYPE(DarcyProblemTypeTag, SubControlVolumeFace);
    using StokesSubControlVolume = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolume);
    using StokesSubControlVolumeFace = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolumeFace);

    using DarcyElementVolumeVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, ElementVolumeVariables);
    using StokesElementVolumeVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, ElementVolumeVariables);

    using DarcyFVElementGeometry = typename GET_PROP_TYPE(DarcyProblemTypeTag, FVElementGeometry);
    using StokesFVElementGeometry = typename GET_PROP_TYPE(StokesProblemTypeTag, FVElementGeometry);

    using StokesGlobalFaceVars = typename GET_PROP_TYPE(StokesProblemTypeTag, GlobalFaceVars);

    using FacePrimaryVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, FacePrimaryVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, CellCenterPrimaryVariables);
    using DarcyPrimaryVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, PrimaryVariables);

    using StokesFluidSystem = typename GET_PROP_TYPE(StokesProblemTypeTag, FluidSystem);
    using DarcyFluidSystem = typename GET_PROP_TYPE(DarcyProblemTypeTag, FluidSystem);

    using StokesFluidState = typename GET_PROP_TYPE(StokesProblemTypeTag, FluidState); // TODO ?

    enum {
        dim = StokesGridView::dimension,
        dimWorld = StokesGridView::dimensionworld
    };

    using StokesElement = typename StokesGridView::template Codim<0>::Entity;
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;

    using DarcyVertex = typename DarcyGrid::template Codim<dim>::Entity;

    using CouplingMapper = Dumux::CouplingMapperStokesDarcy<TypeTag>;

    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    using cellCenterIdx = typename DofTypeIndices::CellCenterIdx;
    using faceIdx = typename DofTypeIndices::FaceIdx;

    using StokesIndices = typename GET_PROP_TYPE(StokesProblemTypeTag, Indices);
    using DarcyIndices = typename GET_PROP_TYPE(DarcyProblemTypeTag, Indices);

public:

    /*!
     * \brief Constructor
     */
    CouplingManagerStokesDarcy(StokesProblem& stokesProblem, DarcyProblem& darcyProblem)
    : stokesProblem_(stokesProblem),
      darcyProblem_(darcyProblem),
      couplingMapper_(stokesProblem, darcyProblem, asImp_())
    {
        massStokesToDarcy_ = 0.0;
        energyStokesToDarcy_ = 0.0;
    }


    /*!
     * \brief Return a reference to the stokes problem
     */
    const StokesProblem& stokesProblem() const
    {
        return stokesProblem_;
    }

    /*!
     * \brief Return a reference to the low dimensional problem
     */
    const DarcyProblem& darcyProblem() const
    {
        return darcyProblem_;
    }

    /*!
     * \brief Return a reference to the stokes gridview
     */
    const StokesGridView& stokesGridView() const
    {
        return stokesProblem().gridView();
    }

    /*!
     * \brief Return a reference to the low dimensional gridview
     */
    const DarcyGridView& darcyGridView() const
    {
        return darcyProblem().gridView();
    }

    void preInit()
    {}


    /*!
     * \brief Compute the coupling maps and stencil after the sub-problem have been initialized
     */
    void postInit()
    {
        couplingMapper_.computeCouplingMaps();
        computeStencils();
    }

    /*!
     * \brief Returns whether a stokes element sub control volume is considered for coupling
     *        with the darcy domain. This function is used for the boundaryType() method
     *        of the stokes problem which handles vertices.
     *
     * \param element The element
     */
    bool isStokesCouplingEntity(const StokesElement &element) const
    {
        return !(couplingStencil(element).empty());
    }


    /*!
     * \brief Returns whether a stokes element sub control volume is considered for coupling
     *        with the darcy domain. This function is used for the boundaryType() method
     *        of the stokes problem which handles vertices.
     *
     * \param vertex The vertex
     */
    bool isStokesCouplingEntity(const StokesElement &element, const StokesSubControlVolumeFace& scvf) const
    {
        return couplingMapper_.stokesFaceToDarcyMap().count(scvf.dofIndex());
    }

    /*!
     * \brief Returns whether a low dim vertex is considered for coupling
     *        with the stokes domain. This function is used for the boundaryType() method
     *        of the low dim problem which handles vertices.
     *
     * \param vertex The vertex
     */
    bool isDarcyCouplingEntity(const DarcyVertex &vertex) const
    {
        const auto darcyDofIdxGlobal = darcyGridView().indexSet().index(vertex);
        return couplingMapper_.darcyToStokesMap().count(darcyDofIdxGlobal);
    }

    /*!
     * \brief Returns whether a low dim vertex is considered for coupling
     *        with the stokes domain. This function is used for the boundaryType() method
     *        of the low dim problem which handles vertices.
     *
     * \param element The element
     */
    bool isDarcyCouplingEntity(const DarcyElement &element) const
    {
        return !(couplingStencil(element).empty());
    }

    /*!
     * \brief Returns whether a low dim vertex is considered for coupling
     *        with the stokes domain. This function is used for the boundaryType() method
     *        of the low dim problem which handles vertices.
     *
     * \param element The element
     */
    bool isDarcyCouplingEntity(const DarcySubControlVolume &scv) const
    {
        return couplingMapper_.darcyToStokesMap().count(scv.dofIndex());
    }

    /*!
     * \brief Returns whether a low dim vertex is considered for coupling
     *        with the stokes domain. This function is used for the boundaryType() method
     *        of the low dim problem which handles vertices.
     *
     * \param element The element
     */
    bool isDarcyCouplingEntity(const DarcySubControlVolumeFace &scvf) const
    {
        return couplingMapper_.darcyToStokesMap().count(scvf.insideScvIdx()); //dofIndex()); TODO
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const CouplingMapper &couplingMapper() const
    {
        return couplingMapper_;
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const auto &darcyToStokesMap() const
    {
        return couplingMapper_.darcyToStokesMap();
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const auto &stokesFaceToDarcyMap() const
    {
        return couplingMapper_.stokesFaceToDarcyMap();
    }

    void computeStencils()
    {
        const auto &stokesTree = stokesProblem_.boundingBoxTree();

        // compute Stokes cell-center coupling stencil using the coupling map
        for (const auto &entry : couplingMapper_.stokesCCToDarcyMap())
        {
            const auto stokesElementIdx = entry.first;
            Dune::dverb << "Stokes element " << stokesElementIdx <<  " is coupled to Darcy element:";

            // get the Darcy DOFs associated with the Stokes element
            auto &darcyInfo = couplingMapper_.stokesCCToDarcyMap().at(stokesElementIdx);
            auto darcyElementIndex = darcyInfo.darcyElementIdx;
            Dune::dverb << " " << darcyElementIndex ;
            stokesCCCouplingStencils_[stokesElementIdx].push_back(darcyElementIndex);

            const auto& stokesElement = stokesTree.entity(stokesElementIdx);
            StokesFVElementGeometry stokesFvGeometry = localView(stokesProblem_.model().globalFvGeometry());
            stokesFvGeometry.bind(stokesElement);

            for(auto&& stokesScvf : scvfs(stokesFvGeometry))
                stokesFaceCouplingStencils_[stokesScvf.dofIndex()].push_back(darcyElementIndex);

            Dune::dverb << std::endl;
        }

        Dune::dverb << std::endl;

        // compute Darcy coupling stencil using the coupling map
        for (const auto &entry : couplingMapper_.darcyToStokesMap())
        {
            const auto darcyElementIdx = entry.first;
            Dune::dverb << "Darcy element " << darcyElementIdx <<  " is coupled to Stokes element:";

            // get the Stokes DOFs associated with the Darcy element
            const auto &stokesElements = couplingMapper_.darcyToStokesMap().at(darcyElementIdx).stokesElementIndices();
            auto stokesElementIndex = stokesElements[0];
            Dune::dverb << " Stokes element " << stokesElementIndex;
            darcyToCCCouplingStencils_[darcyElementIdx].push_back(stokesElementIndex);

            const auto& stokesElement = stokesTree.entity(stokesElementIndex);
            StokesFVElementGeometry stokesFvGeometry = localView(stokesProblem_.model().globalFvGeometry());
            stokesFvGeometry.bind(stokesElement);

            for(auto&& stokesScvf : scvfs(stokesFvGeometry))
                darcyToFaceCouplingStencils_[darcyElementIdx].push_back(stokesScvf.dofIndex());

            Dune::dverb << std::endl;
        }

        auto sortAndMakeUnique = [](auto&& stencils)
        {
            for(auto&& stencil : stencils)
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            }
        };

        // sort and make unique
        sortAndMakeUnique(stokesCCCouplingStencils_);
        sortAndMakeUnique(stokesFaceCouplingStencils_);
        sortAndMakeUnique(darcyToCCCouplingStencils_);
        sortAndMakeUnique(darcyToFaceCouplingStencils_);
    }

    /*!
     * \brief Returns a coupling stencil for a Stokes element
     */
    const std::vector<unsigned int>& couplingStencil(const StokesElement& element) const
    {
        const unsigned int eIdx = stokesProblem().elementMapper().index(element);
        if (stokesCCCouplingStencils_.count(eIdx))
            return stokesCCCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Stokes element
     */
    const std::vector<unsigned int>& couplingStencil(const StokesSubControlVolumeFace& scvf) const
    {
        const unsigned int dofIdx = scvf.dofIndex();
        if (stokesFaceCouplingStencils_.count(dofIdx))
            return stokesFaceCouplingStencils_.at(dofIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Darcy element
     */
    const std::vector<unsigned int>& couplingStencil(const DarcyElement& element, cellCenterIdx) const
    {
        const unsigned int eIdx = darcyProblem().elementMapper().index(element);
        if (darcyToCCCouplingStencils_.count(eIdx))
            return darcyToCCCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Darcy element
     */
    const std::vector<unsigned int>& couplingStencil(const DarcyElement& element, faceIdx) const
    {
        const unsigned int eIdx = darcyProblem().elementMapper().index(element);
        if (darcyToFaceCouplingStencils_.count(eIdx))
            return darcyToFaceCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    //! evaluate coupling residual for the derivative Stokes DOF with respect to Darcy DOF
    auto evalStokesCCCouplingResidual(const StokesElement& element, // TODO not needed
                              const StokesFVElementGeometry& fvGeometry,
                              const StokesElementVolumeVariables& elemVolVars,
                              const StokesGlobalFaceVars& globalFaceVars)
    {
        Scalar couplingScvfIdx = couplingSubControlVolumeFace(fvGeometry);
        Scalar outerNormalScalar = fvGeometry.scvf(couplingScvfIdx).outerNormalScalar();
        Scalar interfaceArea = fvGeometry.scvf(couplingScvfIdx).area();
        // TODO volVars = elemVolVars[scv], only 1 scv --> [0]?
        Scalar densityStokes = elemVolVars[0].density(0); // TODO upwind? // TODO phaseIdx = 0???
        const Scalar velocity = globalFaceVars.faceVars(couplingScvfIdx).velocity();

        Scalar diffusionCoefficient = StokesFluidSystem::binaryDiffusionCoefficient(fluidState, StokesFluidSystem::nPhaseIdx, StokesFluidSystem::wCompIdx, StokesFluidSystem::nCompIdx);
        Scalar yPoscenterStokesElement = fvGeometry.scv(fvGeometry.scvf(couplingScvfIdx).insideScvIdx()).center()[1];
        Scalar yPoscenterDarcyElement = centerInDarcyElement(fvGeometry.scvf(couplingScvfIdx))[1];
        Scalar distanceCenters = yPoscenterStokesElement - yPoscenterDarcyElement ; // TODO

        // TODO stokesProblem_.massBalanceIdx, temperature, components
        // mass coupling condition
        CellCenterPrimaryVariables stokesCCCouplingResidual(0.0);

        std::cout << "** couplingmanager: Stokes numPhases = " << StokesFluidSystem::numPhases << ", numComponents = " << StokesFluidSystem::numComponents
                << ", Darcy numPhases = " << DarcyFluidSystem::numPhases << ", numComponents = " << DarcyFluidSystem::numComponents << std::endl;
        // 1p
        if (StokesFluidSystem::numComponents == 1) // TODO numphases?
        {
            stokesCCCouplingResidual[StokesIndices::pressureIdx] = densityStokes * velocity * outerNormalScalar * interfaceArea;
        }
        // 2p2c
        else if (StokesFluidSystem::numComponents == 2) // TODO numphases?
        {
            // TODO check each variable, rough draft/attempt
            StokesFluidState fluidState;
            Scalar massFractionVapor = elemVolVars[0].massFraction(StokesFluidSystem::nPhaseIdx, StokesFluidSystem::wCompIdx);
            Scalar moleFracWaterInGasDarcy = moleFracInDarcyElement(fvGeometry.scvf(couplingScvfIdx));
            Scalar moleFracStokes = elemVolVars[0].moleFraction(StokesFluidSystem::nPhaseIdx, StokesFluidSystem::wCompIdx); // TODO correct???
            Scalar deltaX = moleFracWaterInGasDarcy - moleFracStokes;
            Scalar gradientMassFractionVapor = deltaX / distanceCenters;
            Scalar diffusiveFlux =  -1.0 * diffusionCoefficient * densityStokes * gradientMassFractionVapor; // j_g^wff = -D_g * rho_g * grad X_g^w
            stokesCCCouplingResidual[StokesIndices::massOrMoleFracIdx] = (densityStokes * massFractionVapor * velocity + diffusiveFlux) * outerNormalScalar * interfaceArea;
        }
        else // TODO error
        {
            std::cout << "** WARNING: coupling not implemented for more than two phases!" << std::endl;
        }

        // TODO
        // if (EnableEnergyBalance)
        {
            Scalar energyFlux(0.0);
            // advective fluxes
            Scalar enthalpyStokes = elemVolVars[0].enthalpy(StokesFluidSystem::nPhaseIdx); // TODO upwinding!
            energyFlux += densityStokes * enthalpyStokes * velocity;

            // diffusive fluxes (for each component, here: water, air)
            Scalar diffusiveFlux = -1.0 * enthalpyStokes * diffusionCoefficient * densityStokes;

            Scalar moleFracAirInGasDarcy = moleFracAirInDarcyElement(fvGeometry.scvf(couplingScvfIdx));
            Scalar moleFracAirInStokes = elemVolVars[0].moleFraction(StokesFluidSystem::wPhaseIdx, StokesFluidSystem::nCompIdx); // TODO correct???
            Scalar deltaXAir = moleFracAirInGasDarcy - moleFracAirInStokes;
            Scalar gradMassFractionAir = deltaXAir / distanceCenters;

            Scalar gradMassFractionAir = (elemVolVars[0].massFraction) / distanceCenters;
            Scalar moleFracWaterInGasDarcy = moleFracWaterInDarcyElement(fvGeometry.scvf(couplingScvfIdx));
            Scalar moleFracWaterInStokes = elemVolVars[0].moleFraction(StokesFluidSystem::nPhaseIdx, StokesFluidSystem::wCompIdx); // TODO correct???
            Scalar deltaXWater = moleFracWaterInGasDarcy - moleFracWaterInStokes;
            Scalar gradMassFractionWater = deltaXWater / distanceCenters;

            energyFlux += diffusiveFlux * (gradMassFractionAir + gradMassFractionWater);

            // TODO compare implem. Fetzer2017 (see following) with equation
//            Scalar diffusiveFluxesAir = -1.0 * diffusiveFluxesWater * BaseFluid::molarMassComponent(StokesIndices::phaseCompIdx)
//                                             / BaseFluid::molarMassComponent(StokesIndices::transportCompIdx);
//            Scalar enthalpyComponentWater_up = enthalpyComponentWater_s; // Stokes
//            Scalar enthalpyComponentAir_up = enthalpyComponentAir_s; // Stokes
//
//            if (diffusiveFluxesWater < 0)
//            {
//                enthalpyComponentWater_up = enthalpyComponentWater_n; // Darcy
//                enthalpyComponentAir_up = enthalpyComponentAir_n; // Darcy
//            }
//
//            Scalar diffusiveEnthalpyFluxAir = enthalpyComponentAir_up * diffusiveFluxesAir_s;
//            Scalar diffusiveEnthalpyFluxWater = enthalpyComponentWater_up * diffusiveFluxesWater_s;

            // conductive fluxes
            Scalar thermalConductivity = elemVolVars[0].thermalConductivity(); // TODO
            Scalar tempDarcy = tempInDarcyElement(fvGeometry.scvf(couplingScvfIdx));
            Scalar gradTemp = (elemVolVars[0].temperature(StokesFluidSystem::nPhaseIdx) - tempDarcy) / distanceCenters;
            energyFlux += -1.0 * thermalConductivity * gradTemp * outerNormalScalar; // TODO = implem. Fetzer2017a, in eq. "whole flux * n"

            stokesCCCouplingResidual[StokesIndices::temperatureIdx]  = energyFlux * interfaceArea;
        }

//        Scalar elemIdx = fvGeometry.scv(fvGeometry.scvf(couplingScvfIdx).insideScvIdx()).elementIndex();
//        std::cout << "** couplingmanager: StokesCCCouplingResidual = " << stokesCCCouplingResidual[0] << " at cc in element " << elemIdx<< std::endl;

        massStokesToDarcy_ = stokesCCCouplingResidual[StokesIndices::pressureIdx];
        energyStokesToDarcy_ = stokesCCCouplingResidual[StokesIndices::temperatureIdx];

        return stokesCCCouplingResidual; // fluxStokesCCToDarcy
    }

    //! evaluate coupling residual for the derivative stokes DOF with respect to low dim DOF
    auto evalStokesFaceCouplingResidual(const StokesElement& element,
                              const StokesFVElementGeometry& fvGeometry,
                              const StokesGlobalFaceVars& globalFaceVars,
                              const StokesSubControlVolumeFace& scvf)
    {
        FacePrimaryVariables stokesFaceCouplingResidual(0.0);

        Scalar couplingScvfIdx = scvf.index();
        Scalar outerNormalScalar = fvGeometry.scvf(couplingScvfIdx).outerNormalScalar();
        Scalar interfaceArea = fvGeometry.scvf(couplingScvfIdx).area();

        if (scvf.directionIndex() == 1 && stokesProblem_.onCouplingInterface(scvf.center())) // horizontal scvf on coupling interface
        {
            Scalar pressureDarcy = pressureInDarcyElement(scvf);

            // normal momentum coupling condition
            stokesFaceCouplingResidual = pressureDarcy * outerNormalScalar * interfaceArea;
        }

        // TODO only for 2D, add loop over tangDims
        // Fetzer2017 "classical BJ condition"
        // tangential momentum coupling condition --> b.c. for stokesProblem?

        else if (scvf.directionIndex() == 0 && !scvf.boundary()) // vertical scvf, not on boundary
        { // TODO works only for horizontal coupling interface on lower boundary of Stokes domain
            if( (scvf.index() % 4 == 0 && stokesProblem_.onCouplingInterface(fvGeometry.scvf(scvf.index()+2).center())) // if left scvf --> right neighbor scvf on interface
             || (scvf.index() % 4 == 1 && stokesProblem_.onCouplingInterface(fvGeometry.scvf(scvf.index()+1).center())) ) // if right scvf --> left neighbor scvf on interface
            {

                Scalar distanceCenters = (fvGeometry.scv(scvf.insideScvIdx()).center() - scvf.center())[0]; // 0.25 TODO direction? 0,1,2?
                Scalar permeability = darcyProblem_.spatialParams().permeabilityAtPos(scvf.center()); // 1e-10
                Scalar alphaBeaversJoseph = darcyProblem_.spatialParams().beaversJosephCoeffAtPos(scvf.center()); // 1
                Scalar beta = -1.0 * std::sqrt(permeability) / alphaBeaversJoseph / distanceCenters * outerNormalScalar; // 4e-5

                auto stokesElemVolVars = localView(stokesProblem_.model().curGlobalVolVars());
                stokesElemVolVars.bind(element, fvGeometry, stokesProblem_.model().curSol());
                const auto stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];

                Scalar dynViscosity = stokesVolVars.viscosity();
//                Scalar effDynViscosity = BaseFluid::dynamicViscosity(stokesVolVars.pressure(), stokesVolVars.temperature(), stokesVolVars.massFrac());
//                                                                + (eddyKinematicViscosity_(elementInsideIdx) * stokesVolVars.density()); TODO
                const Scalar tangXVelocity = globalFaceVars.faceVars(scvf.dofIndex()).velocity();
                Scalar beaversJosephXVelocity = beta * tangXVelocity / (1.0 + beta);
//                Scalar tangZVelocity = faceSolutionVector[velocityZIdx];
//                Scalar beaversJosephZVelocity = beta * tangZVelocity / (1.0 + beta);

                stokesFaceCouplingResidual = -0.5 * dynViscosity * (tangXVelocity - beaversJosephXVelocity)
                    / distanceCenters * outerNormalScalar * interfaceArea;

//                stokesFaceCouplingResidual = -0.5 * effDynViscosity * (tangZVelocity - beaversJosephZVelocity)
//                                                                          / distanceCenters * outerNormalScalar * interfaceArea;
            }
        }
//        if (stokesFaceCouplingResidual != 0.0)
//            std::cout << "** couplingmanager: StokesFaceCouplingResidual = " << stokesFaceCouplingResidual<< " at scvf " <<  scvf.index() << std::endl;

        return stokesFaceCouplingResidual;
    }

    //! evaluate coupling residual for the derivative Darcy DOF with respect to Stokes DOF
    auto evalDarcyCouplingResidual(const DarcyElement& element)
    {
        DarcyPrimaryVariables fluxDarcyToStokes(0.0);
        fluxDarcyToStokes[DarcyIndices::pressureIdx] = -1.0 * massStokesToDarcy_[StokesIndices::pressureIdx];
//        // 2p2c
//        fluxDarcyToStokes[DarcyIndices::contiWEqIdx] = -1.0 * fluxStokesToDarcy_[StokesIndices::massOrMoleFracIdx];
//        // temperature
//        fluxDarcyToStokes[DarcyIndices::temperatureIdx] = -1.0 * fluxStokesToDarcy_[StokesIndices::temeratureIdx];

        // TODO if (non-isothermal) // temperatureIdx not available for isothermal model!
        fluxDarcyToStokes[DarcyIndices::temperatureIdx] = -1.0 * energyStokesToDarcy_[StokesIndices::temperatureIdx];

        return fluxDarcyToStokes; // fluxDarcyToStokesCC
    }

protected:
    // ! Returns the global index of the subcontrolvolumeface at the coupling interface
    const Scalar couplingSubControlVolumeFace(const StokesFVElementGeometry& fvGeometry)
    {
        Scalar couplingScvfIdx(-1);
        // assumption: only one interface scvf (no corners in coupling interface)
        for(const auto &scvf : scvfs(fvGeometry))
        {
            if(stokesProblem_.onCouplingInterface(scvf.center()))
                couplingScvfIdx = scvf.index();
        }
        return couplingScvfIdx;
    }

    // ! Returns the pressure value of the adjacent Darcy element
    const Scalar pressureInDarcyElement(const StokesSubControlVolumeFace& scvf) // TODO remove duplicate code to find Darcy element
    {
        const auto& darcyTree = darcyProblem_.boundingBoxTree();

        // create a vector containing all Darcy elements coupled to the Stokes scvf
        const auto darcyCouplingInfo = stokesFaceToDarcyMap().at(scvf.dofIndex());

        const auto& darcyCouplingElement = darcyTree.entity(darcyCouplingInfo.darcyElementIdx);
        DarcyFVElementGeometry darcyFvGeometry = localView(darcyProblem_.model().globalFvGeometry());
        darcyFvGeometry.bind(darcyCouplingElement);

        const auto darcyDofIdx = darcyCouplingInfo.darcyDofIdx;

        auto darcyElemVolVars = localView(darcyProblem_.model().curGlobalVolVars());
        darcyElemVolVars.bind(darcyCouplingElement, darcyFvGeometry, darcyProblem_.model().curSol());
        const auto darcyVolVars = darcyElemVolVars[darcyDofIdx];

        return darcyVolVars.pressure(0); // TODO phaseIdx?!!!
    }

    const Scalar moleFracWaterInDarcyElement(const StokesSubControlVolumeFace& scvf) // TODO remove duplicate code to find Darcy element
    {
        const auto& darcyTree = darcyProblem_.boundingBoxTree();

        // create a vector containing all Darcy elements coupled to the Stokes scvf
        const auto darcyCouplingInfo = stokesFaceToDarcyMap().at(scvf.dofIndex());

        const auto& darcyCouplingElement = darcyTree.entity(darcyCouplingInfo.darcyElementIdx);
        DarcyFVElementGeometry darcyFvGeometry = localView(darcyProblem_.model().globalFvGeometry());
        darcyFvGeometry.bind(darcyCouplingElement);

        const auto darcyDofIdx = darcyCouplingInfo.darcyDofIdx;

        auto darcyElemVolVars = localView(darcyProblem_.model().curGlobalVolVars());
        darcyElemVolVars.bind(darcyCouplingElement, darcyFvGeometry, darcyProblem_.model().curSol());
        const auto darcyVolVars = darcyElemVolVars[darcyDofIdx];

        return darcyVolVars.moleFraction(DarcyIndices::nPhaseIdx, DarcyIndices::wCompIdx);
    }

    const Scalar moleFracAirInDarcyElement(const StokesSubControlVolumeFace& scvf) // TODO remove duplicate code to find Darcy element
    {
        const auto& darcyTree = darcyProblem_.boundingBoxTree();

        // create a vector containing all Darcy elements coupled to the Stokes scvf
        const auto darcyCouplingInfo = stokesFaceToDarcyMap().at(scvf.dofIndex());

        const auto& darcyCouplingElement = darcyTree.entity(darcyCouplingInfo.darcyElementIdx);
        DarcyFVElementGeometry darcyFvGeometry = localView(darcyProblem_.model().globalFvGeometry());
        darcyFvGeometry.bind(darcyCouplingElement);

        const auto darcyDofIdx = darcyCouplingInfo.darcyDofIdx;

        auto darcyElemVolVars = localView(darcyProblem_.model().curGlobalVolVars());
        darcyElemVolVars.bind(darcyCouplingElement, darcyFvGeometry, darcyProblem_.model().curSol());
        const auto darcyVolVars = darcyElemVolVars[darcyDofIdx];

        return darcyVolVars.moleFraction(DarcyIndices::wPhaseIdx, DarcyIndices::nCompIdx);
    }

    const Scalar tempInDarcyElement(const StokesSubControlVolumeFace& scvf) // TODO remove duplicate code to find Darcy element
    {
        const auto& darcyTree = darcyProblem_.boundingBoxTree();

        // create a vector containing all Darcy elements coupled to the Stokes scvf
        const auto darcyCouplingInfo = stokesFaceToDarcyMap().at(scvf.dofIndex());

        const auto& darcyCouplingElement = darcyTree.entity(darcyCouplingInfo.darcyElementIdx);
        DarcyFVElementGeometry darcyFvGeometry = localView(darcyProblem_.model().globalFvGeometry());
        darcyFvGeometry.bind(darcyCouplingElement);

        const auto darcyDofIdx = darcyCouplingInfo.darcyDofIdx;

        auto darcyElemVolVars = localView(darcyProblem_.model().curGlobalVolVars());
        darcyElemVolVars.bind(darcyCouplingElement, darcyFvGeometry, darcyProblem_.model().curSol());
        const auto darcyVolVars = darcyElemVolVars[darcyDofIdx];

        return darcyVolVars.temperature(DarcyIndices::nPhaseIdx);
    }

    const auto centerInDarcyElement(const StokesSubControlVolumeFace& scvf) // TODO remove duplicate code to find Darcy element
    {
        const auto& darcyTree = darcyProblem_.boundingBoxTree();

        // create a vector containing all Darcy elements coupled to the Stokes scvf
        const auto darcyCouplingInfo = stokesFaceToDarcyMap().at(scvf.dofIndex());

        const auto& darcyCouplingElement = darcyTree.entity(darcyCouplingInfo.darcyElementIdx);
        DarcyFVElementGeometry darcyFvGeometry = localView(darcyProblem_.model().globalFvGeometry());
        darcyFvGeometry.bind(darcyCouplingElement);

        auto center = darcyFvGeometry.scv(scvf.insideScvIdx()).center();
        return center;
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

private:
    StokesProblem& stokesProblem_;
    DarcyProblem& darcyProblem_;
    CouplingMapper couplingMapper_;

    DarcyPrimaryVariables massStokesToDarcy_;
    DarcyPrimaryVariables energyStokesToDarcy_

    std::unordered_map<unsigned int, std::vector<unsigned int> > stokesCCCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > stokesFaceCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > darcyToCCCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > darcyToFaceCouplingStencils_;
    std::vector<unsigned int> emptyStencil_;
};

} // namespace Dumux

#endif
