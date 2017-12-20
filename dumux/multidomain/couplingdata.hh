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
 * \brief Coupling data class.
 */
#ifndef DUMUX_BOUNDARYCOUPLINGDATA_STOKES_DARCY_HH
#define DUMUX_BOUNDARYCOUPLINGDATA_STOKES_DARCY_HH

#include <type_traits>

#include "staggered-ccfv/localresidual.hh" // TODO necessary?


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
 * \brief Manages the coupling data coming from the Stokes domain
 * \ingroup BoundaryCoupling
 */
template<class TypeTag>
class StokesData
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);
    using StokesIndices = typename GET_PROP_TYPE(StokesProblemTypeTag, Indices);
    using StokesSubControlVolumeFace = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolumeFace);
    using DarcySubControlVolumeFace = typename GET_PROP_TYPE(DarcyProblemTypeTag, SubControlVolumeFace);
    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using StokesElement = typename StokesGridView::template Codim<0>::Entity;
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;
    using StokesFVElementGeometry = typename GET_PROP_TYPE(StokesProblemTypeTag, FVElementGeometry);

public:
    StokesData(const CouplingManager &couplingManager) : couplingManager_(couplingManager)
    {}

    /*!
     * \brief Returns the value of a primary variable at the position of a given Darcy vertex.
     *
     * \param darcyDofIdx The darcyDofIndex
     * \param primaryVariableIdx The index of the primary variable
     */
    Scalar valueFromSolVec(const int darcyDofIdx, const int primaryVariableIdx) const
    {
        const auto& stokesSolution  = couplingManager_.stokesProblem().model().curSol();
        const auto& stokesInfo = couplingManager_.couplingMapper().darcyToStokesMap().at(darcyDofIdx);
        const auto& stokesDofIndices = stokesInfo.stokesCCDofIndices();

        Scalar result = 0.0;
        for(const auto dofIdx : stokesDofIndices)
            result += stokesSolution[cellCenterIdx][dofIdx][primaryVariableIdx];
        return result /= stokesDofIndices.size();
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return couplingManager_; }

private:
    const CouplingManager &couplingManager_;
};

/*!
 * \brief Manages the coupling data coming from the Darcy domain
 * \ingroup BoundaryCoupling
 */
template <class TypeTag>
class DarcyData
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;
    using DarcyPrimaryVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, PrimaryVariables);
    using DarcySubControlVolume = typename GET_PROP_TYPE(DarcyProblemTypeTag, SubControlVolume);

    using DarcyFVElementGeometry = typename GET_PROP_TYPE(DarcyProblemTypeTag, FVElementGeometry);
    using DarcyElementVolumeVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, ElementVolumeVariables);

    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using StokesSubControlVolumeFace = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolumeFace);
    enum { stokesDim = StokesGridView::dimension };
    using StokesVertex = typename StokesGridView::template Codim<stokesDim>::Entity;

    using GlobalPosition = Dune::FieldVector<Scalar, stokesDim>;
    using StokesIndices = typename GET_PROP_TYPE(StokesProblemTypeTag, Indices);

    using DarcyIndices = typename GET_PROP_TYPE(DarcyProblemTypeTag, Indices);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(DarcyProblemTypeTag, ElementBoundaryTypes);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);
    using DarcyCouplingLocalResidual = typename GET_PROP_TYPE(DarcyProblemTypeTag, LocalResidual);

public:
    DarcyData(const CouplingManager &couplingManager) : couplingManager_(couplingManager)
    {}

    /*!
     * \brief Returns the value of a primary variable at the position of a given Darcy vertex.
     *
     * \param darcyDofIdx The darcyDofIndex
     * \param primaryVariableIdx The index of the primary variable
     */
    Scalar valueFromSolVec(const int stokesDofIdx, const int primaryVariableIdx) const
    {
        const auto& darcySolution  = couplingManager_.darcyProblem().model().curSol();
        const auto& darcyInfo = couplingManager_.couplingMapper().stokesCCToDarcyMap().at(stokesDofIdx);

        Scalar result = darcySolution[darcyInfo.darcyDofIdx][primaryVariableIdx];

        return result;
    }

    auto boundaryVelocity(const StokesSubControlVolumeFace &scvf, const bool verbose = false) const
    {
        const auto& darcyProblem = couplingManager().darcyProblem();
        const auto& darcyToStokesMap = couplingManager().darcyToStokesMap();
        const auto& darcyTree = darcyProblem.boundingBoxTree();

        // create a vector containing all Darcy elements coupled to the stokes SCV face
        const auto darcyCouplingInfo = couplingManager().stokesFaceToDarcyMap().at(scvf.dofIndex());

        Scalar velocity(0.0);

        const auto& darcyElement = darcyTree.entity(darcyCouplingInfo.darcyElementIdx);
        const auto darcyDofIdx = darcyCouplingInfo.darcyDofIdx;

        const auto darcyScvIdx = darcyToStokesMap.at(darcyDofIdx).darcyScvIdx();

        auto fvGeometry = localView(darcyProblem.model().globalFvGeometry());
        fvGeometry.bind(darcyElement);

        auto curElemVolVars = localView(darcyProblem.model().curGlobalVolVars());
        curElemVolVars.bind(darcyElement, fvGeometry, darcyProblem.model().curSol());
        // TODO eventually update necessary if caching enabled!
        const auto volVars = curElemVolVars[darcyScvIdx];

        auto prevElemVolVars = localView(darcyProblem.model().prevGlobalVolVars());
        prevElemVolVars.bindElement(darcyElement, fvGeometry, darcyProblem.model().prevSol());

        auto elemFluxVarsCache = localView(darcyProblem.model().globalFluxVarsCache());
        elemFluxVarsCache.bindElement(darcyElement, fvGeometry, curElemVolVars);

        ElementBoundaryTypes bcTypes;
        bcTypes.update(darcyProblem, darcyElement, fvGeometry);

        DarcyPrimaryVariables massFlux(0.0);
        DarcyCouplingLocalResidual localResidual;

        //HACK: As long as the localResidual demands problem to be non-const, the const cast has to be used here.
        // This should be ok because a new instance of localResidual is created
        // while the problem's one remains untouched. The init() and evalFluxes() functions do not have any impact one the problem itself.
        localResidual.init(const_cast<DarcyProblem &>(darcyProblem));
        localResidual.computeIntegralFluxAcrossBoundary(darcyElement, fvGeometry, prevElemVolVars, curElemVolVars, bcTypes, elemFluxVarsCache, false);

        // The flux must be substracted:
        // On an inlet boundary, the flux part of the local residual will be positive, since all fluxes will leave the SCV towards to interior domain.
        // For the domain itself, however, the sign has to be negative, since mass is entering the system.
        massFlux -= localResidual.residual(0);

//        const auto couplingArea = scvf.area() * volVars.extrusionFactor();
        // TODO scvf.area() not correctly initialized? (= 10^-310)
        const Scalar xLength = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[1];
        const Scalar xCells = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, StokesGrid, Cells0);
        const auto couplingArea = xLength / xCells;

        // Account for the orientation of the Stokes face.
        // Negative mass fluxes occur when mass is entering the Darcy domain.
        // We have to make sure the resulting velocity for the Stokes face has the correct sign.
        const Scalar directionFactor = scvf.normalInPosCoordDir() ? -1.0 : 1.0;

        velocity = massFlux / (volVars.density() * couplingArea) * directionFactor;
//        velocity = massFlux[DarcyIndices::conti0EqIdx] / (volVars.density(DarcyIndices::nPhaseIdx) * couplingArea) * directionFactor;

//        std::cout << "** couplingdata: massFlux = " << massFlux[DarcyIndices::conti0EqIdx]
//                                                    << ", rho = " << volVars.density(/*DarcyIndices::nPhaseIdx*/)
//                                                    << ", v = " << velocity
//                                                    << ", A = " << couplingArea
//                                                    << " for stokesScvfDofIdx " << scvf.dofIndex()
//                                                    << std::endl;

        return velocity;
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return couplingManager_; }

private:
    const CouplingManager &couplingManager_;
};

} // namespace Dumux

#endif // DUMUX_BOUNDARYCOUPLINGDATA_STOKES_DARCY_HH
