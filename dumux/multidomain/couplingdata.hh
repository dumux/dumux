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

double totalMassFlux = 0.0; // TODO remove

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

    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using StokesFVElementGeometry = typename GET_PROP_TYPE(StokesProblemTypeTag, FVElementGeometry);
    using StokesIndices = typename GET_PROP_TYPE(StokesProblemTypeTag, Indices);
    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);
    using DarcySubControlVolumeFace = typename GET_PROP_TYPE(DarcyProblemTypeTag, SubControlVolumeFace);

public:
    StokesData(const CouplingManager &couplingManager) : couplingManager_(couplingManager)
    {}

    /*!
     * \brief Returns the mass flux across the coupling boundary
     *
     * For the mass coupling, the free flow side of the coupling condition
     * is evalueted, i.e. [rho v]^ff.
     *
     * \param scvf The Darcy sub control volume face (on the coupling interface)
     */
    Scalar massCouplingCondition(const DarcySubControlVolumeFace& scvf) const
    {
        const auto stokesCouplingInfo = couplingManager().darcyToStokesMap().at(scvf.insideScvIdx());
        const Scalar stokesElementIndex = stokesCouplingInfo.stokesElementIndices()[0];
        const auto &stokesTree = couplingManager().stokesProblem().boundingBoxTree();
        const auto& stokesElement = stokesTree.entity(stokesElementIndex);

        StokesFVElementGeometry stokesFvGeometry = localView(couplingManager().stokesProblem().model().globalFvGeometry());
        stokesFvGeometry.bind(stokesElement);
        auto stokesElemVolVars = localView(couplingManager().stokesProblem().model().curGlobalVolVars());
        stokesElemVolVars.bind(stokesElement, stokesFvGeometry, couplingManager().stokesProblem().model().curSol());
        const auto stokesVolVars = stokesElemVolVars[stokesElementIndex]; // stokesScvIdx
        Scalar couplingDofIdx = stokesCouplingInfo.stokesFaceDofIndices()[0];

        // TODO upwind --> rho_ff / rho_pm
        const Scalar densityStokes = stokesVolVars.density(StokesIndices::phaseIdx);
        const auto& stokesSolution = couplingManager().stokesProblem().model().curSol();
        const auto velocityStokes = stokesSolution[faceIdx][couplingDofIdx];

        const Scalar massFlux = -1.0* densityStokes * velocityStokes * scvf.unitOuterNormal()[1];

        // sum up individual mass fluxes
        totalMassFlux += massFlux * scvf.area(); // multiplied by area in flux method (different file)
        std::cout << "** couplingdata: totalMassFlux = " << totalMassFlux << std::endl;

        return massFlux;
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
    using DarcyIndices = typename GET_PROP_TYPE(DarcyProblemTypeTag, Indices);

    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using StokesSubControlVolumeFace = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolumeFace);

public:
    DarcyData(const CouplingManager &couplingManager) : couplingManager_(couplingManager)
    {}

    /*!
     * \brief Returns the momentum flux across the coupling boundary
     *
     * For the momentum coupling, the porous medium side of the coupling condition
     * is evalueted, i.e. -[p n]^pm.
     *
     * \param scvf The Stokes sub control volume face (on the coupling interface)
     */
    Scalar momentumCouplingCondition(const StokesSubControlVolumeFace& scvf) const
    {
        Scalar momentumFlux(0.0);

        const auto darcyCouplingInfo = couplingManager().stokesFaceToDarcyMap().at(scvf.dofIndex());
        const auto& darcySolution = couplingManager().darcyProblem().model().curSol();
        const auto darcyPressure = darcySolution[darcyCouplingInfo.darcyElementIdx][DarcyIndices::pressureIdx];

        // - p_pm n_pm = p_pm n_ff
        momentumFlux = darcyPressure;
        momentumFlux *= scvf.outerNormalScalar() * scvf.area();

        return momentumFlux;
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return couplingManager_; }

private:
    const CouplingManager &couplingManager_;
};

} // namespace Dumux

#endif // DUMUX_BOUNDARYCOUPLINGDATA_STOKES_DARCY_HH
