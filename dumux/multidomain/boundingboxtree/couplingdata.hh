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
#ifndef DUMUX_COUPLINGDATA_STOKES_DARCY_HH
#define DUMUX_COUPLINGDATA_STOKES_DARCY_HH

#include <type_traits>


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
 * \brief Manages the coupling data coming form the stokes domain
 * \ingroup BoundaryCoupling
 */
template<class TypeTag>
class StokesData
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

    // obtain the type tags of the sub problems
    typedef typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag) StokesProblemTypeTag;
    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    StokesData(const CouplingManager &couplingManager) : couplingManager_(couplingManager)
    {}

    /*!
     * \brief Returns the value of a primary variable at the position of a given darcy vertex.
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
 * \brief Manages the coupling data coming from the low dim domain
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

public:
    DarcyData(const CouplingManager &couplingManager) : couplingManager_(couplingManager)
    {}

    /*!
     * \brief Returns the velocity \f$\mathrm{[\frac{m}{s}]}\f$ within a pore throat at a given location on the boundary
     *
     * \param element The pore throats which shall be considered for the calculation of the velocity.
     * \param scalingFactor A geometrical factor (e.g. the stokes coupling volume or area the pore is coupled with).
     * \param verbose If set true, the velocity is printed
     */
    auto boundaryVelocity(const StokesSubControlVolumeFace &scvf, const bool verbose = false)
    {
       const auto& darcyProblem = couplingManager().darcyProblem();
       const auto& darcyToStokesMap = couplingManager().darcyToStokesMap();
       const auto& darcyTree = darcyProblem.boundingBoxTree();

       auto restriction = [&darcyProblem, &darcyToStokesMap] (const DarcySubControlVolume& scv)
       {
           return darcyToStokesMap.count(scv.dofIndex());
       };

       // create a vector containing all darcy elements coupled to the stokes SCV face
       const auto darcyCouplingInfo = couplingManager().stokesFaceToDarcyMap().at(scvf.dofIndex());

       DarcyPrimaryVariables velocity(0); // TODO
//        for(auto x:darcyCouplingInfo)
//        {
//            const auto& darcyElement = darcyTree.entity(x.darcyElementIdx);
//            const auto darcyDofIdx = x.darcyDofIdx;
//
//            const auto darcyScvIdx = darcyToStokesMap.at(darcyDofIdx).darcyScvIdx();
//
//            auto fvGeometry = localView(darcyProblem.model().globalFvGeometry());
//            fvGeometry.bind(darcyElement);
//
//            auto curElemVolVars = localView(darcyProblem.model().curGlobalVolVars());
//            curElemVolVars.bind(darcyElement, fvGeometry, darcyProblem.model().curSol());
//            const auto volVars = curElemVolVars[darcyScvIdx];
//            auto massflux = Functions<DarcyProblemTypeTag>::boundaryFlux(darcyProblem, darcyElement, restriction);
//            const auto couplingArea = darcyToStokesMap.at(darcyDofIdx).couplingArea() * volVars.extrusionFactor();
//
// //            std::cout << "at lowdim : " << darcyDofIdx << ", mass flux : " << flux << " , area: " << couplingArea << std::endl;
//            velocity += massflux / (volVars.density() * couplingArea);
//        }
       return velocity; // TODO
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return couplingManager_; }

private:
    const CouplingManager &couplingManager_;
};

} // namespace Dumux

#endif
