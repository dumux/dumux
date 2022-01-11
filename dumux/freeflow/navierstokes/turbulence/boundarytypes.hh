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
 * \ingroup RANSModel
 * \copydoc Dumux::RANSCCBoundaryTypes
 */
#ifndef FREEFLOW_TURBULENCE_BOUNDARY_TYPES_HH
#define FREEFLOW_TURBULENCE_BOUNDARY_TYPES_HH

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Class to specify the type of a boundary condition for the RANS extension to the Navier-Stokes model.
 */
template <class ModelTraits, int numEq>
class RANSCCBoundaryTypes : public BoundaryTypes<numEq>
{
    using ParentType = BoundaryTypes<numEq>;
    using Indices = typename ModelTraits::Indices;
public:
    RANSCCBoundaryTypes()
    {
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
            resetEq(eqIdx);
    }

    void setWall()
    {
        isWall_ = true;
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
        {
            if constexpr (numTurbulenceEq(ModelTraits::turbulenceModel()) == 1)
            {
                if (eqIdx == Indices::viscosityTildeIdx)
                    BoundaryTypes<numEq>::boundaryInfo_[eqIdx].isDirichlet = true;
            }
            else if constexpr (numTurbulenceEq(ModelTraits::turbulenceModel()) == 2)
            {
                if (eqIdx == Indices::turbulentKineticEnergyEqIdx || eqIdx == Indices::dissipationEqIdx)
                    BoundaryTypes<numEq>::boundaryInfo_[eqIdx].isDirichlet = true;
            }
        }
    }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        Wall boundary condition.
     */
    bool hasWall() const
    { return isWall_; }

    /*!
     * \brief Reset the boundary types for one equation.
     */
    void resetEq(const int eqIdx)
    {
        ParentType::resetEq(eqIdx);
        isWall_ = false;
    }

protected:
    bool isWall_;
};


/*!
 * \ingroup RANSModel
 * \brief Class to specify the type of a boundary condition for the RANS extension to the Navier-Stokes model.
 */
template <class ModelTraits, int numEq>
class RANSFCBoundaryTypes : public NavierStokesMomentumBoundaryTypes<numEq>
{
    using ParentType = NavierStokesMomentumBoundaryTypes<numEq>;
    static constexpr auto dimWorld = ModelTraits::dim();
    using Indices = typename ModelTraits::Indices;
    static_assert(dimWorld > 1, "Wall conditions cannot be set for 1D domains.");
public:
    RANSFCBoundaryTypes()
    {
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
            resetEq(eqIdx);
    }

    void setWall()
    {
        isWall_ = true;
        for (int eqIdx=0; eqIdx < numEq; ++eqIdx)
        {
            if constexpr (dimWorld == 3)
            {
                if ((eqIdx == Indices::velocityXIdx)
                 || (eqIdx == Indices::velocityYIdx)
                 || (eqIdx == Indices::velocityZIdx))
                    BoundaryTypes<numEq>::boundaryInfo_[eqIdx].isDirichlet = true;
            }
            else if constexpr (dimWorld == 2)
            {
                if ((eqIdx == Indices::velocityXIdx)
                 || (eqIdx == Indices::velocityYIdx))
                    BoundaryTypes<numEq>::boundaryInfo_[eqIdx].isDirichlet = true;
            }
            else
                DUNE_THROW(Dune::NotImplemented, "1D Turbulence models are not supported");
        }
    }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        Wall boundary condition.
     */
    bool hasWall() const
    { return isWall_; }

    /*!
     * \brief Reset the boundary types for one equation.
     */
    void resetEq(const int eqIdx)
    {
        ParentType::resetEq(eqIdx);
        isWall_ = false;
    }

protected:
    bool isWall_;
};

} // end namespace Dumux

#endif
