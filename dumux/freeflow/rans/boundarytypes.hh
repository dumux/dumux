// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RANSModel
 * \copydoc Dumux::RANSBoundaryTypes
 */
#ifndef FREEFLOW_RANS_BOUNDARY_TYPES_HH
#define FREEFLOW_RANS_BOUNDARY_TYPES_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Class to specify the type of a boundary condition for the RANS extension to the Navier-Stokes model.
 */
template <class ModelTraits, int numEq>
class RANSBoundaryTypes : public NavierStokesBoundaryTypes<numEq>
{
    using ParentType = NavierStokesBoundaryTypes<numEq>;
    static constexpr auto dimWorld = ModelTraits::dim();
    using Indices = typename ModelTraits::Indices;
    static_assert(dimWorld > 1, "Wall conditions cannot be set for 1D domains.");
public:
    RANSBoundaryTypes()
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

} // end namespace Dumux

#endif
