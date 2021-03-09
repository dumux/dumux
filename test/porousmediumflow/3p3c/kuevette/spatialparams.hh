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
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the kuevette problem.
 */

#ifndef DUMUX_KUEVETTE3P3CNI_SPATIAL_PARAMS_HH
#define DUMUX_KUEVETTE3P3CNI_SPATIAL_PARAMS_HH

#include <dune/common/float_cmp.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>

#include <dumux/material/fluidmatrixinteractions/3p/parkervangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/3p/napladsorption.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the kuevette problem
 */
template<class GridGeometry, class Scalar>
class KuevetteSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         KuevetteSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       KuevetteSpatialParams<GridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using ThreePhasePcKrSw = FluidMatrix::ParkerVanGenuchten3PDefault<Scalar>;
    using AdsorptionModel = FluidMatrix::ThreePNAPLAdsorption<Scalar>;

public:
    using PermeabilityType = Scalar;

    KuevetteSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurveFine_("SpatialParams.Fine")
    , pcKrSwCurveCoarse_("SpatialParams.Coarse")
    , adsorption_("SpatialParams")
    {
        // intrinsic permeabilities
        fineK_ = 6.28e-12;
        coarseK_ = 9.14e-10;

        // porosities
        finePorosity_ = 0.42;
        coarsePorosity_ = 0.42;
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return makeFluidMatrixInteraction(pcKrSwCurveFine_, adsorption_);
        else
            return makeFluidMatrixInteraction(pcKrSwCurveCoarse_, adsorption_);
    }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    {
        return ((Dune::FloatCmp::ge<Scalar>(globalPos[0], 0.13)
                    && Dune::FloatCmp::le<Scalar>(globalPos[0], 1.24)
                    && Dune::FloatCmp::ge<Scalar>(globalPos[1], 0.32)
                    && Dune::FloatCmp::le<Scalar>(globalPos[1], 0.60))
                || (Dune::FloatCmp::ge<Scalar>(globalPos[0], 1.20)
                    && Dune::FloatCmp::le<Scalar>(globalPos[1], 0.15)));
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    const ThreePhasePcKrSw pcKrSwCurveFine_;
    const ThreePhasePcKrSw pcKrSwCurveCoarse_;
    const AdsorptionModel adsorption_;

};

} // end namespace Dumux

#endif
