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
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c mpnc comparison problem.
 */

#ifndef DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH
#define DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedlinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>

namespace Dumux {

/**
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c mpnc comparison problem.
 */
template<class FVGridGeometry, class Scalar>
class TwoPTwoCComparisonSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         TwoPTwoCComparisonSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       TwoPTwoCComparisonSpatialParams<FVGridGeometry, Scalar>>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum {dimWorld=GridView::dimensionworld};

    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;

public:
    //! Export permeability type
    using PermeabilityType = Scalar;
    //! Export the type used for the material law
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;


    TwoPTwoCComparisonSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry)
    {
        // intrinsic permeabilities
        coarseK_ = 1e-12;
        fineK_ = 1e-15;

        // the porosity
        porosity_ = 0.3;

         // residual saturations
        fineMaterialParams_.setSwr(0.2);
        fineMaterialParams_.setSnr(0.0);
        coarseMaterialParams_.setSwr(0.2);
        coarseMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setPe(1e4);
        coarseMaterialParams_.setPe(1e4);
        fineMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);
    }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (isFineMaterial_(scv.dofPosition()))
            return fineK_;
        else
            return coarseK_;
    }

    /*!
     * \brief Defines the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the sub-control volume.
     * \return the material parameters object
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param globalPos The global position of the sub-control volume.
     * \return The material parameters object
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position of the sub-control volume.
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:
    /*!
     * \brief Returns whether a given global position is in the
     *        fine-permeability region or not.
     * \param pos The global position of the sub-control volume.
     */
    static bool isFineMaterial_(const GlobalPosition &pos)
    {
        return
            30 - eps_ <= pos[0] && pos[0] <= 50 + eps_ &&
            20 - eps_ <= pos[1] && pos[1] <= 40 + eps_;
    }

    Scalar coarseK_;
    Scalar fineK_;
    Scalar porosity_;
    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;
    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
