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
 * \ingroup OnePTests
 * \brief The spatial params the incompressible test
 */
#ifndef DUMUX_MULTIDOMAIN_1P_2P_TEST_SPATIAL_PARAMS_HH
#define DUMUX_MULTIDOMAIN_1P_2P_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        incompressible 1p model
 */
template<class FVGridGeometry, class Scalar>
class TestSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar, TestSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = TestSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using MaterialLaw = EffToAbsLaw<VanGenuchten<Scalar>>;
    using MaterialLawParams = typename MaterialLaw::Params;

    using PermeabilityType = Scalar;

    TestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        materialParams_.setSwr(0.05);
        materialParams_.setSnr(0.0);
        materialParams_.setVgAlpha(0.0037);
        materialParams_.setVgn(4.7);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     * \return the intrinsic permeability
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return 1e-12;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 0.3;
    }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law.
     *
     * In this test, we use element-wise distributed material parameters.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

private:
    MaterialLawParams materialParams_;
};

} // end namespace Dumux

#endif
