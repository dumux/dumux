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
 *
 * \brief The spatial parameters for the fully coupled tutorial problem
 *        which uses the twophase box model.
 */
#ifndef DUMUX_EXERCISE_THREE_SPATIAL_PARAMS_HH
#define DUMUX_EXERCISE_THREE_SPATIAL_PARAMS_HH

// include parent spatialparameters
#include <dumux/material/spatialparams/fv.hh>

// include material laws
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux {

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief The spatial parameters for the fully coupled tutorial problem
 *        which uses the twophase box model.
 */
template<class FVGridGeometry, class Scalar>
class ExerciseThreeSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar, ExerciseThreeSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = ExerciseThreeSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Dune::FieldMatrix<Scalar, dim, dim>;
    // get material law from property system
    using MaterialLaw = EffToAbsLaw<RegularizedBrooksCorey<Scalar>>;
    using MaterialLawParams = typename MaterialLaw::Params;

    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     */
    ExerciseThreeSpatialParams(std::shared_ptr<const FVGridGeometry>& fvGridGeometry)
    : ParentType(fvGridGeometry)
    , K_(0)
    , KLens_(0)
    {
        //set main diagonal entries of the permeability tensor to a value
        //setting to one value means: isotropic, homogeneous
        for (int i = 0; i < dim; i++)
        {
            K_[i][i] = 1e-7;
            KLens_[i][i] = 1e-10;
        }

        //set residual saturations
        materialParams_.setSwr(0.0);
        materialParamsLens_.setSwr(0.1);
        materialParams_.setSnr(0.0);
        materialParamsLens_.setSnr(0.1);

        //parameters of Brooks & Corey Law
        materialParams_.setPe(500.0);
        materialParamsLens_.setPe(1000.0);
        materialParams_.setLambda(2);
        materialParamsLens_.setLambda(2);
    }


    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const

    {
        if (isInLens(globalPos))
            return KLens_;
        return K_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens(globalPos))
            return 0.1;
        return 0.2;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param globalPos The global position
     *
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens(globalPos))
            return materialParamsLens_;
        return materialParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The position of the center of the element
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

    //! if we are in the lens
    bool isInLens(const GlobalPosition& globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return (x < 40 && x > 20 && y > 35 && y < 45) ||
               (x < 50 && x > 30 && y < 30 && y > 15);
    }

private:

    Dune::FieldMatrix<Scalar, dim, dim> K_;
    Dune::FieldMatrix<Scalar, dim, dim> KLens_;
    // Object that holds the values/parameters of the selected material law.
    MaterialLawParams materialParams_;
    MaterialLawParams materialParamsLens_;
};
} // end namespace Dumux
#endif
