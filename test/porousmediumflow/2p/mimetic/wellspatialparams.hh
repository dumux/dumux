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
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model
 */
#ifndef DUMUX_MIMETIC_WELL_SPATIAL_PARAMS_HH
#define DUMUX_MIMETIC_WELL_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p/implicit/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class WellSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(WellSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(WellSpatialParams, SpatialParams, WellSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(WellSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<EffectiveLaw>;
};
}
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the WellProblem which uses the
 *        two-phase fully implicit model
 */
template<class TypeTag>
class WellSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

    //get the material law from the property system
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    // export permeability type
    using PermeabilityType = Tensor;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    WellSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        materialLawParams_.setSwr(0.0);
        materialLawParams_.setSnr(0.0);

        // parameters for the Brooks-Corey Law
        // entry pressures
        materialLawParams_.setPe(1.0e3);
        // Brooks-Corey shape parameters
        materialLawParams_.setLambda(2);

        double pi = 4.0*atan(1.0);
        double theta = 0.0;
        double cost = cos(theta);
        double sint = sin(theta);

        double k1 = 1.0e-10;
        double k2 = 1.0e-10;
        double k3 = 1.0e-11;

        permeability_[0][0] = k1;
        permeability_[1][1] = k2;
        permeability_[0][1] = permeability_[1][0] = 0.0;

        // set permeability values
        permeability_[2][2] = k3;
        permeability_[0][2] = permeability_[2][0] = permeability_[1][2] = permeability_[2][1] = 0.0;
    }

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return permeability_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialLawParams_;
    }



private:
    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    PermeabilityType permeability_;
    MaterialLawParams materialLawParams_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
