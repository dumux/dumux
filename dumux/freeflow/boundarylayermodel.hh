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
  * \brief This files contains a class to calculate the boundary layer thickness
  *        of a free flow based on analytical and empirical approaches.
  */
#ifndef DUMUX_BOUNDARY_LAYER_MODEL_HH
#define DUMUX_BOUNDARY_LAYER_MODEL_HH


#include <cassert>

#include <dumux/common/basicproperties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{

/*!
  * \brief This is a class to calculate the boundary layer thickness
  *        of a free flow based on analytical and empirical approaches.
  */
template <class TypeTag>
class BoundaryLayerModel
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

/*!
  * \brief This is a class to calculate the boundary layer thickness
  *        of a free flow based on analytical and empirical approaches.
  *
  * \param velocity Mean Free Flow Velocity
  * \param distance Distance to the current location (including artificial run-up lengths)
  * \param kinematicViscosity Kinematic viscosity of the free fluid
  * \param boundaryLayerModel The index of the chosen boundary layer model
  */
public:
    BoundaryLayerModel(const Scalar velocity, const Scalar distance,
                       const Scalar kinematicViscosity, unsigned int boundaryLayerModel)
        : velocity_(velocity), distance_(distance), kinematicViscosity_(kinematicViscosity),
          boundaryLayerModel_(boundaryLayerModel)
    {
        // distance has to be positive, that the boundary layer can not be zero
        assert (distance_ > 1e-10);

        constThickness_ = 0.0; // dummy default
        roughnessLength_ = 0.0; // dummy default
        yPlus_ = 5.0; // end of viscous sublayer
        hydraulicDiamater_ = 0.0; // dummy default
    }

    //! \brief Sets the constant boundary layer thickness \f$[m]\f$.
    void setConstThickness(Scalar constThickness)
    { constThickness_ = constThickness; }

    //! \brief Sets the roughness length \f$[m]\f$.
    void setRoughnessLength(Scalar roughnessLength)
    { roughnessLength_ = roughnessLength; }

    //! \brief Sets the \f$ y^+\f$-value \f$[-]\f$.
    void setYPlus(Scalar yPlus)
    { yPlus_ = yPlus; }

    //! \brief Sets the hydraulic diameter \f$[m]\f$.
    void setHydraulicDiameter(Scalar hydraulicDiamater)
    { hydraulicDiamater_ = hydraulicDiamater; }

    /*!
     * \brief Returns the momentum boundary layer thickness of the flow.
     *
     * Different models for the approximation of the boundary layer thickness are available:<br>
     * 0 = no boundary layer model<br>
     * 1 = constant boundary layer thickness, which can be set in the parameter file<br>
     * 2 = laminar: Blasius' analytical solution<br>
     * 3 = turbulent (smooth);<br>
     * 4 = viscous sublayer (smooth)<br>
     * 5 = viscous sublayer (rough)<br>
     * 6 = viscous sublayer (rough)<br>
     * 7 = viscous sublayer (smooth, aerodynamically smooth, and aerodynamically rough)<br>
     */
    const Scalar viscousBoundaryLayerThickness()
    {
        Scalar reynoldsX = velocity_ * distance_ / kinematicViscosity_;
        Scalar reynoldsD = velocity_ * hydraulicDiamater_ / kinematicViscosity_;

        // no boundary layer model
        if (boundaryLayerModel_ == 0)
        {
            // has to be empty for the multidomain models
        }
        // Constant Boundary Layer Thickness
        else if (boundaryLayerModel_ == 1)
        {
            // boundary layer thickness has to be positive
            assert (constThickness_ > 1e-10);
            return constThickness_;
        }
        // laminar: Blasius (analytical solution)
        else if (boundaryLayerModel_ == 2)
        {
            return 5.0 * distance_ / std::sqrt(reynoldsX);
        }
        // turbulent, smooth: turbulent boundary layer thickness
        // source: http://en.wikipedia.org/wiki/Boundary_layer_thickness
        else if (boundaryLayerModel_ == 3)
        {
            return 0.37 * distance_ / std::pow(reynoldsX, 0.2);
        }
        // turbulent, smooth: viscous sublayer thickness via friction coefficient (after Schultz-Grunow)
        // source: Pope, S. B. Turbulent flows Cambridge University Press, 2006, XXXIV, p. 307
        else if (boundaryLayerModel_ == 4)
        {
            Scalar cf = 0.37 * std::pow(std::log10(reynoldsX), -2.584);
            return yPlus_ * distance_ / (reynoldsX * std::sqrt(cf / 2.0));
        }
        // turbulent, fully rough: viscous sublayer thickness via friction coefficient
        // source: Truckenbrodt, E. Elementare Strömungsvorgänge dichteveränderlicher Fluide Fluidmechanik,
        //         Springer Berlin Heidelberg, 2008, doi: 10.1007/978-3-540-79024-2_1, p. 333
        else if (boundaryLayerModel_ == 5)
        {
            // roughness length has to be positive
            assert (roughnessLength_ > 1e-10);
            // application is bounded to specific roughness length
            assert (1e-6 < roughnessLength_ / distance_ && roughnessLength_ / distance_ < 1e-2);
            Scalar cf = std::pow(1.89 - 1.62 * std::log10(roughnessLength_ / distance_), -2.5);
            // application is bounded to rough cases, indicated by the line in the chart in Truckenbrodt
            // NOTE: disabling the assertion assumes that the cf of the hydrodynamically
            //       rough region is a good approximation of the cf in the hydrodynamically
            //       smooth case
            assert (130.0e-3 * std::pow(reynoldsX, -0.1872) < cf);
            return yPlus_ * distance_ / (reynoldsX * std::sqrt(cf / 2.0));
        }
        // turbulent, fully rough: viscous sublayer thickness via friction coefficient
        // source: Truckenbrodt, E. Elementare Strömungsvorgänge dichteveränderlicher Fluide Fluidmechanik,
        //         Springer Berlin Heidelberg, 2008, doi: 10.1007/978-3-540-79024-2_1, p. 333
        else if (boundaryLayerModel_ == 6)
        {
            // roughness length has to be positive
            assert (roughnessLength_ > 1e-10);
            Scalar cf = 0.024 * std::pow(roughnessLength_ / distance_, 1.0/6.0);
            return yPlus_ * distance_ / (reynoldsX * std::sqrt(cf / 2.0));
        }
        // turbulent, smooth, aerodynamically smooth and aerodynamically rough
        // estimation of the viscous sublayer thickness via friction factor
        // derived for flow in pipes or ducts
        // see: White, F. M. Fluid mechanics McGraw-Hill, 2011, 7th ed. (eq. 6.35)
        // source: Colebrook, C. F. and White, C. M. (1937). "Experiments with Fluid Friction in Roughened Pipes"
        //         Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences 161 (906): 367–381.
        // explicit form: Haaland, SE (1983). "Simple and Explicit Formulas for the Friction Factor in Turbulent Flow".
        //                Journal of Fluids Engineering (ASME) 105 (1): 89–90. doi:10.1115/1.3240948.
        else if (boundaryLayerModel_ == 7)
        {
            // roughness length has to be positive
            assert (roughnessLength_ > 1e-10);
            // hydraulic diameter has to be positive
            assert (hydraulicDiamater_ > 1e-10);
            Scalar sqrtF = 1.0 / (-1.8 * std::log10(std::pow(roughnessLength_ / hydraulicDiamater_ / 3.7, 1.11)
                                                    + 6.9 / reynoldsD));
            Scalar viscousSublayer = yPlus_ * distance_
                                     / (reynoldsX * std::sqrt(0.25 * sqrtF * sqrtF / 2.0));
            // assert we are not in the fully rough case
            assert (roughnessLength_ < viscousSublayer);
            return viscousSublayer;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "This boundary layer model is not implemented");

        return 0.0;
    }

    /*!
     * \brief Convenience function to return the mass boundary layer thickness
     *
     * The mass and momentum boundary layer thickness could different.
     */
    const Scalar massBoundaryLayerThickness()
    { return viscousBoundaryLayerThickness(); }

    /*!
     * \brief Convenience function to return the thermal boundary layer thickness
     *
     * The thermal and momentum boundary layer thickness could different.
     */
    const Scalar thermalBoundaryLayerThickness()
    { return viscousBoundaryLayerThickness(); }

private:
    Scalar velocity_;
    Scalar distance_;
    Scalar kinematicViscosity_;
    unsigned int boundaryLayerModel_;

    Scalar constThickness_;
    Scalar yPlus_;
    Scalar roughnessLength_;
    Scalar hydraulicDiamater_;
};

} //end namespace

#endif // DUMUX_BOUNDARY_LAYER_MODEL_HH
