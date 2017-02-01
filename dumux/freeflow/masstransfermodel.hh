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
  * \brief This files contains a class to adapt the mass transfer in the
  *        boundary layer of the free flow, near to a porous-interface.
  */
#ifndef DUMUX_MASS_TRANSFER_MODEL_HH
#define DUMUX_MASS_TRANSFER_MODEL_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{
/*!
  * \brief This files contains a class to adapt the mass transfer in the
  *        boundary layer of the free flow, near to a porous-interface.
  */
template <class TypeTag>
class MassTransferModel
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

/*!
  * \brief This files contains a class to adapt the mass transfer in the
  *        boundary layer of the free flow, near to a porous-interface.
  *
  * \param saturation Saturation of the component's phase inside the porous medium
  * \param porosity Porosity of the porous medium
  * \param blThickness Free flow boundary layer thickness
  * \param massTransferModel The index of the chosen mass transfer model
  */
public:
    MassTransferModel(Scalar saturation, Scalar porosity,
                      Scalar blThickness, unsigned int massTransferModel)
        : saturation_(saturation), porosity_(porosity),
          blThickness_(blThickness), massTransferModel_(massTransferModel)
    {
        moistureContent_ = saturation_ * porosity_;

        massTransferCoeff_ = 0.0; // dummy default
        charPoreRadius_ = 1e10; // dummy default
        capillaryPressure_ = 1e10; // dummy default
    }

    //! \brief Sets a mass transfer coefficient \f$[-]\f$.
    void setMassTransferCoeff(Scalar massTransferCoeff)
    { massTransferCoeff_ = massTransferCoeff; }

    //! \brief Sets the characteristic pore diameter \f$[m]\f$.
    void setCharPoreRadius(Scalar charPoreRadius)
    { charPoreRadius_ = charPoreRadius; }

    //! \brief Sets the capillary pressure \f$[Pa]\f$.
    void setCapillaryPressure(Scalar capillaryPressure)
    { capillaryPressure_ = capillaryPressure; }


    /*!
     * \brief Returns a mass transfer coefficient to weaken to diffusive mass transfer
     *        in vicinity of a porous surface.
     *
     * Provides different options for the determining the mass transfer coefficient:<br>
     * 0) no mass transfer model<br>
     * 1) exponential law with a saturation as base and a mass transfer coefficient as exponent<br>
     * 2) the conventional Schlünder model with one characteristic pore size<br>
     * 3) the Schlünder model with a variable characteristic pore size deduced from the
     *    two-phase relations (Van Genuchten curve)<br>
     * 4) a manually adapted Schlünder model
     */
    const Scalar massTransferCoefficient()
    {
        Scalar massTransferCoeff = 1.0;

        using std::sqrt;
        using std::pow;

        // no mass transfer model
        if (massTransferModel_ == 0)
        {
            massTransferCoeff = 1.0;
        }
        // exponential mass transfer law
        else if (massTransferModel_ == 1)
        {
            // use meaningful transfer coefficients
            assert (0 < massTransferCoeff_ && massTransferCoeff_ < 1.0);
            massTransferCoeff = pow(saturation_, massTransferCoeff_);
        }
        // Schlünder model (Schlünder, CES 1988)
        else if (massTransferModel_ == 2)
        {
            // check if characteristic pore radius was set
            assert (charPoreRadius_ < 9e9);
            massTransferCoeff = 1. + 2./M_PI * charPoreRadius_ / blThickness_
                                     * sqrt(M_PI/(4.*moistureContent_))
                                     * (sqrt(M_PI/(4.*moistureContent_)) - 1.);

            massTransferCoeff = 1./massTransferCoeff;
        }
        // Schlünder model (Schlünder, CES 1988) with variable char. pore diameter depending on Pc
        else if (massTransferModel_ == 3)
        {
            // check if capillary pressure was set
            assert (capillaryPressure_ < 9e9);
            Scalar surfaceTension = 72.85e-3; // source: Wikipedia
            Scalar charPoreRadius = 2 * surfaceTension / capillaryPressure_;

            massTransferCoeff = 1. + 2./M_PI * charPoreRadius / blThickness_
                                     * sqrt(M_PI/(4.*moistureContent_))
                                     * (sqrt(M_PI/(4.*moistureContent_)) - 1.);

            massTransferCoeff = 1./massTransferCoeff;
        }
        // modified Schlünder model
        else if (massTransferModel_ == 4)
        {
            // check if characteristic pore radius was set
            assert (charPoreRadius_ < 9e9);
            massTransferCoeff = 1. + 2./M_PI * charPoreRadius_ / blThickness_
                                     * (1./moistureContent_) * (1./moistureContent_ - 1.);

            massTransferCoeff = 1./massTransferCoeff;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "This mass transfer model is not implemented");

        // assert the massTransferCoeff is inside the validity range
        assert(!(massTransferCoeff > 1.0 || massTransferCoeff < 0.0));
        return massTransferCoeff;
    }

private:
    Scalar saturation_;
    Scalar porosity_;
    Scalar moistureContent_;
    Scalar blThickness_;
    unsigned int massTransferModel_;

    Scalar massTransferCoeff_;
    Scalar charPoreRadius_;
    Scalar capillaryPressure_;
};

} //end namespace

#endif // DUMUX_MASS_TRANSFER_MODEL_HH
