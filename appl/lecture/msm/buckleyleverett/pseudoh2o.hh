/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Properties of pure water \f$H_2O\f$.
 */
#ifndef DUMUX_PSEUDOH2O_HH
#define DUMUX_PSEUDOH2O_HH


#include <dumux/material/components/component.hh>
//#include "interface_BL.xml"

namespace Dumux
{
/*!
 * \brief Rough estimate for testing purposes of some oil.
 */
template <class ScalarT>
class PseudoH2O : public Component<ScalarT, PseudoH2O<ScalarT> >
{
    typedef Component<ScalarT, PseudoH2O<ScalarT> > ParentType;
    typedef ScalarT Scalar;
public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "H2O"; }

    /*!
     * \brief Rough estimate of the density of oil [kg/m^3].
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return density_;
    }

    /*!
     * \brief Rough estimate of the viscosity of oil kg/(ms).
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return viscosity_;
    };

    static void setViscosity(Scalar viscosity)
    {
        viscosity_ =  viscosity;
    }

    static void setDensity(Scalar density)
    {
        density_ =  density;
    }

private:
    static Scalar viscosity_;
    static Scalar density_;
};
template <class ScalarT>
typename PseudoH2O<ScalarT>::Scalar PseudoH2O<ScalarT>::viscosity_ = 0.01;
template <class ScalarT>
typename PseudoH2O<ScalarT>::Scalar PseudoH2O<ScalarT>::density_ = 1000;
} // end namepace

#endif
