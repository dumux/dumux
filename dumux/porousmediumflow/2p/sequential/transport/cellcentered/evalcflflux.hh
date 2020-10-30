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
 * \ingroup SequentialTwoPModel
 * \brief Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition.
 */
#ifndef DUMUX_EVALCFLFLUX_HH
#define DUMUX_EVALCFLFLUX_HH

#include <dumux/porousmediumflow/2p/sequential/properties.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPModel
 * \brief Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition.
 *
 *  Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition of the form
 *
 *  \f[
 *  \frac{F_i \Delta t}{V_{p_i}} < 1
 *  \f]
 *
 *  where \f$ V_{p_i} \f$ is the pore volume of cell i and
 *
 *  \f[
 *  F_i = \sum f_{ij}
 *  \f]
 *
 * with \f$ f_{ij} \f$ being the CFL-flux over edge \f$ ij \f$.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class EvalCflFlux
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Intersection = typename GridView::Intersection;
    using Element = typename GridView::Traits::template Codim<0>::Entity;

public:

    //! For initialization
    void initialize()
    {}

    /*!
     * \brief adds a flux to the cfl-criterion evaluation
     *
     *  \param lambdaW        wetting phase mobility
     *  \param lambdaNw       nonwetting phase mobility
     *  \param viscosityW     wetting phase viscosity
     *  \param viscosityNw    nonwetting phase viscosity
     *  \param flux           flux to add
     *  \param intersection   intersection corresponding to the flux
     *  \param phaseIdx       index of the phase (wetting, nonwetting)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNw, Scalar& viscosityW, Scalar& viscosityNw, Scalar flux,
                 const Intersection& intersection, int phaseIdx)
    {}

    /*!
     * \brief adds a flux to the cfl-criterion evaluation
     *
     *  \param lambdaW        wetting phase mobility
     *  \param lambdaNw       nonwetting phase mobility
     *  \param viscosityW     wetting phase viscosity
     *  \param viscosityNw    nonwetting phase viscosity
     *  \param flux           flux to add
     *  \param element        element corresponding to the flux
     *  \param phaseIdx       index of the phase (wetting, nonwetting)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNw, Scalar& viscosityW, Scalar& viscosityNw, Scalar flux,
                 const Element& element, int phaseIdx)
    {}

    /*!
     * \brief Returns the CFL flux-function
     *
     *  \param element       element on which the CFL-criterion is evaluated
     *  \return fluxFunction for the calculation of the CFL time-step (\f$ 1/F_i \f$)
     */
    Scalar getCflFluxFunction(const Element& element)
    {
        return 0.0;
    }

    /*!
     * \brief  Returns the CFL time-step
     *
     *  \param element       element on which the CFL-criterion is evaluated
     *  \return CFL time-step
     */
    Scalar getDt(const Element& element)
    {
        return 0.0;
    }

    //! reset function
    void reset()
    {}

    /*! \brief Constructs a EvalCflFlux instance */
    EvalCflFlux ()
    {}
};
} // end namespace Dumux

#endif
