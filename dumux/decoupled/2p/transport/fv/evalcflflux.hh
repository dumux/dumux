// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                 *
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
#ifndef DUMUX_EVALCFLFLUX_HH
#define DUMUX_EVALCFLFLUX_HH

/**
 * @file
 * @brief  Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition
 * @author Markus Wolff
 */
namespace Dumux
{
/*!\ingroup Saturation2p
 * @brief  Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition
 *
 *  Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition of the form
 *
 *  \f[\frac{F_i \Delta t}{V_{p_i}} < 1\f]
 *
 *  where \f$ V_{p_i} \f$ is the pore volume of cell i and
 *
 *  \f[F_i = \sum f_{ij}\f]
 *
 * with \f$f_{ij}\f$ being the CFL-flux over edge \f$ij\f$.
 *
 * Template parameters are:

 - TypeTag PropertyTag of the problem implementation
 */
template<class TypeTag>
class EvalCflFlux
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    //! adds a flux to the cfl-criterion evaluation
    /*!
     *  \param lambdaW        wetting phase mobility
     *  \param lambdaNW       non-wetting phase mobility
     *  \param viscosityW     wetting phase viscosity
     *  \param viscosityNW    non-wetting phase viscosity
     *  \param flux           flux to add
     *  \param intersection   intersection corresponding to the flux
     *  \param phaseIdx       index of the phase (wetting, non-wetting)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux, const Intersection& intersection, int phaseIdx)
    {}

    //! adds a flux to the cfl-criterion evaluation
    /*!
     *  \param lambdaW        wetting phase mobility
     *  \param lambdaNW       non-wetting phase mobility
     *  \param viscosityW     wetting phase viscosity
     *  \param viscosityNW    non-wetting phase viscosity
     *  \param flux           flux to add
     *  \param element        element corresponding to the flux
     *  \param phaseIdx       index of the phase (wetting, non-wetting)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux, const Element& element, int phaseIdx)
    {}

    //! adds a flux to the cfl-criterion evaluation
    /*!
     *  \param globalPos     global position
     *  \param element       element on which the CFL-criterion is evaluated
     *  \return fluxFunction for the calculation of the CFL-time step (\f$ 1/F_i\f$)
     */
    Scalar getCFLFluxFunction(const GlobalPosition& globalPos, const Element& element)
    {
        return 0.0;
    }

    /*! @brief Constructs a EvalCflFlux instance */
    EvalCflFlux ()
    {}
};
}

#endif
