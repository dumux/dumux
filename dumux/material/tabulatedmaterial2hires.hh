/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A generic template for tabulated material laws that depend
 *        on two parameters.
 *
 * This class is used for rectengular areas of higher resolution and
 * is required to make sure that the material law is continuous.
 */
#ifndef DUMUX_TABULATED_MATERIAL2_HIRES_HH
#define DUMUX_TABULATED_MATERIAL2_HIRES_HH

#include "tabulatedmaterial2.hh"

namespace Dune {
    /*!
     * \brief A generic template for tabulated material laws that depend
     *        on two parameters.
     *
     * This class is used for rectengular areas of higher resolution and
     * is required to make sure that the material law is continuous.
     *
     * In order to achieve this, we use "transition zones" at the
     * edges of the high resolution area. If we need a value inside
     * such a transition zone, we interpolate between the value of the
     * low-resolution table and the one of the high-resolution table,
     * depending on how far the point is away from the edge.
     */
    template <class Traits>
    class TabulatedMaterial2HiRes : public TabulatedMaterial2<Traits>
    {
        typedef typename Traits::Scalar Scalar;
        typedef TabulatedMaterial2<Traits> ParentType;
        enum { numX = Traits::numX, numY = Traits::numY };

    public:
        TabulatedMaterial2HiRes()
            {
            };

        /*!
         * \brief Returns the total weight of the value of the hires
         *        area compared to the one of the lowres area.
         */
        Scalar hiresWeight(Scalar x, Scalar y) const
            {
                return northernWeight_(y)*southernWeight_(y)*
                       westernWeight_(x)*easternWeight_(x);
            };
        
    protected:
        /*!
         * \brief Returns the weighting factor of the high resolution
         *        area relative to the low resolution area on the
         *        northern border.
         */
        Scalar northernWeight_(Scalar y) const
            {
                assert(Traits::transitionNorth >= 0);

                // make sure that there is a transition zone
                if (Traits::transitionNorth == 0)
                    return 1.0;
                
                // we need the tWidth variable to avoid a bogus
                // division by zero warning in gcc
                Scalar tWidth = Traits::transitionNorth;

                return std::min(1.0, (y - Traits::yMin)/tWidth);
            }

        /*!
         * \brief Returns the weighting factor of the high resolution
         *        area relative to the low resolution area on the
         *        southern border.
         */
        Scalar southernWeight_(Scalar y) const
            {
                assert(Traits::transitionSouth >= 0);

                // make sure that there is a transition zone
                if (Traits::transitionSouth == 0)
                    return 1.0;

                // we need the tWidth variable to avoid a bogus
                // division by zero warning in gcc
                Scalar tWidth = Traits::transitionSouth;

                return std::min(1.0, (Traits::yMax - y)/tWidth);
            }

        /*!
         * \brief Returns the weighting factor of the high resolution
         *        area relative to the low resolution area on the
         *        western border.
         */
        Scalar westernWeight_(Scalar x) const
            {
                assert(Traits::transitionWest >= 0);

                // make sure that there is a transition zone
                if (Traits::transitionWest == 0)
                    return 1.0;

                // we need the tWidth variable to avoid a bogus
                // division by zero warning in gcc
                Scalar tWidth = Traits::transitionWest;
                return std::min(1.0, (x - Traits::xMin)/tWidth);
            }

        /*!
         * \brief Returns the weighting factor of the high resolution
         *        area relative to the low resolution area on the
         *        eastern border.
         */
        Scalar easternWeight_(Scalar x) const
            {
                assert(Traits::transitionEast >= 0);

                // make sure that there is a transition zone
                if (Traits::transitionEast == 0)
                    return 1.0;

                // we need the tWidth variable to avoid a bogus
                // division by zero warning in gcc
                Scalar tWidth = Traits::transitionEast;

                return std::min(1.0, (Traits::xMax - x)/tWidth);
            }
    };
}

#endif
