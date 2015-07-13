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
#ifndef DUMUX_IMPLICIT_GRIDADAPTINDICATORDEFAULT_HH
#define DUMUX_IMPLICIT_GRIDADAPTINDICATORDEFAULT_HH

#include <dumux/common/propertysystem.hh>

/**
 * @file
 * @brief  Class defining a default indicator for grid adaptation
 */
namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Problem);
}

/*!\ingroup ImplicitGridAdaptIndicator
 * @brief  Class defining a default indicator for grid adaptation
 *
 *Default implementation
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class ImplicitGridAdaptIndicatorDefault
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:
    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     */
    void calculateIndicator()
    {}

    /*! \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool refine(const Element& element)
    {
        return false;
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool coarsen(const Element& element)
    {
        return false;
    }

    /*! \brief Initializes the adaptation indicator class*/
    void init()
    {};

    /*! \brief Constructs a GridAdaptationIndicator for initialization of an adaptive grid
     *
     * Default implementation
     *
     * \param problem The problem object
     * \param adaptationIndicator Indicator whether a be adapted
     */
    ImplicitGridAdaptIndicatorDefault(Problem& problem)
    {}
};
}

#endif
