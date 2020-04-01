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
#ifndef DUMUX_GRIDADAPTINITIALIZATIONINDICATORDEFAULT_HH
#define DUMUX_GRIDADAPTINITIALIZATIONINDICATORDEFAULT_HH

#include "properties.hh"
#include <dumux/common/properties.hh>

#include <dune/common/dynvector.hh>

/**
 * @file
 * @brief  Class defining a start indicator for grid adaption
 */
namespace Dumux
{
/*!\ingroup IMPES
 * @brief  Class defining a start indicator for grid adaption
 *
 *Default implementation
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridAdaptInitializationIndicatorDefault
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using AdaptionIndicator = GetPropType<TypeTag, Properties::AdaptionIndicator>;

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

    int maxLevel()
    {
        return maxLevel_;
    }


    /*! \brief Initializes the adaption indicator class*/
    void init()
    {}

    /*! \brief Returns true if the IMPET-Model needs to be initialized*/
    bool initializeModel()
    {
        return false;
    }

    /*! \brief Constructs a GridAdaptionIndicator for initialization of an adaptive grid
     *
     * Default implementation
     *
     * \param problem The problem object
     * \param adaptionIndicator Indicator whether a be adapted
     */
    GridAdaptInitializationIndicatorDefault(Problem& problem, AdaptionIndicator& adaptionIndicator)
    {
        maxLevel_ = getParam<int>("GridAdapt.MaxLevel");
    }

private:
    int maxLevel_; //!< maximum allowed level
};
}

#endif
