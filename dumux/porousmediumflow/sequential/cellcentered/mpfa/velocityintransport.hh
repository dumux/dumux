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
#ifndef DUMUX_FVMPFAVELOCITYINTRANSPORT_HH
#define DUMUX_FVMPFAVELOCITYINTRANSPORT_HH

/**
 * @file
 * @brief  Implementation of a interface velocity class for MPFA models.
 */

#include <dune/grid/common/gridenums.hh>
#include <dumux/porousmediumflow/sequential/properties.hh>
#include "properties.hh"

namespace Dumux
{
//! \ingroup IMPET
/*! \brief Implementation of a interface velocity class for MPFA models.
 *
 * Allows to calculate the MPFA-velocity in the transport model. (Less efficient!)
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FvMpfaVelocityInTransport
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using Intersection = typename GridView::Intersection;


public:
    //! Constructs a FvMpfaVelocityInTransport object
    /*!
     * \param problem A problem class object
     */
    FvMpfaVelocityInTransport(Problem& problem):
        problem_(problem)
    {
        calcVelocityInTransport_ = getParam<bool>("MPFA.CalcVelocityInTransport");
    }

    //! For initialization
    void initialize()
    {}

    //! Local velocity calculation
    void calculateVelocity(const Intersection& intersection, CellData& cellData)
    {
        problem_.pressureModel().calculateVelocity(intersection, cellData);
    }

    //! Local velocity calculation
    void calculateVelocityOnBoundary(const Intersection& intersection, CellData& cellData)
    {
        problem_.pressureModel().calculateVelocityOnBoundary(intersection, cellData);
    }

    //! \brief Indicates if velocity is reconstructed the transport step
    bool calculateVelocityInTransport()
    {
        return calcVelocityInTransport_;
    }

    /*! \brief Adds velocity output to the output file
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {}

private:
    Problem& problem_;
    bool calcVelocityInTransport_;
};
}
#endif
