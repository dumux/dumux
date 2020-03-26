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
 * \ingroup SequentialTwoPTwoCModel
 * \brief Class including the variables and data of discretized data of the constitutive relations.
 */
#ifndef DUMUX_FLUXDATA2P2C_HH
#define DUMUX_FLUXDATA2P2C_HH

#include <dumux/porousmediumflow/sequential/properties.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Class including the variables and data of discretized data of the constitutive relations.
 *
 * The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the sequential two-phase flow is stored,
 * as well as discretized data of constitutive relationships like mobilities, fractional flow functions
 * and capillary pressure. Thus, they have to be calculated just once in every time step or every iteration step.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FluxData2P2C
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        dim = GridView::dimension
    };

    enum
    {
        numEquations = getPropValue<TypeTag, Properties::NumEq>()
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    enum
    {
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };

    typename Dune::BlockVector<typename Dune::FieldVector<bool, numEquations>> isUpwindCell_;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using VelocityVector = Dune::FieldVector<DimVector, 2*dim>;
    VelocityVector velocity_[numPhases];

public:

    //! Constructor
    FluxData2P2C()
    {
        isUpwindCell_.resize(2*dim);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int fIdx = 0; fIdx < 2*dim; ++fIdx)
                velocity_[phaseIdx][fIdx] = 0.0;
    }

    /*!
     * \brief Returns the phase velocity vector at a cell-cell interface
     *
     * \param phaseIdx Index of a fluid phase
     * \param indexInInside Index of the cell-cell interface in this cell
     */
    const DimVector& velocity(int phaseIdx, int indexInInside)
    {
        return velocity_[phaseIdx][indexInInside];
    }

    /*!
     * \brief Returns the phase velocity vector at a cell-cell interface
     *
     * \param phaseIdx Index of a fluid phase
     * \param indexInInside Index of the cell-cell interface in this cell
     */
    const DimVector& velocity(int phaseIdx, int indexInInside) const
    {
        return velocity_[phaseIdx][indexInInside];
    }

    /*!
     * \brief Sets the phase velocity vector at a cell-cell interface
     *
     * \param phaseIdx Index of a fluid phase
     * \param indexInInside Index of the cell-cell interface in this cell
     * \param velocity Phase velocity vector which is stored
     */
    void setVelocity(int phaseIdx, int indexInInside, const DimVector& velocity)
    {
        velocity_[phaseIdx][indexInInside] = velocity;
    }

    /*!
     * \brief Returns the total velocity vector at a cell-cell interface
     *
     * \param indexInInside Index of the cell-cell interface in this cell
     */
    DimVector velocityTotal(int indexInInside)
    {
        return velocity_[wPhaseIdx][indexInInside]
                + velocity_[nPhaseIdx][indexInInside];
    }

    /*!
     * \brief Returns the total velocity vector at a cell-cell interface
     *
     * \param indexInInside Index of the cell-cell interface in this cell
     */
    DimVector velocityTotal(int indexInInside) const
    {
        return velocity_[wPhaseIdx][indexInInside]
                + velocity_[nPhaseIdx][indexInInside];
    }

    //! resizes the upwind vector for the case of hanging nodes
    void resize(int size)
    {
        isUpwindCell_.resize(size);
    }
    //! returns the size of the upwind vector which equals number of faces
    int size()
    {
        return isUpwindCell_.size();
    }

    /*!
     * \brief Functions returning upwind information
     *
     * \param indexInInside The local inside index of the intersection
     * \param equationIdx The equation index
     */
    const bool& isUpwindCell(int indexInInside, int equationIdx) const
    {
        return isUpwindCell_[indexInInside][equationIdx];
    }
    /*!
     * \brief Sets the upwind information
     *
     * \param indexInInside The local inside index of the intersection
     * \param equationIdx The equation index
     * \param value set true or false
     */
    void setUpwindCell(int indexInInside, int equationIdx, bool value)
    {
        isUpwindCell_[indexInInside][equationIdx] = value;
    }

    //! Console output for the FluxData
    void outputFluxData()
    {
        for(int banana=0; banana<isUpwindCell_.size(); banana++)
            printvector(std::cout, isUpwindCell_, "upwindInformation", "row", 3);
    }


};
} // end namespace Dumux
#endif
