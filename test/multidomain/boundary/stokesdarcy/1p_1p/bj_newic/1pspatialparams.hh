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
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the 1p cc model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePModel
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p cc model
 */
template<class GridGeometry, class Scalar>
class OnePSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar, OnePSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar, OnePSpatialParams<GridGeometry, Scalar>>;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    static constexpr int dim = GridView::dimension;
    using PermeabilityType = Dune::FieldMatrix<Scalar, dim, dim>;

    OnePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        K_ = PermeabilityType(0.0);
        M_ = PermeabilityType(0.0);

        K_[0][0] =           getParam<Scalar>("Darcy.SpatialParams.Permeability_xx");
        K_[0][1] = K_[1][0] = getParam<Scalar>("Darcy.SpatialParams.Permeability_xy");
        K_[1][1] =           getParam<Scalar>("Darcy.SpatialParams.Permeability_yy");

        porosity_ = getParam<Scalar>("Darcy.SpatialParams.Porosity");
        alphaBJ_ = getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");

        //for new ic
        epsInterface_ = getParam<Scalar>("Darcy.InterfaceParams.epsInterface");
        N_s_bl = getParam<Scalar>("Darcy.InterfaceParams.N_s_bl");
        N_1_bl = getParam<Scalar>("Darcy.InterfaceParams.N_1_bl");

        M_[0][0]= getParam<Scalar>("Darcy.InterfaceParams.M_1_bl");
        M_[1][1]= getParam<Scalar>("Darcy.InterfaceParams.M_2_bl");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param globalPos The global position
     * \return the intrinsic permeability
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return K_; }

    /*! \brief Define the porosity in [-].
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*! \brief Define the Beavers-Joseph coefficient in [-].
     *
     * \param globalPos The global position
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition& globalPos) const
    { return alphaBJ_; }

    Scalar epsInterfaceAtPos(const GlobalPosition& globalPos) const
    { return epsInterface_;}

    Scalar factorNMomentumAtPos(const GlobalPosition& globalPos) const
    { return N_s_bl;}

    Scalar factorNTangentialAtPos(const GlobalPosition& globalPos) const
    { return N_1_bl;}

    Dune::FieldMatrix<Scalar, dim, dim> matrixNTangentialAtPos(const GlobalPosition& globalPos) const
    { return M_; }



private:
    PermeabilityType K_;

    Scalar porosity_;
    Scalar alphaBJ_;

    //for new ic
    Scalar epsInterface_;
    Scalar N_s_bl;
    Scalar N_1_bl;
    Dune::FieldMatrix<Scalar, dim, dim> M_;

};

} // end namespace

#endif
