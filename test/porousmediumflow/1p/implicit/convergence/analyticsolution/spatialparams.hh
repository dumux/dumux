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
 * \ingroup OnePTests
 * \brief The spatial params of the incompressible single-phase convergence test with analytic solution
 */

#ifndef DUMUX_CONVERGENCE_TEST_ONEP_SPATIALPARAMS_HH
#define DUMUX_CONVERGENCE_TEST_ONEP_SPATIALPARAMS_HH

#include <cmath>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dune/common/fmatrix.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial params of the incompressible single-phase convergence test with analytic solution
 */
template<class GridGeometry, class Scalar>
class ConvergenceTestSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             ConvergenceTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           ConvergenceTestSpatialParams<GridGeometry, Scalar>>;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dimWorld = GridView::dimensionworld;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;


public:
    // export permeability type
    using PermeabilityType = DimWorldMatrix;

    ConvergenceTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        c_ = getParam<Scalar>("Problem.C");
    }

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position where we evaluate
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        PermeabilityType K(0.0);

        using std::cos;
        using std::sin;
        using std::exp;

        const Scalar x = globalPos[0];
        K[0][0] = 1.0;
        K[0][1] = -c_/(2*omega_) * sin(omega_*x);
        K[1][0] = K[0][1];
        K[1][1] = exp(-2)*(1 + c_*cos(omega_*x));

        return K;
    }

    /*! \brief Defines the porosity in [-].
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

private:
    static constexpr Scalar omega_ = M_PI;
    Scalar permeability_;
    Scalar c_;
};

} // end namespace Dumux

#endif
