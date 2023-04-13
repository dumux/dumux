// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial params of the incompressible single-phase convergence test with analytic solution
 */
template<class GridGeometry, class Scalar>
class ConvergenceTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         ConvergenceTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;

    using ThisType =  ConvergenceTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

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
