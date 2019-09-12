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

#include <dune/geometry/quadraturerules.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p cc model
 */
template<class FVGridGeometry, class Scalar>
class OnePConvSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                             OnePConvSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                           OnePConvSpatialParams<FVGridGeometry, Scalar>>;

    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dimWorld = GridView::dimensionworld;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = DimWorldMatrix;

    OnePConvSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        omega_ = getParam<Scalar>("Problem.FreqFactor")*M_PI;
        c_ = getParam<Scalar>("Problem.PermFactor");
        alphaBJ_ = getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        PermeabilityType perm(0.0);

        using std::cos;
        using std::sin;
        using std::exp;

        // It seems that evaluating at the cell center, (equivalent to order 1 integration)
        // results in better rates but slightly larger errors.
        const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(element.geometry().type(), 1);
        for (auto&& qp : quad)
        {
            auto globalPos = element.geometry().global(qp.position());
            Scalar x = globalPos[0];

            auto integrationElement = element.geometry().integrationElement(qp.position());

            perm[1][1] += exp(-2)*(1 + c_*cos(omega_*x)) * qp.weight()*integrationElement;
            perm[0][1] += -c_/(2*omega_) * sin(omega_*x) * qp.weight()*integrationElement;
        }
        perm[1][0] = perm[0][1];
        perm /= scv.volume();

        perm[0][0] = 1.0;

        return perm;
    }

    /*! \brief Define the porosity in [-].
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    /*! \brief Define the Beavers-Joseph coefficient in [-].
     *
     * \param globalPos The global position
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition& globalPos) const
    { return alphaBJ_; }


private:
    Scalar alphaBJ_;
    Scalar c_;
    Scalar omega_;
};

} // end namespace

#endif
