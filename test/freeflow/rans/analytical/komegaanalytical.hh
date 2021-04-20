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
* \ingroup RANSTests
* \brief A analytic 1/2-D solution for the k-omega turbulence model.
*
* The analytic solution is given by
* \f[ v_\text{x} = 2 \cdot x^3 \f]
* \f[ v_\text{y} = 2 - 2 \cdot y^3 \f]
* \f[ p = 2 - 2 \cdot x y \f]
* \f[ k = 1 + x^2 y^2 \f]
* \f[ \omega = 2 - x^2 y^2 \f]
*/

#ifndef DUMUX_FREEFLOW_RANS_KOMEGA_ANALYTICAL_HH
#define DUMUX_FREEFLOW_RANS_KOMEGA_ANALYTICAL_HH

#include <math.h>
#include <dune/common/exceptions.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

template <class TypeTag>
class KOmegaAnalytical
{
public:

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    static constexpr auto dim = GridGeometry::GridView::dimensionworld;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ScalarGradientVector = Dune::FieldVector<Scalar, dim>;
    using ScalarGradientMatrix = Dune::FieldVector<ScalarGradientVector, dim>;
    using VelocityVector = Dune::FieldVector<Scalar, dim>;
    using VelocityGradientMatrix = Dune::FieldVector<VelocityVector, dim>;
    using VelocitySecondGradientMatrix = Dune::FieldVector<VelocityGradientMatrix, dim>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

public:

    /*!
    * \brief Returns the analytical solution of the problem at a given position.
    * \param globalPos The global position
    */
    PrimaryVariables analyticalSolutionAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = pressure_(globalPos);
        values[Indices::velocityXIdx] = velocity_(globalPos)[0];
        values[Indices::velocityYIdx] = velocity_(globalPos)[1];
        values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_(globalPos);
        values[Indices::dissipationIdx] = turbulentDissipation_(globalPos);
        return values;
    }

    PrimaryVariables analyticalSourceTermAtPos(const GlobalPosition &globalPos, const Scalar& rho, const Scalar& viscosity) const
    {
        PrimaryVariables source(0.0);
        source[Indices::pressureIdx] = pressureSource_(globalPos, rho, viscosity);
        source[Indices::velocityXIdx] = momentumSource_(globalPos, rho, viscosity)[0];
        source[Indices::velocityYIdx] = momentumSource_(globalPos, rho, viscosity)[1];
        source[Indices::turbulentKineticEnergyIdx] = tkeSource_(globalPos, rho, viscosity);
        source[Indices::dissipationIdx] = dissipationSource_(globalPos, rho, viscosity);
        return source;
    }

private:

    Scalar pressure_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        return 2.0 - 2.0 * x * y;
    }

    ScalarGradientVector pressureGradient_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        // index [derivative]
        ScalarGradientVector gradients(0.0);

        gradients[0] = -2.0 * y;
        gradients[1] = -2.0 * x;
        return gradients;
    }

    VelocityVector velocity_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        // index [velocity]
        VelocityVector v(0.0);

        v[0] = 2.0 * x * x * x;
        v[1] = 2.0 - 2.0 * y * y * y;
        return v;
    }

    //! \brief The gradients
    VelocityGradientMatrix velGradient_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        // index [velocity] [derivative]
        VelocityGradientMatrix gradients;

        gradients[0][0] = 6.0 * x * x;
        gradients[0][1] = 0.0;
        gradients[1][0] = 0.0;
        gradients[1][1] = -6.0 * y * y;
        return gradients;
    }

    //! \brief The second derivatives
    VelocitySecondGradientMatrix velGradient2_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        // index [velocity] [first derivative] [second derivative]
        VelocitySecondGradientMatrix gradients;

        gradients[0][0][0] = 12.0 * x;
        gradients[0][0][1] = 0.0;
        gradients[0][1][0] = gradients[0][0][1];
        gradients[0][1][1] = 0.0;

        gradients[1][0][0] = 0.0;
        gradients[1][0][1] = 0.0;
        gradients[1][1][0] = gradients[1][0][1];
        gradients[1][1][1] = -12.0 * y;

        return gradients;
    }

    VelocityGradientMatrix velocityProductGradient_(const GlobalPosition& globalPos) const
    {
        // index [velocity] [derivative]
        VelocityGradientMatrix products;
        // product rule d(v*v)
        // 2 vx dvxdx
        products[0][0] = 2.0 * velocity_(globalPos)[0] * velGradient_(globalPos)[0][0];
        // vx dvydy + vy dvxdy
        products[0][1] = (velocity_(globalPos)[0] * velGradient_(globalPos)[1][1]) + (velocity_(globalPos)[1] * velGradient_(globalPos)[0][1]);
        // vx dvydx + vy dvxdx
        products[1][0] = (velocity_(globalPos)[0] * velGradient_(globalPos)[1][0]) + (velocity_(globalPos)[1] * velGradient_(globalPos)[0][0]);
        // 2 vy dvydy
        products[1][1] = 2.0 * velocity_(globalPos)[1] * velGradient_(globalPos)[1][1];
        return products;
    }

    // shearStressTensorProduct
    Scalar shearStressTensorProduct_(const GlobalPosition& globalPos) const
    {
        Scalar shearStressTensorProduct = 0.0;

        ScalarGradientMatrix stressTensor;
        stressTensor[0][0] = 0.5 * (velGradient_(globalPos)[0][0] + velGradient_(globalPos)[0][0]);
        stressTensor[0][1] = 0.5 * (velGradient_(globalPos)[0][1] + velGradient_(globalPos)[1][0]);
        stressTensor[1][0] = 0.5 * (velGradient_(globalPos)[0][1] + velGradient_(globalPos)[1][0]);
        stressTensor[1][1] = 0.5 * (velGradient_(globalPos)[1][1] + velGradient_(globalPos)[1][1]);

        shearStressTensorProduct += stressTensor[0][0] * stressTensor[0][0];
        shearStressTensorProduct += stressTensor[0][1] * stressTensor[0][1];
        shearStressTensorProduct += stressTensor[1][0] * stressTensor[1][0];
        shearStressTensorProduct += stressTensor[1][1] * stressTensor[1][1];

        return shearStressTensorProduct;
    }

    Scalar turbulentKineticEnergy_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        return 1.0 + x * x * y * y;
    }

    ScalarGradientVector turbulentKineticEnergyGradient_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        // index [derivative]
        ScalarGradientVector gradients(0.0);

        gradients[0] = 2.0 * x * y * y;
        gradients[1] = 2.0 * y * x * x;
        return gradients;
    }

    ScalarGradientMatrix turbulentKineticEnergyGradient2_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        // index [first derivative] [second derivative]
        ScalarGradientMatrix gradients;

        gradients[0][0] = 2.0 * y * y;
        gradients[0][1] = 4.0 * x * y;
        gradients[1][0] = gradients[0][1];
        gradients[1][1] = 2.0 * x * x;
        return gradients;
    }

    Scalar turbulentDissipation_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        return 2.0 - x * x * y * y;
    }

    ScalarGradientVector turbulentDissipationGradient_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        // index [derivative]
        ScalarGradientVector gradients(0.0);

        gradients[0] = -2.0 * x * y * y;
        gradients[1] = -2.0 * x * x * y;
        return gradients;
    }

    ScalarGradientVector turbulentDissipationGradient2_(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        // index [first derivative] [second derivative]
        ScalarGradientVector gradients(0.0);

        gradients[0][0] = -2.0 * y * y;
        gradients[0][1] = -4.0 * x * y;
        gradients[1][0] = gradients[0][1];
        gradients[1][1] = -2.0 * x * x;
        return gradients;
    }

    // eddy viscosity and its gradient
    Scalar nut_(const GlobalPosition& globalPos) const
    {
        static const auto enableKOmegaDissipationLimiter = getParam<bool>("KOmega.EnableDissipationLimiter", true);
        if (enableKOmegaDissipationLimiter)
            DUNE_THROW(Dune::InvalidStateException, "This limiter is not available for the analytical solution");

        return turbulentKineticEnergy_(globalPos) / turbulentDissipation_(globalPos);
    }

    ScalarGradientVector nutGradient_(const GlobalPosition& globalPos) const
    {
        static const auto enableKOmegaDissipationLimiter = getParam<bool>("KOmega.EnableDissipationLimiter", true);
        if (enableKOmegaDissipationLimiter)
            DUNE_THROW(Dune::InvalidStateException, "This limiter is not available for the analytical solution");

        ScalarGradientVector gradients(0.0);
        gradients[0] = ( (turbulentKineticEnergyGradient_(globalPos)[0] * turbulentDissipation_(globalPos))
                       - (turbulentKineticEnergy_(globalPos) * turbulentDissipationGradient_(globalPos)[0]))
                     / (turbulentDissipation_(globalPos) * turbulentDissipation_(globalPos));
        gradients[1] = ( (turbulentKineticEnergyGradient_(globalPos)[1] * turbulentDissipation_(globalPos))
                       - (turbulentKineticEnergy_(globalPos) * turbulentDissipationGradient_(globalPos)[1]))
                     / (turbulentDissipation_(globalPos) * turbulentDissipation_(globalPos));
        return gradients;
    }

    Scalar turbulentTermGradientProduct_(const GlobalPosition& globalPos) const
    {
        Scalar turbTermProduct = 0.0;
        turbTermProduct += turbulentKineticEnergyGradient_(globalPos)[0] * turbulentDissipationGradient_(globalPos)[0];
        turbTermProduct += turbulentKineticEnergyGradient_(globalPos)[1] * turbulentDissipationGradient_(globalPos)[1];
        return turbTermProduct;
    }

    //////////////////////////////////
    /// Constants and Coefficients ///
    //////////////////////////////////

    Scalar sigmaK_() const
    { return 0.6; }

    Scalar sigmaOmega_() const
    { return 0.5; }

    Scalar betaK_() const
    { return 0.09; }

    Scalar betaOmega_() const
    { return 0.0708; }

    Scalar alpha_() const
    { return 0.52; }

    //! \brief Returns the \$f \sigma_{d} \$f constant
    Scalar sigmaD_(const Scalar& product) const
    { return nabProductBool_(product) ? 0.125 : 0.0; }

    //! \brief Returns if the product of \$f  \nabla K \$f and \$f  \nabla omega \$f is more than 0
    bool nabProductBool_(const Scalar& product) const
    { return product > 0.0; }

    /////////////////////////////
    /// Compiled Source Terms ///
    /////////////////////////////

    Scalar pressureSource_(const GlobalPosition& globalPos, const Scalar& density, const Scalar& viscosity) const
    { return density * (velGradient_(globalPos)[0][0] + velGradient_(globalPos)[1][1]); }

    VelocityVector momentumSource_(const GlobalPosition& globalPos, const Scalar& density, const Scalar& viscosity) const
    {
        VelocityVector momentum(0.0);

        // inertial term
        momentum[0] += density * velocityProductGradient_(globalPos)[0][0];
        momentum[0] += density * velocityProductGradient_(globalPos)[0][1];
        momentum[1] += density * velocityProductGradient_(globalPos)[1][0];
        momentum[1] += density * velocityProductGradient_(globalPos)[1][1];

        // pressure term
        momentum[0] += pressureGradient_(globalPos)[0];
        momentum[1] += pressureGradient_(globalPos)[1];

        // viscous term (molecular)
        momentum[0] -= 2.0 * viscosity * velGradient2_(globalPos)[0][0][0];
        momentum[0] -=       viscosity * velGradient2_(globalPos)[0][1][1];
        momentum[0] -=       viscosity * velGradient2_(globalPos)[1][0][1];

        momentum[1] -= 2.0 * viscosity * velGradient2_(globalPos)[1][1][1];
        momentum[1] -=       viscosity * velGradient2_(globalPos)[1][0][0];
        momentum[1] -=       viscosity * velGradient2_(globalPos)[0][1][0];

        // viscous term (turbulent)
        momentum[0] -= 2 * ((nutGradient_(globalPos)[0] * velGradient_(globalPos)[0][0]) + (nut_(globalPos) * velGradient2_(globalPos)[0][0][0]));
        momentum[0] -=     ((nutGradient_(globalPos)[1] * velGradient_(globalPos)[0][1]) + (nut_(globalPos) * velGradient2_(globalPos)[0][1][1]));
        momentum[0] -=     ((nutGradient_(globalPos)[1] * velGradient_(globalPos)[1][0]) + (nut_(globalPos) * velGradient2_(globalPos)[1][0][1]));

        momentum[1] -= 2 * ((nutGradient_(globalPos)[1] * velGradient_(globalPos)[1][1]) + (nut_(globalPos) * velGradient2_(globalPos)[1][1][1]));
        momentum[1] -=     ((nutGradient_(globalPos)[0] * velGradient_(globalPos)[1][0]) + (nut_(globalPos) * velGradient2_(globalPos)[1][0][0]));
        momentum[1] -=     ((nutGradient_(globalPos)[0] * velGradient_(globalPos)[1][1]) + (nut_(globalPos) * velGradient2_(globalPos)[0][1][0]));

        momentum[0] -= (2.0 / 2.0) * density * turbulentKineticEnergy_(globalPos);
        momentum[1] -= (2.0 / 2.0) * density * turbulentKineticEnergy_(globalPos);

        return momentum;
    }

    Scalar tkeSource_(const GlobalPosition& globalPos, const Scalar& density, const Scalar& viscosity) const
    {
        using std::min;
        Scalar tkeSource = 0.0;

        // inertial term
        tkeSource += (velGradient_(globalPos)[0][0] * turbulentKineticEnergy_(globalPos) + velocity_(globalPos)[0] * turbulentKineticEnergyGradient_(globalPos)[0]);
        tkeSource += (velGradient_(globalPos)[1][1] * turbulentKineticEnergy_(globalPos) + velocity_(globalPos)[1] * turbulentKineticEnergyGradient_(globalPos)[1]);

        // viscous term (molecular + turbulent)
        tkeSource -= ( ((viscosity + sigmaK_() * nut_(globalPos)) * turbulentKineticEnergyGradient2_(globalPos)[0][0])
                     + ((sigmaK_() * nutGradient_(globalPos)[0]) * turbulentKineticEnergyGradient_(globalPos)[0]) );
        tkeSource -= ( ((viscosity + sigmaK_() * nut_(globalPos)) * turbulentKineticEnergyGradient2_(globalPos)[1][1])
                     + ((sigmaK_() * nutGradient_(globalPos)[1]) * turbulentKineticEnergyGradient_(globalPos)[1]) );

        // Production term (plus 2006 K-limiter via P term)
        static const auto enableKOmegaProductionLimiter = getParam<bool>( "KOmega.EnableProductionLimiter", false);
        Scalar productionTerm = 2.0 * nut_(globalPos) * shearStressTensorProduct_(globalPos);
        if (enableKOmegaProductionLimiter)
        {
            Scalar productionAlternative = 20.0 * density * betaK_() * turbulentKineticEnergy_(globalPos) * turbulentDissipation_(globalPos);
            productionTerm = min(productionTerm, productionAlternative);
        }
        tkeSource -= productionTerm;

        // Destruction term
        tkeSource += betaK_() * density * turbulentKineticEnergy_(globalPos) * turbulentDissipation_(globalPos);

        return tkeSource;
    }

    Scalar dissipationSource_(const GlobalPosition& globalPos, const Scalar& density, const Scalar& viscosity) const
    {
        using std::min;
        Scalar dissipationSource = 0.0;

        // inertial term
        dissipationSource += (velGradient_(globalPos)[0][0] * turbulentDissipation_(globalPos) + velocity_(globalPos)[0] * turbulentDissipationGradient_(globalPos)[0]);
        dissipationSource += (velGradient_(globalPos)[1][1] * turbulentDissipation_(globalPos) + velocity_(globalPos)[1] * turbulentDissipationGradient_(globalPos)[1]);

        // viscous term (molecular + turbulent)
        dissipationSource -= ( ((viscosity + sigmaK_() * nut_(globalPos)) * turbulentKineticEnergyGradient2_(globalPos)[0][0])
                             + ((sigmaK_() * nutGradient_(globalPos)[0]) * turbulentKineticEnergyGradient_(globalPos)[0]) );
        dissipationSource -= ( ((viscosity + sigmaK_() * nut_(globalPos)) * turbulentKineticEnergyGradient2_(globalPos)[1][1])
                             + ((sigmaK_() * nutGradient_(globalPos)[1]) * turbulentKineticEnergyGradient_(globalPos)[1]) );

        // Production term (plus 2006 K-limiter via P term)
        static const auto enableKOmegaProductionLimiter = getParam<bool>("KOmega.EnableProductionLimiter", false);
        Scalar productionTerm = 2.0 * nut_(globalPos) * shearStressTensorProduct_(globalPos);
        if (enableKOmegaProductionLimiter)
        {
            Scalar productionAlternative = 20.0 * density * betaK_() * turbulentKineticEnergy_(globalPos) * turbulentDissipation_(globalPos);
            productionTerm = min(productionTerm, productionAlternative);
        }
        dissipationSource -= alpha_() * productionTerm * (turbulentDissipation_(globalPos) / turbulentKineticEnergy_(globalPos));

        // Destruction term
        dissipationSource += betaOmega_() * density * turbulentDissipation_(globalPos) * turbulentDissipation_(globalPos);

        // CrossDiffusion term
        dissipationSource -= sigmaD_(turbulentTermGradientProduct_(globalPos)) * density * turbulentTermGradientProduct_(globalPos) / turbulentDissipation_(globalPos);

        return dissipationSource;
    }

};
} // end namespace Dumux
#endif
