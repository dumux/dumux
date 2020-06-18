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
 * \ingroup ShallowWaterTests
 * \brief A test for the Shallow water model (bowl).
 */
#ifndef DUMUX_BOWL_TEST_PROBLEM_HH
#define DUMUX_BOWL_TEST_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief A wetting and drying test with sloshing water in a bowl.
 *
 * The domain is 4 meters long and 4 meters wide. The center of the domain is loacted at
 * x = 0 and y = 0. There is no flow over the boundaries and no friction is considered.
 *
 * This example is demanding for the implicit model if a high mesh resolution is applied
 * (e.g. 150x150 cells) in combination with a large time step size. Using the new limiting
 * (UpwindFluxLimiting = true) will help to improve the convergence for such cases.
 *
 * This test uses a low mesh resoultion and only ensures that UpwindFluxLimiting for the mobility
 * works. For low mesh resolution the solution is very diffusive so that the oscillation is dampened.
 * This gets better with grid refinement (not tested here).
 *
 * The results are checked against a analytical solution which is based on the "Thacker-Solution"
 * (William Thacker, "Some exact solutions to the nonlinear shallow-water wave equations", Journal
 * of Fluid Mechanics, 107:499–508, 1981, doi: https://doi.org/10.1017/S0022112081001882).
 * This implements the oscillating solution in a circular bowl (Section 4 in the paper).
 * Further examples and details on the solution are given
 * in SWASHES (Shallow Water Analytic Solutions for Hydraulic and Environmental Studies,
 * https://www.idpoisson.fr/swashes/).
 *
 * This problem uses the \ref ShallowWaterModel
 */
template <class TypeTag>
class BowlProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    using NeumannFluxes = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;


public:
    BowlProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {

        name_ = getParam<std::string>("Problem.Name");
        bowlDepthAtCenter_ =  getParam<Scalar>("Problem.BowlDepthAtCenter");
        bowlParaboloidRadius_ =  getParam<Scalar>("Problem.BowlParaboloidRadius");

        // Thacker (1981) Eq. (43)
        using std::sqrt;
        bowlAnalyticParameterOmega_ = sqrt(8.0 * getParam<Scalar>("Problem.Gravity") * bowlDepthAtCenter_) / bowlParaboloidRadius_;
        std::cout << "One full oscillation period of the water table is: "
                  << oscillationPeriodInSeconds() << " seconds." << std::endl;

        // Thacker (1981) Eq. (50)
        const auto D0PlusEta = bowlDepthAtCenter_ + getParam<Scalar>("Problem.BowlInitialWaterElevationAtCenter");
        const auto D0PlusEtaSquared = D0PlusEta*D0PlusEta;
        const auto D0Squared = bowlDepthAtCenter_*bowlDepthAtCenter_;
        bowlAnalyticParameterA_ = (D0PlusEtaSquared - D0Squared)/(D0PlusEtaSquared + D0Squared);

        // check constraint Thacker (1981) text after Eq. (49)
        if (bowlAnalyticParameterA_ >= 1.0)
            DUNE_THROW(Dune::InvalidStateException, "Parameter A has to be smaller than unity!");
    }

    //! One oscillation period of the water table (analytically this goes on forever)
    Scalar oscillationPeriodInSeconds() const
    { return 2*M_PI/bowlAnalyticParameterOmega_; }

    //! Get the analytical water depth at time t and position pos
    PrimaryVariables analyticalSolution(const Scalar t, const GlobalPosition& pos) const
    {
        using std::sqrt;
        using std::cos;
        using std::sin;

        // see Thacker (1981) Eq. (51) for formula
        const auto radiusSquared = pos[0]*pos[0] + pos[1]*pos[1];
        const auto LSquared = bowlParaboloidRadius_*bowlParaboloidRadius_;
        const auto A = bowlAnalyticParameterA_;
        const auto omega = bowlAnalyticParameterOmega_;
        const auto D0 = bowlDepthAtCenter_;

        const auto oneMinusASq = 1.0 - A*A;
        const auto oneMinusACosOmegaT = 1.0 - A*cos(omega*t);
        const auto ratioSq = oneMinusASq / (oneMinusACosOmegaT*oneMinusACosOmegaT);
        const auto localRadiusSq = radiusSquared / LSquared;
        // bowl depth function (cf. D in Thacker (1981))
        const auto D = D0*(1.0 - localRadiusSq);
        // height above equilibrium water level (cf. h in Thacker (1981))
        const auto h = D0*(sqrt(ratioSq) - 1.0 - localRadiusSq*(ratioSq - 1.0));
        // see remark about total water depth in Thacker (1981) beginning section 2
        const auto analyticalWaterDepth = h + D;

        const auto halfOmegaASinOmegaT = 0.5*omega*A*sin(omega*t);
        const auto analyticalVelocityX = pos[0]*halfOmegaASinOmegaT/oneMinusACosOmegaT;
        const auto analyticalVelocityY = pos[1]*halfOmegaASinOmegaT/oneMinusACosOmegaT;

        // The radius of the shoreline (where h + D = 0), Eq. (48)
        const auto h0 = D0*(sqrt(ratioSq) - 1.0); // h in the middle of the bowl (r=0)
        const auto shoreLineRadiusSquared = LSquared*(D0/(D0 + h0));

        // outside shoreline the water height and velocity is zero
        if (radiusSquared > shoreLineRadiusSquared)
            return { 0.0, 0.0, 0.0 };
        else
            return { analyticalWaterDepth, analyticalVelocityX, analyticalVelocityY };
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Specifies the neumann boundary
     *
     *  We need the Riemann invariants to compute the values depending of the boundary type.
     *  Since we use a weak imposition we do not have a dirichlet value. We impose fluxes
     *  based on q, h, etc. computed with the Riemann invariants
     */
    template<class ElementFluxVariablesCache>
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = this->spatialParams().gravity(scvf.center());
        std::array<Scalar, 3> boundaryStateVariables;

        //no flow with zero normal velocity and tangential velocity
        const auto vNormalGhost = -(nxy[0] * insideVolVars.velocity(0) +  nxy[1] * insideVolVars.velocity(1));
        const auto vTangentialGhost = -nxy[1] * insideVolVars.velocity(0) + nxy[0] * insideVolVars.velocity(1);

        boundaryStateVariables[0] = insideVolVars.waterDepth();
        boundaryStateVariables[1] =  nxy[0] * vNormalGhost - nxy[1] * vTangentialGhost;
        boundaryStateVariables[2] =  nxy[1] * vNormalGhost + nxy[0] * vTangentialGhost;

        const auto riemannFlux =
            ShallowWater::riemannProblem(insideVolVars.waterDepth(), boundaryStateVariables[0],
                                         insideVolVars.velocity(0), boundaryStateVariables[1],
                                         insideVolVars.velocity(1), boundaryStateVariables[2],
                                         insideVolVars.bedSurface(), insideVolVars.bedSurface(),
                                         gravity, nxy);

        values[Indices::massBalanceIdx] = riemannFlux[0];
        values[Indices::velocityXIdx]   = riemannFlux[1];
        values[Indices::velocityYIdx]   = riemannFlux[2];

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values at position globalPos
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        using std::max; // regularize so that we are virtually dry but not completely dry
        return { max(analyticalSolution(0, globalPos)[0], 1e-5), 0.0, 0.0 };
    };

    // \}

private:
    Scalar bowlDepthAtCenter_;
    Scalar bowlParaboloidRadius_;
    Scalar bowlAnalyticParameterOmega_;
    Scalar bowlAnalyticParameterA_;
    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
