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
 * \brief A test for the Shallow water model (Poiseuille flow).
 */
#ifndef DUMUX_POISEUILLE_FLOW_TEST_PROBLEM_HH
#define DUMUX_POISEUILLE_FLOW_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include "spatialparams.hh"
#include <dumux/common/parameters.hh>

#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>

namespace Dumux {

template <class TypeTag>
class PoiseuilleFlowProblem;

// Specify the properties for the problem
namespace Properties {

// Create new type tags
namespace TTag {
struct PoiseuilleFlow { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::PoiseuilleFlow>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PoiseuilleFlow>
{ using type = Dumux::PoiseuilleFlowProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoiseuilleFlow>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;

public:
    using type = PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PoiseuilleFlow>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PoiseuilleFlow>
{ static constexpr bool value = false; };
} // end namespace Properties

/*!
 * \ingroup ShallowWaterTests
 * \brief A simple test for the 2D flow in a channel with rough side walls (Poiseuille flow).
 *
 * The domain has a length L =  400 meters long and a width W = 100 meters.
 * The domain extent is from (x,y) = (0.0, -50.0) to (400.0, 50.0), i.e. the channel centreline is at y=0.
 * The bed level is sloped from z = -9.98 (x=0) to z = -10.0 meters (x=L).
 * The initital water depth corresponds to the analytical solution:
 * having a slope equal to ib = dTheta / L, where dTheta is the water level difference between upstream and downstream: dTheta = 0.02 m (positive downwards).
 * With L = 400 m, the slope is ib = 0.00005 m/m, resulting in a constant water depth along the channel of H = 10.0 m.

 * At the west/left  (inflow)  boundary a discharge is prescribed of Q_in = -408.75 m^3/s or q_in 4.0875 m^2/s per meter.
 * At the east/right (outflow) boundary a fixed water level is prescribed of theta = 0.0 m.
 * The south and north boundaries are set to roughwall type boundaries,
 * with a coefficient alphaWall = 1.0, where:
 *      alphaWall = 0.0: full    slip (smooth wall)
 * 0.0 <alphaWall < 1.0: partial slip (partially-rough wall)
 *      alphaWall = 1.0: no      slip (fully-rough wall)
 * Additionally these (south and north) boundaries are (automatically) set to no-flow boundaries.

 * The flow in the channel experiences two forces:
 * 1) the pressure gradient, driving the flow downstream
 * 2) the wall roughness,

 * It can be verified that this force balance reduces the momentum equation in X-direction to:
 *
 * \f[
 * \frac{\partial p}{\partial x} = nu_T \frac{\partial^2 u}{\partial y^2}
 * \f]
 * where nu_T is the turbulent viscosity.
 *
 * This ordinary differential equation can be solved (applying the boundary conditions to obtain the integration constants),
 * resulting in a parabolic velocity profile in lateral direction (over the width of the channel):
 *
 * \f[
 * u(y) = \frac{g*dTheta}{8*nu_T*L} * \left(4*y^2 - W^2\right)
 * \f]
 *
 * where y is the coordinate in lateral direction.
 * The velocity has a maximum value at y = 0:
 * \f[
 * u_{max} = \frac{g*dTheta*W^2}{8*nu_T*L}
 * \f]
 *
 * Therefore u_{max} can be calculated to be:
 *
 * \f[
 * u_{max} = \frac{9.81*0.02*100^2}{8*0.1*400} = 6.13125 m/s
 * \f]
 *
 * The formula for u(y) is also used to calculate the analytic solution.
 * It should be noted that u /= f(x,t) and the v = 0 in the whole domain.
 * Therefore momentumn advection/convection should be zero for this test, as u*du/dx = v*du/dy = u*dv/dx = v*dv/dy = 0
 *
 * This problem uses the \ref ShallowWaterModel
 */
template <class TypeTag>
class PoiseuilleFlowProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NeumannFluxes = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    PoiseuilleFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityY_.resize(gridGeometry->numDofs(), 0.0);
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        discharge_ = getParam<Scalar>("Problem.Discharge");
        hBoundary_ = getParam<Scalar>("Problem.HBoundary");
        turbViscosity_ = getParam<Scalar>("Problem.TurbViscosity");
        alphaWall_ = getParam<Scalar>("Problem.AlphaWall");
        ksWall_ = getParam<Scalar>("Problem.KsWall");
        wallFrictionLawType_ = getParam<std::string>("Problem.WallFrictionLaw");
        // Make the wallFrictionLawType_ lower case
        transform(wallFrictionLawType_.begin(), wallFrictionLawType_.end(), wallFrictionLawType_.begin(), ::tolower);
    }

    //! Get the analytical water depth
    const std::vector<Scalar>& getExactWaterDepth()
    {
        return exactWaterDepth_;
    }

    //! Get the analytical U-velocity
    const std::vector<Scalar>& getExactVelocityX()
    {
        return exactVelocityX_;
    }

    //! Get the analytical V-velocity
    const std::vector<Scalar>& getExactVelocityY()
    {
        return exactVelocityY_;
    }

    //! Update the analytical solution
    void updateAnalyticalSolution()
    {
        using std::pow;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto eIdx = this->gridGeometry().elementMapper().index(element);
            GlobalPosition globalPos = element.geometry().center();
            auto gravity = this->spatialParams().gravity(globalPos);
            Scalar y = globalPos[1];
            Scalar W = 100.0; // for now. Should be read from the grid input
            Scalar h = hBoundary_; // hBoundary_- bedSurface;
            Scalar u = -(gravity*bedSlope_/(8.0*turbViscosity_)) * (4.0*pow(y,2.0) - pow(W,2.0));
            exactWaterDepth_[eIdx] = h;
            exactVelocityX_[eIdx] = u;
            exactVelocityY_[eIdx] = 0.0;
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the boundary type is set
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
     *
     * \param element
     * \param fvGeometry
     * \param elemVolVars
     * \param scvf
     */
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

        // impose discharge at the left side
        if (scvf.center()[0] < 0.0 + eps_)
        {
            // Prescribe the exact q-distribution long the inflow boundaryStateVariables
            // based on the parabolic u-profile:
            // q(y) = (g*H*ib/(8*nu_T)) * (4*y^2^-W^2)
            // Note the opposite sign to the velocity u, due to the fact that it is an inflow discharge
            auto gravity = this->spatialParams().gravity(scvf.center());
            Scalar y = scvf.center()[1];
            Scalar W = 100.0; // for now. Should be read from the grid input
            Scalar h = hBoundary_; // hBoundary_- bedSurface;
            // Now compute a weighted average between the constant q (from input) and the parabolic q from the analytical solution
            // based on the prescribed alphaWall
            Scalar q_in0 = (gravity*h*bedSlope_/(8.0*turbViscosity_)) * (4.0*pow(y,2.0) - pow(W,2.0));
            Scalar q_in1 = discharge_;
            Scalar q_in  = (1.0-alphaWall_)*q_in1 + alphaWall_*q_in0;
            boundaryStateVariables = ShallowWater::fixedDischargeBoundary(q_in,
                                                                          insideVolVars.waterDepth(),
                                                                          insideVolVars.velocity(0),
                                                                          insideVolVars.velocity(1),
                                                                          gravity,
                                                                          nxy);
        }
        // impose water depth at the right side
        else if (scvf.center()[0] > 400.0 - eps_)
        {
            boundaryStateVariables =  ShallowWater::fixedWaterDepthBoundary(hBoundary_,
                                                                            insideVolVars.waterDepth(),
                                                                            insideVolVars.velocity(0),
                                                                            insideVolVars.velocity(1),
                                                                            gravity,
                                                                            nxy);
        }
        // no flow boundary
        else
        {
            // For the rough side walls of type no-slip we prescribe a slip condition based on alphaWall
            // for smooth closed walls we prescribed a full-slip condition (zero-roughness)

            // Get inside velocity components in cell face coordinates (t,n) using normal vector nxy
            // Note that the first component is the normal component
            // since we are rotating to the normal vector coordinate system
            Scalar insideVelocityNWall =  insideVolVars.velocity(0)*nxy[0] + insideVolVars.velocity(1)*nxy[1];
            Scalar insideVelocityTWall = -insideVolVars.velocity(0)*nxy[1] + insideVolVars.velocity(1)*nxy[0];

            // Initialisation of outside velocities
            auto outsideVelocityNWall = insideVelocityNWall;
            auto outsideVelocityTWall = insideVelocityTWall;

            // Now set the outside (ghost cell) velocities based on the chosen slip condition
            if ((scvf.center()[1] < -50.0 + eps_ || scvf.center()[1] > 50.0 - eps_) && (wallFrictionLawType_ == "noslip" || wallFrictionLawType_ == "nikuradse"))
            {
                // Set the outside state using the no-slip wall roughness conditions based on alphaWall
                // alphaWall = 0.0: full slip (wall tangential velocity equals inside tangential velocity)
                // alphaWall = 1.0: no slip (wall tangential velocity=0.0: point mirroring of the velocity)
                // 0.0 < alphaWall < 1.0: 'partial' slip.
                // e.g. for alphaWall = 0.5 the wall tangential velocity is half the inside tangential velocity.
                outsideVelocityNWall = -insideVelocityNWall;
                outsideVelocityTWall =  (1.0 - 2.0*alphaWall_)*insideVelocityTWall;
            }
            else
            {
                // Set the outside state using the full-slip wall roughness conditions (line mirroring in the boundary face)
                // Only mirror the normal component
                outsideVelocityNWall = -insideVelocityNWall;
                outsideVelocityTWall =  insideVelocityTWall;
            }
            // Rotate back to cartesian coordinate system
            Scalar outsideVelocityXWall = outsideVelocityNWall*nxy[0] - outsideVelocityTWall*nxy[1];
            Scalar outsideVelocityYWall = outsideVelocityNWall*nxy[1] + outsideVelocityTWall*nxy[0];

            boundaryStateVariables[0] = insideVolVars.waterDepth();
            boundaryStateVariables[1] = outsideVelocityXWall;
            boundaryStateVariables[2] = outsideVelocityYWall;

        }

        auto riemannFlux = ShallowWater::riemannProblem(insideVolVars.waterDepth(),
                                                        boundaryStateVariables[0],
                                                        insideVolVars.velocity(0),
                                                        boundaryStateVariables[1],
                                                        insideVolVars.velocity(1),
                                                        boundaryStateVariables[2],
                                                        insideVolVars.bedSurface(),
                                                        insideVolVars.bedSurface(),
                                                        gravity,
                                                        nxy);

        values[Indices::massBalanceIdx] = riemannFlux[0];
        values[Indices::velocityXIdx]   = riemannFlux[1];
        values[Indices::velocityYIdx]   = riemannFlux[2];

        // Addition of viscosity/diffusive flux rough wall boundaries
        // No-slip wall (with coefficient alphaWall):
        // Compute wall shear stress using turbulent viscosity and local velocity gradient
        // Assume velocity gradient (in cell adjacent to wall) equal to alphaWall*(0 - u_c)
        std::array<Scalar, 3> roughWallFlux;
        roughWallFlux = {0.0};
        if (scvf.center()[1] < -50.0 + eps_ || scvf.center()[1] > 50.0 - eps_)
        {
            // Distancess from the cell center to the wall
            auto dx = scvf.center()[0]-insideScv.center()[0];
            auto dy = scvf.center()[1]-insideScv.center()[1];
            if (wallFrictionLawType_ == "noslip")
            {
                roughWallFlux = ShallowWater::noslipWallBoundary(alphaWall_,
                                                                 insideVolVars.waterDepth(),
                                                                 insideVolVars.velocity(0),
                                                                 insideVolVars.velocity(1),
                                                                 turbViscosity_,
                                                                 dx,
                                                                 dy,
                                                                 nxy,
                                                                 fvGeometry,
                                                                 scvf);
            }
            else if (wallFrictionLawType_ == "nikuradse")
            {
                roughWallFlux = ShallowWater::nikuradseWallBoundary(ksWall_,
                                                                    insideVolVars.waterDepth(),
                                                                    insideVolVars.velocity(0),
                                                                    insideVolVars.velocity(1),
                                                                    dx,
                                                                    dy,
                                                                    nxy);
            }
        }

        values[Indices::massBalanceIdx] += roughWallFlux[0];
        values[Indices::velocityXIdx]   += roughWallFlux[1];
        values[Indices::velocityYIdx]   += roughWallFlux[2];

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param globalPos The position for which the boundary type is set
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        // the water depth is 10.0 m everywhere
        // the bed level runs from -9.98 m at the inflow boundary
        // to -10.0 m at the outflow boundary

        values[0] = hBoundary_;
        values[1] = -(9.81*bedSlope_/(8.0*turbViscosity_)) * (4.0*std::pow(globalPos[1],2.0) - std::pow(100.0,2.0));
        values[2] = 0.0;

        return values;
    };

    // \}

private:

    std::vector<Scalar> exactWaterDepth_;
    std::vector<Scalar> exactVelocityX_;
    std::vector<Scalar> exactVelocityY_;
    Scalar hBoundary_;                        // water level at the outflow boundary
    Scalar bedSlope_;                         // bed slope (positive downwards)
    Scalar discharge_;                        // discharge at the inflow boundary
    Scalar alphaWall_;                        // wall roughness coefficient for no-slip type wall roughness
    Scalar ksWall_;                           // Nikuradse wall roughness height for Nikuradse type wall roughness
    Scalar turbViscosity_;                    // turbulent viscosity
    static constexpr Scalar eps_ = 1.0e-6;
    std::string wallFrictionLawType_;         // wall friction law type
    std::string name_;
};

} //end namespace Dumux

#endif
