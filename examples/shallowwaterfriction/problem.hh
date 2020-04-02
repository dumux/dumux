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

#ifndef DUMUX_EXAMPLE_SHALLOWWATER_FRICTION_PROBLEM_HH
#define DUMUX_EXAMPLE_SHALLOWWATER_FRICTION_PROBLEM_HH

// ## The file `problem.hh`
// We start with includes
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>

// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a shallow water problem, we inherit from the basic ShallowWaterProblem.
namespace Dumux {

template <class TypeTag>
class RoughChannelProblem : public ShallowWaterProblem<TypeTag>
{
    // We use convenient declarations that we derive from the property system.
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
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    // This is the constructor of our problem class.
    RoughChannelProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // We read the parameters from the params.input file.
        name_ = getParam<std::string>("Problem.Name");
        constManningN_ = getParam<Scalar>("Problem.ManningN");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        discharge_ = getParam<Scalar>("Problem.Discharge");
        // We calculate the outflow boundary condition using the Gauckler-Manning-Strickler formula.
        hBoundary_ = this->gauklerManningStrickler(discharge_,constManningN_,bedSlope_);
        // We initialize the analytic solution to a verctor of the appropriate size filled with zeros.
        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(gridGeometry->numDofs(), 0.0);
    }

    // Get the analytical water depth
    const std::vector<Scalar>& getExactWaterDepth()
    {
        return exactWaterDepth_;
    }

    // Get the analytical velocity
    const std::vector<Scalar>& getExactVelocityX()
    {
        return exactVelocityX_;
    }

    // Get the water depth with Gauckler-Manning-Strickler
    Scalar gauklerManningStrickler(Scalar discharge, Scalar manningN, Scalar bedSlope)
    {
        using std::pow;
        using std::abs;
        using std::sqrt;

        return pow(abs(discharge)*manningN/sqrt(bedSlope), 0.6);
    }

    // Get the analytical solution
    void analyticalSolution()
    {
        using std::abs;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const Scalar h = this->gauklerManningStrickler(discharge_,constManningN_,bedSlope_);
            const Scalar u = abs(discharge_)/h;

            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            exactWaterDepth_[eIdx] = h;
            exactVelocityX_[eIdx] = u;
        }
    }

    // Get the problem name. It is used as a prefix for files generated by the simulation.
    const std::string& name() const
    {
        return name_;
    }

    // Get the source term.
     NumEqVector source(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolume &scv) const
    {

        NumEqVector source (0.0);

        // In this model the bottom friction is the only source.
        source += bottomFrictionSource(element, fvGeometry, elemVolVars, scv);

        return source;
    }

     // Get the source term due to bottom friction.
     NumEqVector bottomFrictionSource(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const SubControlVolume &scv) const
     {
        NumEqVector bottomFrictionSource(0.0);
        const auto& volVars = elemVolVars[scv];

        // For the calculation of the source term due to bottom friction the two-dimensional bottom shear stess vector is needed. This is the force per area, which works between the flow and the bed. It is calculated within the `FrictionLaw`, which is a spatialParameter. In this model the `FrictionLawManning` is used (see `params.input`).
        Dune::FieldVector<Scalar, 2> bottomShearStress = this->spatialParams().frictionLaw(element, scv).shearStress(volVars);

        // The bottom shear stress causes a pure loss of momentum. Thus the first entry of the `bottomFrictionSource`, which is related to the mass balance equation is zero. The second entry of the `bottomFricitonSource` corresponds to the momentum equation in x-direction and is therefore equal to the first, the x-component, of the `bottomShearStress`. Accordingly the third entry of the `bottomFrictionSource` is equal to the second component of the `bottomShearStress`.
        bottomFrictionSource[0] = 0.0;
        bottomFrictionSource[1] = bottomShearStress[0];
        bottomFrictionSource[2] = bottomShearStress[1];

        return bottomFrictionSource;
     }

    // We specify the boundary condition type.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        // Since we use a weak imposition all boundary conditions are of Neumann type.
        bcTypes.setAllNeumann();
        return bcTypes;
    }

     // We specify the neumann boundary. Due to the weak imposition we calculate the flux at the
     // boundary, with a  Rieman solver. For this the state of a virtual cell outside of the boundary
     // is needed (`boundaryStateVariables`), wich is calculated with the Riemann invariants
     // (see Yoon and Kang, Finite Volume Model for Two-Dimensional Shallow Water Flows on Unstructured Grids).
     // The calculation of the Riemann invariants differ depending on the type of the boundary (h, q or no-flow boundary).
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

        // Calculate the rieman invariants for imposed discharge at the left side.
        if (scvf.center()[0] < 0.0 + eps_)
        {
            boundaryStateVariables = ShallowWater::fixedDischargeBoundary(discharge_,
                                                                          insideVolVars.waterDepth(),
                                                                          insideVolVars.velocity(0),
                                                                          insideVolVars.velocity(1),
                                                                          gravity,
                                                                          nxy);
        }
        // Calculate the rieman invariants for impose water depth at the right side.
        else if (scvf.center()[0] > 100.0 - eps_)
        {
            boundaryStateVariables =  ShallowWater::fixedWaterDepthBoundary(hBoundary_,
                                                                            insideVolVars.waterDepth(),
                                                                            insideVolVars.velocity(0),
                                                                            insideVolVars.velocity(1),
                                                                            gravity,
                                                                            nxy);
        }
        // Calculate the rieman invarianty for the no-flow boundary.
        else
        {
            boundaryStateVariables[0] = insideVolVars.waterDepth();
            boundaryStateVariables[1] = -insideVolVars.velocity(0);
            boundaryStateVariables[2] = -insideVolVars.velocity(1);
        }
        // We calculate the boundary fluxes based on a Riemann problem.
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

        return values;
    }

    // We set the initial conditions. In this example constant initial conditions are used. Therefore the argument `globalPos` is not needed. If you want to impose spatial variable initial conditions, you have to use the `globalPos`.
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        // We set the initial water depth to one meter.
        values[0] = 1.0;
        // We set the x-component of the initial velocity to zero.
        values[1] = 0.0;
        // We set the y-component of the initial velocity to zero.
        values[2] = 0.0;

        return values;
    };

    // \}

private:
    // We declare the private variables of the problem. They are initialized in the problems constructor.
    // We declare the variable for the analytic solution.
    std::vector<Scalar> exactWaterDepth_;
    std::vector<Scalar> exactVelocityX_;
    // constant friction value. An analytic solution is only available for const friction. If you want to run the simulation with a non constant friciton value (specified in the spatialParams) you have to remove the analytic solution.
    Scalar constManningN_;
    // The constant bed slope.
    Scalar bedSlope_;
    // The water depth at the outflow boundary.
    Scalar hBoundary_;
    // The discharge at the inflow boundary.
    Scalar discharge_;
    // eps is used as a small value for the definition of the boundry conditions
    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};

} // end namespace Dumux

#endif
