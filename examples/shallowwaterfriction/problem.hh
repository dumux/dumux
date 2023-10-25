// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLE_SHALLOWWATER_FRICTION_PROBLEM_HH
#define DUMUX_EXAMPLE_SHALLOWWATER_FRICTION_PROBLEM_HH

// ## The problem file (`problem.hh`)
//
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the shallow water flow simulation with bottom friction.
// In addition, the analytical solution is defined here.
//
// [[content]]
//
// ### Include files
// [[details]] includes
//
// The first include we need here is the `ShallowWaterProblem` class, the base
// class from which we will derive.
#include <dumux/freeflow/shallowwater/problem.hh>
// In addition, we need the boundaryflux header, which handles the flux over
// the model boundaries.
#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>
// Include the `BoundaryTypes` class which specifies the boundary types set in this problem.
#include <dumux/common/boundarytypes.hh>
// Include the `NumEqVector` class which specifies a field vector with size number of equations in this problem.
#include <dumux/common/numeqvector.hh>
// [[/details]]

//
// ### The problem class
//
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// In addition the analytical solution of the problem is calculated.
// As this is a shallow water problem, we inherit from the basic ShallowWaterProblem.
// [[codeblock]]
namespace Dumux {

template <class TypeTag>
class RoughChannelProblem : public ShallowWaterProblem<TypeTag>
{
    // [[/codeblock]]
    // [[details]] convenience aliases
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using NeumannFluxes = NumEqVector;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    // [[/details]]

    // In the constructor, we retrieve all required parameters from the input file
    // [[codeblock]]
public:
    RoughChannelProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // We read the parameters from the params.input file.
        constManningN_ = getParam<Scalar>("Problem.ManningN");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        discharge_ = getParam<Scalar>("Problem.Discharge");
        // We calculate the outflow boundary condition using the Gauckler-Manning-Strickler formula.
        hBoundary_ = this->gaucklerManningStrickler(discharge_,constManningN_,bedSlope_);
        // We initialize the analytic solution to a vector of the appropriate size filled with zeros.
        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(gridGeometry->numDofs(), 0.0);
    }
    // [[/codeblock]]

    // #### Analytical Solution
    //
    // The analytical solution is calculated using the equation of Gauckler, Manning and Strickler.
    // [[codeblock]]

    // Equation of Gauckler, Manning and Strickler
    Scalar gaucklerManningStrickler(Scalar discharge, Scalar manningN, Scalar bedSlope)
    {
        using std::pow;
        using std::abs;
        using std::sqrt;

        return pow(abs(discharge)*manningN/sqrt(bedSlope), 0.6);
    }

    // Calculate the analytical solution
    void analyticalSolution()
    {
        using std::abs;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const Scalar h = this->gaucklerManningStrickler(discharge_,constManningN_,bedSlope_);
            const Scalar u = abs(discharge_)/h;

            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            exactWaterDepth_[eIdx] = h;
            exactVelocityX_[eIdx] = u;
        }
    }

    // [[exclude]]
    // Getter function for the analytical solution of the water depth
    const std::vector<Scalar>& getExactWaterDepth()
    {
        return exactWaterDepth_;
    }

    // Getter function for the analytical solution of the velocity in x-direction
    const std::vector<Scalar>& getExactVelocityX()
    {
        return exactVelocityX_;
    }
    // [[/exclude]]
    // [[/codeblock]]

    // #### Bottom friction
    //
    // The bottom friction is a source term and therefore handled by the `source` function.
    // [[codeblock]]
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {

        NumEqVector source (0.0);

        // Since the bed slope source term is handles within the flux computation,
        // in this model the bottom friction is the only source term.
        source += bottomFrictionSource(element, fvGeometry, elemVolVars, scv);

        return source;
    }
    // [[/codeblock]]

    // The calculation of the source term due to bottom friction needs the bottom shear stess.
    // This is the force per area, which works between the flow and the channel bed
    // (1D vector with two entries) and is calculated within the `FrictionLaw` class.
    // The bottom friction causes a loss of momentum. Thus, the first entry of the `bottomFrictionSource`,
    // which is related to the mass balance equation is zero.
    // The second entry of the `bottomFricitonSource` corresponds to the momentum equation in x-direction
    // and is therefore equal to the first, the x-component, of the `bottomShearStress`.
    // Accordingly, the third entry of the `bottomFrictionSource` is equal to the second component of the `bottomShearStress`.
    // [[codeblock]]
    NumEqVector bottomFrictionSource(const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const SubControlVolume &scv) const
    {
        NumEqVector bottomFrictionSource(0.0);
        const auto& volVars = elemVolVars[scv];

        // bottom shear stress vector
        Dune::FieldVector<Scalar, 2> bottomShearStress = this->spatialParams().frictionLaw(element, scv).bottomShearStress(volVars);

        // source term due to bottom friction
        bottomFrictionSource[Indices::massBalanceIdx] = 0.0;
        bottomFrictionSource[Indices::momentumXBalanceIdx] = -bottomShearStress[0] / volVars.density();
        bottomFrictionSource[Indices::momentumYBalanceIdx] = -bottomShearStress[1] / volVars.density();

        return bottomFrictionSource;
    }
    // [[/codeblock]]

    // #### Boundary conditions
    //
    // We use Neumann boundary conditions on the entire boundary.
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }
    // [[/codeblock]]

    // In the following function we implement the __Neumann boundary conditions__.
    // Due to the weak imposition we calculate the flux at the boundary with a Riemann solver.
    // This needs the state of a virtual cell outside of the boundary (`boundaryStateVariables`),
    // which is calculated with the Riemann invariants
    // (see: Yoon and Kang, "Finite Volume Model for Two-Dimensional Shallow Water Flows on Unstructured Grids").
    // [[codeblock]]
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

        // Calculate the Riemann invariants for imposed discharge at the left side.
        if (scvf.center()[0] < this->gridGeometry().bBoxMin()[0] + eps_)
        {
            boundaryStateVariables = ShallowWater::fixedDischargeBoundary(discharge_,
                                                                          insideVolVars.waterDepth(),
                                                                          insideVolVars.velocity(0),
                                                                          insideVolVars.velocity(1),
                                                                          gravity,
                                                                          nxy);
        }
        // Calculate the Riemann invariants for imposed water depth at the right side.
        else if (scvf.center()[0] > this->gridGeometry().bBoxMax()[0] - eps_)
        {
            boundaryStateVariables =  ShallowWater::fixedWaterDepthBoundary(hBoundary_,
                                                                            insideVolVars.waterDepth(),
                                                                            insideVolVars.velocity(0),
                                                                            insideVolVars.velocity(1),
                                                                            gravity,
                                                                            nxy);
        }
        // Calculate the Riemann invariants for the no-flow boundary.
        else
        {
            boundaryStateVariables[0] = insideVolVars.waterDepth();
            boundaryStateVariables[1] = -insideVolVars.velocity(0);
            boundaryStateVariables[2] = -insideVolVars.velocity(1);
        }
        // Calculate the boundary fluxes based on a Riemann problem.
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
    // [[/codeblock]]

    // #### Initial conditions
    //
    // We specify the initial conditions for the primary variables (water depth, velocity in y-direction
    // and velocity in x-direction). In this example constant initial conditions are used. Therefore the
    //argument `globalPos` is not needed. If you want to impose spatial variable initial conditions,
    // you have to use the `globalPos` argument.
    // [[codeblock]]
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);
        initialValues[Indices::waterdepthIdx] = 1.0;
        initialValues[Indices::velocityXIdx] = 0.0;
        initialValues[Indices::velocityYIdx] = 0.0;
        return initialValues;
    }
    // [[/codeblock]]

// [[details]] private variables
// [[codeblock]]
private:
    // variables for the analytic solution.
    std::vector<Scalar> exactWaterDepth_;
    std::vector<Scalar> exactVelocityX_;
    // constant friction value (an analytic solution is only available for const friction).
    Scalar constManningN_;
    // The constant channel bed slope.
    Scalar bedSlope_;
    // The water depth at the outflow boundary.
    Scalar hBoundary_;
    // The discharge at the inflow boundary.
    Scalar discharge_;
    // We assign a private global variable for the epsilon:
    static constexpr Scalar eps_ = 1.0e-6;

}; // end class definition RoughChannelProblem
} // end namespace Dumux
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
