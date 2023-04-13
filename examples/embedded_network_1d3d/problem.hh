// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_TISSUE_NETWORK_TRANSPORT_PROBLEMS_HH
#define DUMUX_TISSUE_NETWORK_TRANSPORT_PROBLEMS_HH

// # Initial, boundary conditions and sources (`problem.hh`)
//
// This file contains the __problem classes__ which defines the initial and boundary
// conditions and implements the coupling source terms. The file contains two
// problem classes: `TissueTransportProblem` and `NetworkTransportProblem` for the
// respective subdomains. The subdomain problem classes specify boundary and initial
// conditions for the subdomains separately. For this setup, we specify boundary fluxes
// via the `neumann` function (note that despite its name, the function allows you to
// implement both Neumann or Robin-type boundary conditions, weakly imposed).
//
// The subdomain problems are coupled to each other. This is evident from the
// coupling manager pointer that is stored in each subdomain problem (and the
// interface function `couplingManager()` giving access to the coupling manager object).
// This access is used in the `addPointSources()` and `pointSource(...)` functions.
// One point source (in this context) represents a quadrature point for the coupling
// condition integral. This means, we implement
//
// ```math
// \vert P \vert C_M D \bar{\varrho} (x^\bigcirc_\mathrm{T} - x_\mathrm{B}) w_i \mathrm{d}s_i
// ```
//
// where $`w_i`$ is the weight and $`\mathrm{d}s_i`$ (units of m) the integration element at
// the quadrature point. The units of the resulting term is mol/s (tracer amount per meter and second),
// which are the units of the one-dimensional balance equation integrated over one dimensional control volumes.
// We are integrating conceptually over the surface of the vessel. The quantity
// $`x^\bigcirc_\mathrm{T}`$ is the perimeter average of the 3D mole fraction field. This is not
// directly visible in the code since the chosen coupling manager determined what to return here
// for the 1D and 3D mole fraction depending on the coupling scheme.
//
// [[content]]
//
// ### Include headers
// We use properties (compile-time parameterization) and parameters (run-time parameterization).
// Moreover, we need array types and the problem base class for porous medium problems.
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>
//
// # The `TissueTransportProblem` class
// We start with the tissue transport with boundary, initial condition and sources.
// To set the initial condition, the problem reads some parameters from the configuration file
// `params.input`. Concentration is defined by $`c = \varrho x`$. This is usually a more
// convenient quantity which is why we read this from the configuration file.
// (What we read from the configuration file is entirely up to us. You can simply swap
// the parameter name in the `getParam` command and adapt the configuration file to provide the parameter.)
//
// [[codeblock]]
namespace Dumux {
template <class TypeTag>
class TissueTransportProblem : public PorousMediumFlowProblem<TypeTag>
{
// [[/codeblock]]
// [[details]] alias definitions and local variables
// [[codeblock]]
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ResidualVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimworld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr int phaseIdx = 0;
    static constexpr int compIdx = 0;

    static constexpr auto lowDimIdx = typename CouplingManager::MultiDomainTraits::template SubDomain<1>::Index();
// [[/codeblock]]
// [[/details]]
// [[codeblock]]
public:
    TissueTransportProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                           std::shared_ptr<CouplingManager> couplingManager,
                           const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_3d_transport";

        freeD_ = getParam<Scalar>("Tracer.DiffusionCoefficient");

        initialPeakConcentration_ = getParam<Scalar>("Tissue.Problem.InitialPeakConcentration", 1.0); // mmol/l
        initialCenter_ = getParam<GlobalPosition>("Tissue.Problem.InitialCenter");
        initialStddev_ = getParam<GlobalPosition>("Tissue.Problem.InitialStddev");
    }

    // name (used for output)
    const std::string& name() const
    { return name_; }
    // [[/codeblock]]
    //
    // We set the initial condition as a Gaussian
    //
    // ```math
    // x(t=0) = \frac{c}{\varrho} \operatorname{exp}\left( -\frac{x-x_0}{2 \sigma_x^2}
    // -\frac{y-y_0}{2 \sigma_y^2}-\frac{z-z_0}{2 \sigma_z^2} \right)
    //```
    //
    // if the Gaussian is very wide, this essentially corresponds to a constant value.
    //
    // [[codeblock]]
    // initial conditions (called when setting initial solution)
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables({
            initialPeakConcentration_/1000*0.018*std::exp(
                -(globalPos[0] - initialCenter_[0])*(globalPos[0] - initialCenter_[0])/(2*initialStddev_[0]*initialStddev_[0])
                -(globalPos[1] - initialCenter_[1])*(globalPos[1] - initialCenter_[1])/(2*initialStddev_[1]*initialStddev_[1])
                -(globalPos[2] - initialCenter_[2])*(globalPos[2] - initialCenter_[2])/(2*initialStddev_[2]*initialStddev_[2])
            )
        });
    }
    // [[/codeblock]]

    // We set the boundary condition to boundary flux (or Neumann conditions)
    // and specify zero flux everywhere for the tissue domain
    //
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    // flux boundaries (called for Neumann boundaries)
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    ResidualVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& fluxCache,
                           const SubControlVolumeFace& scvf) const
    {
        // no-flow conditions / symmetry conditions
        ResidualVector values(0.0);
        return values;
    }
    // [[/codeblock]]

    // Now follows the implementation of the point sources. The point sources
    // are obtained from (and computed by) the coupling manager. The function
    // `addPointSources` is called when we call `problem->computePointSourceMap()`
    // in the `main.cc`. The mechanism is provided by the base problem we inherit
    // from. Essentially all point sources computed by the coupling manager are
    // located and assigned to the respective control volume. Finally, the
    // function `pointSource` implements the actual value of the point source.
    // (While this is linear in the primary variable here, note that you can
    // also implement non-linear relations here, provided that a Newton solver
    // is used in the main program.) The `source(...)` function allows to
    // specify additional volumetric source terms (in units of mol/s/m^3),
    // for example, some metabolic demand.
    //
    // [[codeblock]]
    // this is needed for the 1d-3d coupling (realized via point source interface)
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().bulkPointSources(); }

    // called for every point source (coupling term integration point)
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // compute source at every integration point
        const Scalar x1D = this->couplingManager().lowDimPriVars(source.id())[compIdx];
        const Scalar x3D = this->couplingManager().bulkPriVars(source.id())[compIdx];

        // get the segment radius (outer radius)
        const Scalar radius = this->couplingManager().radius(source.id());

        // molar density (mol/m^3) (we assemble a mole balance equation instead of mass balance)
        const Scalar meanRhoMolar = 1000/0.018;

        // trans-membrane diffusive permeability
        const auto lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar lc = freeD_
            * this->couplingManager().problem(lowDimIdx).spatialParams().membraneFactor(lowDimElementIdx); // m/s

        // diffusive flux across membrane
        auto transMembraneFlux = -2*M_PI*radius*lc*meanRhoMolar*(x3D - x1D);
        transMembraneFlux *= source.quadratureWeight()*source.integrationElement();

        source = transMembraneFlux;
    }

    // additional volumetric source term in mol/(s*m^3), i.e. metabolic consumption
    // positive value means production into tissue, negative value means extraction from tissue
    template<class ElementVolumeVariables>
    ResidualVector source(const Element &element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolume &scv) const
    {
        return { 0.0 };
    }
    // [[/codeblock]]

    // Compute the current contrast agent amount in mole in the whole tissue
    // This function is used from the main file to write output.
    //
    // [[codeblock]]
    template<class SolutionVector, class GridVariables>
    Scalar computeTracerAmount(const SolutionVector& sol, const GridVariables& gridVars)
    {
        Scalar totalAmount(0.0);
        auto fvGeometry = localView(this->gridGeometry());
        auto elemVolVars = localView(gridVars.curGridVolVars());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                const auto localAmount = volVars.porosity()
                    *volVars.moleFraction(phaseIdx, compIdx)*volVars.molarDensity(phaseIdx)
                    *scv.volume()*volVars.extrusionFactor();
                totalAmount += localAmount;
            }
        }

        return totalAmount;
    }
    // [[/codeblock]]

// [[details]] coupling manager interface function and private variables
// [[codeblock]]
    // Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    Scalar freeD_;
    Scalar initialPeakConcentration_;
    GlobalPosition initialCenter_, initialStddev_;

    std::shared_ptr<CouplingManager> couplingManager_;
};
// [[/codeblock]]
// [[/details]]
//
// # The `NetworkTransportProblem` class
// We continue with the network transport with boundary, initial condition and sources.
// The implementation is quite similar to the tissue domain.
// This time the boundary function is more complex as we also take into account
// advective transport over the network boundaries to blood flow.
//
// [[codeblock]]
template <class TypeTag>
class NetworkTransportProblem : public PorousMediumFlowProblem<TypeTag>
{
// [[/codeblock]]
// [[details]] alias definitions and local variables
// [[codeblock]]
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ResidualVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimworld = GridView::dimensionworld;
    static constexpr int compIdx = 0;
    static constexpr int phaseIdx = 0;
// [[/codeblock]]
// [[/details]]
// [[codeblock]]
public:
    // initialize problem instance
    template<class GridData>
    NetworkTransportProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                          std::shared_ptr<CouplingManager> couplingManager,
                          std::shared_ptr<GridData> gridData,
                          const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_1d_transport";

        // initial conditions
        initialNetworkConcentration_ = getParam<Scalar>("Network.Problem.InitialConcentration", 0.0); // mmol/l

        // cross end-feet transport
        freeD_ = getParam<Scalar>("Tracer.DiffusionCoefficient");

        // clearance at sides by diffusion
        // a boundaryMembraneCoefficient_ of zero means zero gradient / zero diffusive flux
        // we can still have advective transport
        farFieldConcentration_ = getParam<Scalar>("Network.Problem.FarFieldConcentration", 0.0); // mmol/l
        boundaryMembraneCoefficient_ = getParam<Scalar>("Network.Problem.BoundaryMembraneCoefficient", 0.0); // dimensionless

        this->spatialParams().readGridParams(*gridData);
    }

    // name (used for output)
    const std::string& name() const
    { return name_; }
    // [[/codeblock]]
    // initial conditions (called when setting initial solution)
    // convert concentration to mole fraction
    // [[codeblock]]
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return { initialNetworkConcentration_/1000*0.018 }; }
    // [[/codeblock]]
    //
    // Set the boundary conditions. Note that we generalize the boundary conditions
    // for diffusion to allow for a Robin-type boundary condition with a "far-field"
    // concentration. This can be used to relax the constraint that there is absolutely
    // no diffusive flux over the boundary (`coeff` $`= 0`$).
    //
    // [[codeblock]]
    // type of boundary condition
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann(); // Robin / flux boundary conditions
        return bcTypes;
    }

    // flux boundaries (called for Neumann boundaries)
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    ResidualVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& fluxCache,
                           const SubControlVolumeFace& scvf) const
    {
        ResidualVector values(0.0);

        // Robin type boundary condition with far field concentration c
        // and conductance coefficient (dimensionless)
        // -> 1.0 means the concentration c reached in 1 µm distance from the boundary
        // -> if this coefficient is low this essentially means no-flow, if high this means direct outflow/inflow
        const auto robinFlux = [&](const auto& vv, const Scalar c, const Scalar coeff)
        {
            return - (c/vv.molarDensity() - vv.moleFraction(0, compIdx))
                     *vv.molarDensity()
                     *coeff * vv.effectiveDiffusionCoefficient(phaseIdx, phaseIdx, compIdx);
        };

        const auto& volVars = elemVolVars[scvf.insideScvIdx()];

        // diffusive part
        values[0] = robinFlux(volVars, farFieldConcentration_, boundaryMembraneCoefficient_);

        // advection
        const auto volumeFlux = this->spatialParams().volumeFlux(scvf.index());
        if (volumeFlux > 0) // outflow
            values[0] += volumeFlux*volVars.molarDensity()*volVars.moleFraction(0, compIdx)/(volVars.extrusionFactor()*scvf.area());

        return values;
    }
    // [[/codeblock]]
    //
    // The point source implementation is symmetric to the implementation of the
    // tissue domain. This is important such that the coupling condition is implemented
    // in a mass-conservative way.
    //
    // [[codeblock]]
    // this is needed for the 1d-3d coupling (realized via point source interface)
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().lowDimPointSources(); }

    // called for every point source (coupling term integration point)
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume& scv) const
    {
        // compute source at every integration point
        const Scalar x1D = this->couplingManager().lowDimPriVars(source.id())[compIdx];
        const Scalar x3D = this->couplingManager().bulkPriVars(source.id())[compIdx];

        // get the segment radius (outer PVS radius)
        const Scalar radius = this->couplingManager().radius(source.id());

        // molar density (mol/m^3) (we assemble a mole balance equation instead of mass balance)
        const Scalar meanRhoMolar = 1000/0.018;

        // trans-membrane diffusive permeability
        const Scalar lc = freeD_
            * this->spatialParams().membraneFactor(scv.elementIndex()); // m/s

        // diffusive flux across membrane
        auto transMembraneFlux = 2*M_PI*radius*lc*meanRhoMolar*(x3D - x1D);
        transMembraneFlux *= source.quadratureWeight()*source.integrationElement();

        source = transMembraneFlux;
    }

    // additional volumetric source term in mol/(s*m^3), i.e. metabolic consumption
    // positive value means production into blood, negative value means extraction from blood
    template<class ElementVolumeVariables>
    ResidualVector source(const Element &element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolume &scv) const
    {
        return { 0.0 };
    }
    // [[/codeblock]]
// [[details]] coupling manager function and internal variables (visibility: private)
// [[codeblock]]
    // The coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::string name_;

    // parameters
    Scalar freeD_;

    Scalar initialNetworkConcentration_; // mmol/l
    Scalar farFieldConcentration_; // mmol/l
    Scalar boundaryMembraneCoefficient_; // dimensionless

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
