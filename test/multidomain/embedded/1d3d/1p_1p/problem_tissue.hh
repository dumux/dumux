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
/**
 * \file
 * \ingroup EmbeddedTests
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of oxygen in interstitial fluid.
 */

#ifndef DUMUX_TISSUE_PROBLEM_HH
#define DUMUX_TISSUE_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh> // for coupling mode

#include "spatialparams_tissue.hh"

namespace Dumux {

template <class TypeTag>
class TissueProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct Tissue { using InheritsFrom = std::tuple<OneP>; };
struct TissueCC { using InheritsFrom = std::tuple<Tissue, CCTpfaModel>; };
struct TissueBox { using InheritsFrom = std::tuple<Tissue, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Tissue> { using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 3> >; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Tissue> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Tissue> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Tissue> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Tissue> { using type = TissueProblem<TypeTag>; };

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Tissue> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Tissue>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Tissue>
{
    using type = TissueSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                     GetPropType<TypeTag, Properties::Scalar>>;
};
} // end namespace Properties


/*!
 * \ingroup EmbeddedTests
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of oxygen in interstitial fluid.
 */
template <class TypeTag>
class TissueProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    TissueProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                  std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Tissue")
    , couplingManager_(couplingManager)
    {
        // read parameters from input file
        name_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        exactPressure_.resize(this->gridGeometry().numDofs());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
                exactPressure_[scv.dofIndex()] = exactSolution(scv.dofPosition());
        }
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 37 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 37; } // in [K]

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[2] > this->gridGeometry().bBoxMax()[2] - eps_  || globalPos[2] < this->gridGeometry().bBoxMin()[2] + eps_)
            values.setAllNeumann();
        else
            values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values = exactSolution(globalPos);
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables, class ElemFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElemFluxVarsCache&,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        // integrate over the scvf to compute the flux
        const auto geometry = scvf.geometry();
        Scalar derivative = 0.0;
        const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension-1>::rule(geometry.type(), 4);
        for(auto&& qp : quad)
        {
            auto globalPos = geometry.global(qp.position());
            globalPos[2] = 0; // the derivative in z-direction is the exact solution evaluated with z=0
            derivative += exactSolution(globalPos)*qp.weight()*geometry.integrationElement(qp.position());
        }

        const auto globalPos = scvf.ipGlobal();
        if (globalPos[2] > this->gridGeometry().bBoxMax()[2] - eps_ )
            flux[0] = -derivative/scvf.area();
        else
            flux[0] = derivative/scvf.area();
        return flux;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Applies a vector of point sources which are possibly solution dependent.
     *
     * \param pointSources A vector of Dumux::PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().bulkPointSources(); }

    /*!
     * \brief Evaluates the point sources (added by addPointSources)
     *        for all phases within a given sub control volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param source A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub-control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilated in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        if constexpr (CouplingManager::couplingMode != Embedded1d3dCouplingMode::kernel)
        {
            // compute source at every integration point
            const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];
            const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];

            // calculate the source
            const Scalar radius = this->couplingManager().radius(source.id());
            const Scalar beta = 2*M_PI/(2*M_PI + std::log(radius));
            const Scalar sourceValue = beta*(pressure1D - pressure3D);
            source = sourceValue*source.quadratureWeight()*source.integrationElement();
        }
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilated per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        if constexpr(CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const auto& sourceIds = this->couplingManager().bulkSourceIds(eIdx, scv.indexInElement());
            const auto& sourceWeights = this->couplingManager().bulkSourceWeights(eIdx, scv.indexInElement());

            for (int i = 0; i < sourceIds.size(); ++i)
            {
                const auto id = sourceIds[i];
                const auto weight = sourceWeights[i];
                const auto xi = this->couplingManager().fluxScalingFactor(id);

                const Scalar radius = this->couplingManager().radius(id);
                const Scalar pressure3D = this->couplingManager().bulkPriVars(id)[Indices::pressureIdx];
                const Scalar pressure1D = this->couplingManager().lowDimPriVars(id)[Indices::pressureIdx];

                // calculate the source
                static const Scalar beta = 2*M_PI/(2*M_PI + std::log(radius));
                const Scalar sourceValue = beta*(pressure1D - pressure3D)*xi;
                source[Indices::conti0EqIdx] += sourceValue*weight;
            }

            const auto volume = scv.volume()*elemVolVars[scv].extrusionFactor();
            source[Indices::conti0EqIdx] /= volume;
        }

        return source;
    }

    //! compute the flux scaling factor (xi) for a distance r for a vessel with radius R and kernel width rho
    Scalar fluxScalingFactor(const Scalar r, const Scalar R, const Scalar rho) const
    {
        using std::log;
        static const Scalar beta = 2*M_PI/(2*M_PI + log(R));
        static const Scalar kernelWidthFactorLog = log(rho/R);
        return r < rho ? 1.0/(1.0 + R*beta/(2*M_PI*R)*(r*r/(2*rho*rho) + kernelWidthFactorLog - 0.5))
                       : 1.0/(1.0 + R*beta/(2*M_PI*R)*log(r/R));
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    Scalar p1DExact(const GlobalPosition &globalPos) const
    { return 1.0 + globalPos[2]; }

    //! The exact solution
    Scalar exactSolution(const GlobalPosition &globalPos) const
    {
        const auto r = std::hypot(globalPos[0], globalPos[1]);
        static const auto R = getParam<Scalar>("SpatialParams.Radius");

        if (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
        {
            static const auto rho = getParam<Scalar>("MixedDimension.KernelWidthFactor")*R;
            if (r > rho)
                return -1.0*p1DExact(globalPos)/(2*M_PI)*std::log(r);
            else
                return -1.0*p1DExact(globalPos)/(2*M_PI)*(r*r/(2.0*rho*rho) + std::log(rho) - 0.5);
        }
        else if (CouplingManager::couplingMode == Embedded1d3dCouplingMode::surface)
        {
            if (r > R)
                return -1.0*p1DExact(globalPos)/(2*M_PI)*std::log(r);
            else
               return -1.0*p1DExact(globalPos)/(2*M_PI)*std::log(R);
        }

        return -1.0*p1DExact(globalPos)/(2*M_PI)*std::log(r);
    }

    //! Called after every time step
    //! Output the total global exchange term
    void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars)
    {
        PrimaryVariables source(0.0);
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                pointSources += this->source(element, fvGeometry, elemVolVars, scv);
                pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                source += pointSources;
            }
        }

        std::cout << "Global integrated source (3D): " << source << '\n';
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter.
     *
     * Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& vtk) const
    {
        vtk.addField(exactPressure_, "exact pressure");
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::vector<Scalar> exactPressure_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
