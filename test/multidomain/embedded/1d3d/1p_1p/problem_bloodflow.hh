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
 * \ingroup EmbeddedTests
 * \brief A test problem for the one-phase blood flow model:
 * Blood is flowing through a 1d network grid.
 */

#ifndef DUMUX_BLOOD_FLOW_PROBLEM_HH
#define DUMUX_BLOOD_FLOW_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh> // for coupling mode

#include "spatialparams_bloodflow.hh"

namespace Dumux {
// forward declaration
template <class TypeTag> class BloodFlowProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct BloodFlow { using InheritsFrom = std::tuple<OneP>; };
struct BloodFlowCC { using InheritsFrom = std::tuple<BloodFlow, CCTpfaModel>; };
struct BloodFlowBox { using InheritsFrom = std::tuple<BloodFlow, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::BloodFlow> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::BloodFlow> { using type = BloodFlowProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BloodFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BloodFlow> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BloodFlow>
{
    using type = BloodFlowSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
};
} // end namespace Properties

/*!
 * \ingroup EmbeddedTests
 * \brief Exact solution 1D-3D.
 */
template <class TypeTag>
class BloodFlowProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    BloodFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                     std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Vessel")
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        p_in_ = getParam<Scalar>("BoundaryConditions1D.PressureInput");
        delta_p_ = getParam<Scalar>("BoundaryConditions1D.DeltaPressure");
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
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 37.0; } // Body temperature

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The global position
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        if(globalPos[2] > 0.5)
            values[Indices::pressureIdx] = p_in_;
        else
            values[Indices::pressureIdx] = p_in_-delta_p_;
        return values;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

     /*!
     * \brief Applies a vector of point sources which are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().lowDimPointSources(); }

    /*!
     * \brief Evaluates the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
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
        // compute source at every integration point
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];

        // calculate the source
        const Scalar radius = this->couplingManager().radius(source.id());
        const Scalar beta = 2*M_PI/(2*M_PI + std::log(radius));
        Scalar sourceValue = beta*(pressure3D - pressure1D);
        if constexpr(CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
            sourceValue *= this->couplingManager().fluxScalingFactor(source.id());

        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    // \}

    //! The exact pressure solution
    Scalar exactSolution(const GlobalPosition &globalPos) const
    { return 1 + globalPos[2]; }

    //! Called after every time step
    //! Output the total global exchange term
    Scalar computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars)
    {
        SolutionVector source(sol.size());
        SolutionVector sourceExact(sol.size());
        Scalar volume = 0.0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                source[scv.dofIndex()] += pointSources;
                auto a = scv.corner(0)[2];
                auto b = scv.corner(1)[2];
                if (a > b) std::swap(a, b);
                sourceExact[scv.dofIndex()] += a*(1.0+0.5*a) - b*(1.0+0.5*b);
                volume += scv.volume();
            }
        }
        std::cout << "Global integrated source (1D): " << std::accumulate(source.begin(), source.end(), 0.0) << '\n';
        std::cout << "Global integrated source exact (1D): " << std::accumulate(sourceExact.begin(), sourceExact.end(), 0.0) << '\n';

        source -= sourceExact;
        const auto norm = source.two_norm()/sourceExact.two_norm();
        std::cout << "relative L2 Norm source (1D): " << norm << '\n';
        return source.two_norm()/volume;
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

    Scalar p_in_;
    Scalar delta_p_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
