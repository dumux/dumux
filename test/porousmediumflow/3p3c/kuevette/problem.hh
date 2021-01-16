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
 * \ingroup ThreePThreeCTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 */

#ifndef DUMUX_KUEVETTE3P3CNIPROBLEM_HH
#define DUMUX_KUEVETTE3P3CNIPROBLEM_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3p3c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

#define ISOTHERMAL 0

namespace Dumux {

/*!
 * \ingroup ThreePThreeCTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 */
template <class TypeTag>
class KuevetteProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Kuevette { using InheritsFrom = std::tuple<ThreePThreeCNI>; };
struct KuevetteBox { using InheritsFrom = std::tuple<Kuevette, BoxModel>; };
struct KuevetteCCTpfa { using InheritsFrom = std::tuple<Kuevette, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Kuevette> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Kuevette> { using type = KuevetteProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Kuevette>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = KuevetteSpatialParams<GridGeometry, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Kuevette>
{ using type = FluidSystems::H2OAirMesitylene<GetPropType<TypeTag, Properties::Scalar>>; };
} // end namespace Properties

/*!
 * \ingroup ThreePThreeCTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 *
 * The domain is a quasi-two-dimensional container (kuevette). Its dimensions
 * are 1.5m x 0.74m. The top and bottom boundaries are closed, the right
 * boundary is a Dirichlet boundary allowing fluids to escape. From the left,
 * an injection of a hot water-air mixture is applied (Neumann boundary condition
 * for the mass components and the enthalpy), aimed at remediating an initial
 * NAPL (Non-Aquoeus Phase Liquid) contamination in the heterogeneous domain.
 * The contamination is initially placed partly into the coarse sand
 * and partly into a fine sand lens.
 *
 * This simulation can be varied through assigning different boundary conditions
 * at the left boundary as described in \cite Class2001.
 *
 * This problem uses the \ref ThreePThreeCModel and \ref NIModel model.
 *
 * To see the basic effect and the differences to scenarios with pure steam or
 * pure air injection, it is sufficient to simulate for about 2-3 hours (10000 s).
 * Complete remediation of the domain requires much longer (about 10 days simulated time).
 * To adjust the simulation time it is necessary to edit the input file.
 *
 * To run the simulation execute:
 * <tt>./test_box3p3cnikuevette test_box3p3cnikuevette.input</tt> or
 * <tt>./test_cc3p3cnikuevette test_cc3p3cnikuevette.input</tt>
 */
template <class TypeTag >
class KuevetteProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum {

        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        temperatureIdx = Indices::temperatureIdx,
        contiWEqIdx = Indices::conti0EqIdx + FluidSystem::wCompIdx,
        contiGEqIdx = Indices::conti0EqIdx + FluidSystem::gCompIdx,
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::nCompIdx,
        energyEqIdx = Indices::energyEqIdx,

        // phase states
        threePhases = Indices::threePhases,
        wgPhaseOnly = Indices::wgPhaseOnly,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    KuevetteProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        name_ = getParam<std::string>("Problem.Name");
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
    const std::string name() const
    { return name_; }

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
        BoundaryTypes bcTypes;
        if(globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
         return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos( const GlobalPosition &globalPos) const
    {
       return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a N eumann boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = scvf.ipGlobal();

        // negative values for injection
        if (globalPos[0] < eps_)
        {
            values[contiWEqIdx] = -0.1435; // 0.3435 [mol/(s m)] in total
            values[contiGEqIdx] = -0.2;
            values[contiNEqIdx] =  0.0;
            values[energyEqIdx] = -6929.;
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos( const GlobalPosition &globalPos) const
    {
       return initial_(globalPos);
    }

    /*!
     * \brief Appends all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * Adjust this in case of anisotropic permeabilities.
     */
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
        const auto& gg = this->gridGeometry();
        Kxx_.resize(gg.numDofs());
        vtk.addField(Kxx_, "permeability");

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(gg);
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
                Kxx_[scv.dofIndex()] = this->spatialParams().intrinsicPermeabilityAtPos(scv.dofPosition());
        }
    }

private:
    // checks, whether a point is located inside the contamination zone
    bool isInContaminationZone(const GlobalPosition &globalPos) const
    {
        return (Dune::FloatCmp::ge<Scalar>(globalPos[0], 0.2)
               && Dune::FloatCmp::le<Scalar>(globalPos[0], 0.8)
               && Dune::FloatCmp::ge<Scalar>(globalPos[1], 0.4)
               && Dune::FloatCmp::le<Scalar>(globalPos[1], 0.65));
    }

    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        if (isInContaminationZone(globalPos))
            values.setState(threePhases);
        else
            values.setState(wgPhaseOnly);

        values[pressureIdx] = 1e5 ;
        values[switch1Idx] = 0.12;
        values[switch2Idx] = 1.e-6;
        values[temperatureIdx] = 293.0;

        if (isInContaminationZone(globalPos))
        {
            values[switch2Idx] = 0.07;
        }
        return values;
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    std::vector<Scalar> Kxx_;
};

} // end namespace Dumux

#endif
