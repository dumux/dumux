// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later vesion.                                      *
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
 * \brief Non-isothermal injection problem where water is injected into a
 *        sand column with a NAPL contamination.
 */

#ifndef DUMUX_COLUMNXYLOLPROBLEM_HH
#define DUMUX_COLUMNXYLOLPROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCTests
 * \brief Non-isothermal injection problem where a water is injected into a
 *        sand column with a NAPL contamination.
 *
 * The 2D domain of this test problem is 0.1m long and 1.2m high.
 * Initially the column is filled with NAPL, gas and water, the NAPL saturation
 * increases to the bottom of the columns, the water saturation is constant.
 * Then water is injected from the top with a rate of 0.395710 mol/(s m) and
 * an enthalpy of 17452.97 [J/(s m)]. *
 *
 * Left, right and top boundaries are Neumann boundaries. Left and right are
 * no-flow boundaries, on top the injection takes places.
 * The bottom is a Dirichlet boundary.
 *
 * This problem uses the \ref ThreePThreeCModel and \ref NIModel model.
 *
 * This problem should typically be simulated for 200 days.
 * A good choice for the initial time step size is 1 s.
 * To adjust the simulation time it is necessary to edit the input file.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box3p3cnicolumnxylol test_columnxylol_fv.input</tt> or
 * <tt>./test_cc3p3cnicolumnxylol test_columnxylol_fv.input</tt>
 */
template <class TypeTag >
class ColumnProblem : public PorousMediumFlowProblem<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum
    {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        temperatureIdx = Indices::temperatureIdx,

        // equation indices
        energyEqIdx = Indices::energyEqIdx,
        contiWEqIdx = Indices::conti0EqIdx + FluidSystem::wCompIdx, //!< Index of the mass conservation equation for the water component,
        contiGEqIdx = Indices::conti0EqIdx + FluidSystem::gCompIdx, //!< Index of the mass conservation equation for the gas component,
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::nCompIdx, //!< Index of the mass conservation equation for the contaminant component

        // Phase State
        threePhases = Indices::threePhases,
        wgPhaseOnly = Indices::wgPhaseOnly
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:
    ColumnProblem(std::shared_ptr<const GridGeometry> gridGeometry)
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
        const auto yMax = this->gridGeometry().bBoxMax()[1];
        if (globalPos[1] < eps_ || globalPos[1] > yMax - eps_)
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
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wgPhaseOnly);

        priVars = initial_(globalPos);
        priVars[switch2Idx] = 0.0;
        return priVars;
//        return initial_(globalPos);
    }

    /*!
    * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

//         // negative values for injection
//         if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
//         {
//             values[contiWEqIdx] = -0.395710;
//             values[contiGEqIdx] = -0.000001;
//             values[contiNEqIdx] = -0.00;
//             values[energyEqIdx] = -17452.97;
//         }
        return values;
    }

    // \}


    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
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
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
//         const auto y = globalPos[1];
//         const auto yMax = this->gridGeometry().bBoxMax()[1];
        bool hasInitialXylenePhase = getParam<bool>("Problem.HasInitialXylenePhase", true);
        if (hasInitialXylenePhase)
            values.setState(threePhases);
        else
            values.setState(wgPhaseOnly);

        values[temperatureIdx] = getParam<Scalar>("Problem.Temperature", 296.15);
        values[pressureIdx] = 1.e5;
        values[switch1Idx] = getParam<Scalar>("Problem.InitialAirSaturation", 0.5);
        values[switch2Idx] = getParam<Scalar>("Problem.InitialXylene", 0.001);

        return values;
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    std::vector<Scalar> Kxx_;
};

} // end namespace Dumux

#endif
