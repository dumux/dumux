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
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the two-phase flow
 *        sub-problem in the coupled poro-mechanical elp problem.
 */

#ifndef DUMUX_2P_SUB_PROBLEM_HH
#define DUMUX_2P_SUB_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The two-phase sub problem in the el2p coupled problem.
 */
template <class TypeTag>
class TwoPSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // copy pressure index for convenience
    enum {
          pressureIdx = GetPropType<TypeTag, Properties::ModelTraits>::Indices::pressureIdx,
          saturationNIdx = GetPropType<TypeTag, Properties::ModelTraits>::Indices::saturationIdx,
          waterPhaseIdx = FluidSystem::phase0Idx,
          gasPhaseIdx = FluidSystem::phase1Idx,
          dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
public:
    TwoPSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams,
                   const std::string& paramGroup = "TwoP")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    {
        FluidSystem::init();
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Evaluates the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }

    //! Evaluates the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
      PrimaryVariables values;

      values[pressureIdx] = 1e6 * ((- 500 - globalPos[1])/2000*(24.6-5) + 5);
      values[saturationNIdx] = 0.0;
      return values;
    }

    //! Evaluates source terms.
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector values(0.0);
        if (!initialized_)
            return values;
        static const Scalar sourceG = getParam<Scalar>("Problem.InjectionRateGas");

        const GlobalPosition globalPos = scv.center();
        //std::cout << "Global Position: " << globalPos << std::endl;
        if(std::abs(globalPos[1] + 1500) < 3 && globalPos[0] < 80  )
        {
            //std::cout << "injection takes place" << std::endl;
            values[gasPhaseIdx] = sourceG; // kg/s
            values[gasPhaseIdx] /= scv.volume(); // kg/m3s
            values[gasPhaseIdx] /= FluidSystem::molarMass(gasPhaseIdx); // mol/m3s
        }
        return values;
    }

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

        if (globalPos[dimWorld-1] < eps_)
            values.setAllNeumann();

        return values;
    }

    void isInitialized(const bool& isInitialized)
    {
        initialized_ = isInitialized;
    }

private:
    static constexpr Scalar eps_ = 1.0e-6;
    std::string problemName_;
    bool initialized_ = false;
};

} // end namespace Dumux

#endif
