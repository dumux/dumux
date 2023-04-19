// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Chemistry
 * \brief Electrochemical model for a fuel cell application.
 */
#ifndef DUMUX_ELECTROCHEMISTRY_NI_HH
#define DUMUX_ELECTROCHEMISTRY_NI_HH

#include <dumux/material/constants.hh>
#include <dumux/material/chemistry/electrochemistry/electrochemistry.hh>

namespace Dumux {

/*!
 * \ingroup Chemistry
 * \brief Class calculating source terms and current densities for fuel cells
 * with the electrochemical models suggested by Ochs (2008) \cite ochs2008 or Acosta (2006) \cite A3:acosta:2006
 * for the non-isothermal case.
 * \todo TODO: Scalar type should be extracted from VolumeVariables!
 * \todo TODO: This shouldn't depend on discretization and grid!!
 */
template <class Scalar, class Indices, class FluidSystem, class GridGeometry, ElectroChemistryModel electroChemistryModel>
class ElectroChemistryNI : public ElectroChemistry<Scalar, Indices, FluidSystem, GridGeometry, electroChemistryModel>
{
    using ParentType = ElectroChemistry<Scalar, Indices, FluidSystem, GridGeometry, electroChemistryModel>;
    using Constant = Constants<Scalar>;

    enum {
        //equation indices
        contiH2OEqIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        contiO2EqIdx = Indices::conti0EqIdx + FluidSystem::O2Idx,
        energyEqIdx = Indices::energyEqIdx, //energy equation
    };

    using GridView = typename GridGeometry::GridView;
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;
    using GlobalPosition = typename Dune::FieldVector<typename GridView::ctype, GridView::dimensionworld>;
    using CellVector = typename Dune::FieldVector<typename GridView::ctype, GridView::dimension>;

public:
    /*!
     * \brief Calculates reaction sources with an electrochemical model approach.
     *
     * \param values The primary variable vector
     * \param currentDensity The current density
     * \param paramGroup The group containing the needed parameters
     *
     * For this method, the \a values parameter stores source values
     */
    template<class SourceValues>
    static void reactionSource(SourceValues &values, Scalar currentDensity,
                               const std::string& paramGroup = "")
    {
        //correction to account for actually relevant reaction area
        //current density has to be divided by the half length of the box
        //\todo Do we have multiply with the electrochemically active surface area (ECSA) here instead?
        static Scalar gridYMax = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.UpperRight")[1];
        static Scalar nCellsY = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.Cells")[1];

        // Warning: This assumes the reaction layer is always just one cell (cell-centered) or half a box (box) thick
        const auto lengthBox = gridYMax/nCellsY;
        if (isBox)
            currentDensity *= 2.0/lengthBox;
        else
            currentDensity *= 1.0/lengthBox;

        static Scalar transportNumberH2O = getParam<Scalar>("ElectroChemistry.TransportNumberH20");
        static Scalar thermoneutralVoltage = getParam<Scalar>("ElectroChemistry.ThermoneutralVoltage");
        static Scalar cellVoltage = getParam<Scalar>("ElectroChemistry.CellVoltage");

        //calculation of flux terms with faraday equation
        values[contiH2OEqIdx] = currentDensity/(2*Constant::F);                  //reaction term in reaction layer
        values[contiH2OEqIdx] += currentDensity/Constant::F*transportNumberH2O;  //osmotic term in membrane
        values[contiO2EqIdx]  = -currentDensity/(4*Constant::F);                 //O2-equation
        values[energyEqIdx] = (thermoneutralVoltage - cellVoltage)*currentDensity; //energy equation
    }
};

} // end namespace Dumux

#endif
