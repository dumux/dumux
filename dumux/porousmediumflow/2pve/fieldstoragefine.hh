// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVE
 * \brief I/O container for storing fine-level quantities, used for vtk output
 */

#ifndef DUMUX_TWOPVE_FINE_LEVEL_FIELDSTORAGE_HH
#define DUMUX_TWOPVE_FINE_LEVEL_FIELDSTORAGE_HH

#include <cstddef>
#include <vector>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/porousmediumflow/2pve/elementstatefine.hh>

namespace Dumux {

template<class TypeTag>
struct TwoPVEFineLevelFieldStorage
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FineLevelElementState = TwoPVEFineLevelElementState<TypeTag>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum {
        wettingPhaseIdx = FluidSystem::phase0Idx,
        nonwettingPhaseIdx = FluidSystem::phase1Idx,
    };
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using VTKSeqWriter = typename Dune::VTKSequenceWriter<GridView>;

public:
    std::vector<Scalar> saturationW;
    std::vector<Scalar> saturationNw;
    std::vector<Scalar> mobilityW;
    std::vector<Scalar> mobilityNw;
    std::vector<Scalar> pressureW;
    std::vector<Scalar> pressureNw;
    std::vector<Scalar> capillaryPressure;
    std::vector<Scalar> permeability;
    std::vector<Scalar> porosity;
    std::vector<Scalar> densityW;
    std::vector<Scalar> densityNw;
    std::vector<Scalar> gasPlumeDistance;

    explicit TwoPVEFineLevelFieldStorage(std::size_t size)
    : saturationW(size),
      saturationNw(size),
      mobilityW(size),
      mobilityNw(size),
      pressureW(size),
      pressureNw(size),
      capillaryPressure(size),
      permeability(size),
      porosity(size),
      densityW(size),
      densityNw(size),
      gasPlumeDistance(size)
    {}

    /*!
     * \brief Sets the fine-level fields to the provided fine-level element state
     *
     * \param fineIdx          index of fine-level element
     * \param fineElementState state of fine-level element
     */
    void set(const int fineIdx,
             const FineLevelElementState& fineElementState)
    {
        this->saturationW[fineIdx] = fineElementState.saturation(wettingPhaseIdx);
        this->saturationNw[fineIdx] = fineElementState.saturation(nonwettingPhaseIdx);
        this->mobilityW[fineIdx] = fineElementState.mobility(wettingPhaseIdx);
        this->mobilityNw[fineIdx] = fineElementState.mobility(nonwettingPhaseIdx);
        this->pressureW[fineIdx] = fineElementState.pressure(wettingPhaseIdx);
        this->pressureNw[fineIdx] = fineElementState.pressure(nonwettingPhaseIdx);
        this->capillaryPressure[fineIdx] = fineElementState.capillaryPressure();
        this->permeability[fineIdx] = fineElementState.permeability();
        this->porosity[fineIdx] = fineElementState.porosity();
        this->densityW[fineIdx] = fineElementState.density(wettingPhaseIdx);
        this->densityNw[fineIdx] = fineElementState.density(nonwettingPhaseIdx);
        this->gasPlumeDistance[fineIdx] = fineElementState.gasPlumeDist();
    }

    /*!
     * \brief Registers fine-level fields to vtk writer, should be called when initializing the vtk writer
     *
     * \param vtkSequenceWriter time-dependant vtkWriter
     */
    void registerFields(VTKSeqWriter& vtkSequenceWriter) const
    {
        vtkSequenceWriter.addCellData(this->saturationW, "S_liq");
        vtkSequenceWriter.addCellData(this->saturationNw, "S_gas");
        vtkSequenceWriter.addCellData(this->mobilityW, "mob_liq");
        vtkSequenceWriter.addCellData(this->mobilityNw, "mob_gas");
        vtkSequenceWriter.addCellData(this->pressureW, "p_liq");
        vtkSequenceWriter.addCellData(this->pressureNw, "p_gas");
        vtkSequenceWriter.addCellData(this->capillaryPressure, "pc");
        vtkSequenceWriter.addCellData(this->permeability, "permeability");
        vtkSequenceWriter.addCellData(this->porosity, "porosity");
        vtkSequenceWriter.addCellData(this->densityW, "rho_liq");
        vtkSequenceWriter.addCellData(this->densityNw, "rho_gas");
        vtkSequenceWriter.addCellData(this->gasPlumeDistance, "zp");
    }
};

} // end namespace Dumux

#endif
