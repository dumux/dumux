// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVEModel
 * \brief I/O container for storing fine-level quantities, used for vtk output
 */

#ifndef DUMUX_TWOPVE_FINE_LEVEL_FIELDSTORAGE_HH
#define DUMUX_TWOPVE_FINE_LEVEL_FIELDSTORAGE_HH

#include <cstddef>
#include <vector>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/porousmediumflow/2pve/elementstatefine.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVEModel
 * \brief Stores the reconstructed quantities of all fine-level elements for output
 *
 * \tparam GridGeometry the fine-level grid geometry
 * \tparam Scalar the scalar type
 * \tparam FluidSystem the immiscible two-phase fluid system
 */
template<class GridGeometry, class Scalar, class FluidSystem>
struct TwoPVEFineLevelFieldStorage
{
private:
    using FineLevelElementState = TwoPVEFineLevelElementState<GridGeometry, Scalar, FluidSystem>;
    static constexpr int wettingPhaseIdx = FluidSystem::phase0Idx;
    static constexpr int nonwettingPhaseIdx = FluidSystem::phase1Idx;
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
    void set(std::size_t fineIdx,
             const FineLevelElementState& fineElementState)
    {
        saturationW[fineIdx] = fineElementState.saturation(wettingPhaseIdx);
        saturationNw[fineIdx] = fineElementState.saturation(nonwettingPhaseIdx);
        mobilityW[fineIdx] = fineElementState.mobility(wettingPhaseIdx);
        mobilityNw[fineIdx] = fineElementState.mobility(nonwettingPhaseIdx);
        pressureW[fineIdx] = fineElementState.pressure(wettingPhaseIdx);
        pressureNw[fineIdx] = fineElementState.pressure(nonwettingPhaseIdx);
        capillaryPressure[fineIdx] = fineElementState.capillaryPressure();
        permeability[fineIdx] = fineElementState.permeability();
        porosity[fineIdx] = fineElementState.porosity();
        densityW[fineIdx] = fineElementState.density(wettingPhaseIdx);
        densityNw[fineIdx] = fineElementState.density(nonwettingPhaseIdx);
        gasPlumeDistance[fineIdx] = fineElementState.gasPlumeDist();
    }

    /*!
     * \brief Registers fine-level fields to vtk writer, should be called when initializing the vtk writer
     *
     * \param vtkSequenceWriter time-dependent vtkWriter
     */
    void registerFields(VTKSeqWriter& vtkSequenceWriter) const
    {
        vtkSequenceWriter.addCellData(saturationW, "S_liq");
        vtkSequenceWriter.addCellData(saturationNw, "S_gas");
        vtkSequenceWriter.addCellData(mobilityW, "mob_liq");
        vtkSequenceWriter.addCellData(mobilityNw, "mob_gas");
        vtkSequenceWriter.addCellData(pressureW, "p_liq");
        vtkSequenceWriter.addCellData(pressureNw, "p_gas");
        vtkSequenceWriter.addCellData(capillaryPressure, "pc");
        vtkSequenceWriter.addCellData(permeability, "permeability");
        vtkSequenceWriter.addCellData(porosity, "porosity");
        vtkSequenceWriter.addCellData(densityW, "rho_liq");
        vtkSequenceWriter.addCellData(densityNw, "rho_gas");
        vtkSequenceWriter.addCellData(gasPlumeDistance, "zp");
    }
};

} // end namespace Dumux

#endif
