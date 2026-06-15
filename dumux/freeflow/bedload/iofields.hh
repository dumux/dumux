// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransportModel
 * \copydoc Dumux::BedloadIOFields
 */
#ifndef DUMUX_BEDLOAD_IO_FIELDS_HH
#define DUMUX_BEDLOAD_IO_FIELDS_HH

#include <dumux/io/name.hh>
#include <string>

namespace Dumux {

/*!
 * \ingroup BedloadTransportModel
 * \brief Adds vtk output fields for the bedload transport model
 *
 * The bedload model contains the shallow water variables as well as the bedload
 * variables. Therefore, it can be used to write both.
 *
 * Only the variables of the active layer are added, as the layers beneath are
 * handled by the layer model and not by the volume variables. If needed, the
 * variables of the layers beneath can be added to the vtk output using the addField
 * function of the vtkwriter in the main file.
 */
class BedloadIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;

        out.addVolumeVariable([](const VolumeVariables& v){ return v.bedSurface(); }, "bedSurface");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.bedSurface() + v.waterDepth(); }, "freeSurface");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.sedimentMass(); }, "massActiveLayer");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.layerThickness(); }, "ticknessActiveLayer");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.velocity(0); }, "velocityX");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.velocity(1); }, "velocityY");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.waterDepth(); }, "waterDepth");

        static const int nGrainClasses = Dumux::getParam<int>("Sediment.NumberGrainClasses");
        for (int i = 0; i < nGrainClasses; i++) {
            out.addVolumeVariable([i](const VolumeVariables& v){ return v.bedloadDischarge(i)[0]; }, "bedloadDischargeXGrainClass"+std::to_string(i+1));
            out.addVolumeVariable([i](const VolumeVariables& v){ return v.bedloadDischarge(i)[1]; }, "bedloadDischargeYGrainClass"+std::to_string(i+1));
            out.addVolumeVariable([i](const VolumeVariables& v){ return v.massFraction(i); }, "massFractionActiveLayerGrainClass"+std::to_string(i+1));
        }
    }
};

} // end namespace Dumux

#endif
