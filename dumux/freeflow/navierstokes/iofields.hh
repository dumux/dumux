// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIOFields
 */
#ifndef DUMUX_NAVIER_STOKES_IO_FIELDS_HH
#define DUMUX_NAVIER_STOKES_IO_FIELDS_HH

#include <dumux/common/parameters.hh>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief helper function to determine the names of cell-centered primary variables of a model with staggered grid discretization
 * \note use this as input for the load solution function
 */
template<class IOFields, class PrimaryVariables, class ModelTraits, class FluidSystem>
std::function<std::string(int,int)> createCellCenterPVNameFunction(const std::string& paramGroup = "")
{
    static constexpr auto offset = ModelTraits::numEq() - PrimaryVariables::dimension;

    if (hasParamInGroup(paramGroup, "LoadSolution.CellCenterPriVarNames"))
    {
        const auto pvName = getParamFromGroup<std::vector<std::string>>(paramGroup, "LoadSolution.CellCenterPriVarNames");
        return [n = std::move(pvName)](int pvIdx, int state = 0){ return n[pvIdx]; };
    }
    else
        // add an offset to the pvIdx so that the velocities are skipped
        return [](int pvIdx, int state = 0){ return IOFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx + offset, state); };
}

/*!
 * \ingroup NavierStokesModel
 * \brief helper function to determine the names of primary variables on the cell faces of a model with staggered grid discretization
 * \note use this as input for the load solution function
 */
template<class IOFields, class PrimaryVariables, class ModelTraits, class FluidSystem>
std::function<std::string(int,int)> createFacePVNameFunction(const std::string& paramGroup = "")
{
    if (hasParamInGroup(paramGroup, "LoadSolution.FacePriVarNames"))
    {
        const auto pvName = getParamFromGroup<std::vector<std::string>>(paramGroup, "LoadSolution.FacePriVarNames");
        return [n = std::move(pvName)](int pvIdx, int state = 0){ return n[pvIdx]; };
    }
    else
        // there is only one scalar-type primary variable called "v" living on the faces
        return [](int pvIdx, int state = 0){ return IOFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx , state); };
}

// forward declare
template<class T, class U>
class StaggeredVtkOutputModule;

/*!
 * \ingroup NavierStokesModel
 * \brief Adds I/O fields for the Navier-Stokes model
 */
class NavierStokesIOFields
{
    //! Helper strcuts to determine whether a staggered grid discretization is used
    template<class T>
    struct isStaggered : public std::false_type {};

    template<class... Args>
    struct isStaggered<StaggeredVtkOutputModule<Args...>>
    : public std::true_type {};

public:
    //! Initialize the Navier-Stokes specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addVolumeVariable([](const auto& v){ return v.pressure(); }, IOName::pressure());
        out.addVolumeVariable([](const auto& v){ return v.density(); }, IOName::density());

        // add discretization-specific fields
        additionalOutput_(out, isStaggered<OutputModule>());
    }

    //! return the names of the primary variables
    template <class ModelTraits, class FluidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        if (pvIdx < ModelTraits::dim())
            return "v";
        else
            return IOName::pressure();
    }

private:

    //! Adds discretization-specific fields (nothing by default).
    template <class OutputModule>
    static void additionalOutput_(OutputModule& out, std::false_type)
    { }

    //! Adds discretization-specific fields (velocity vectors on the faces for the staggered discretization).
    template <class OutputModule>
    static void additionalOutput_(OutputModule& out, std::true_type)
    {
        const bool writeFaceVars = getParamFromGroup<bool>(out.paramGroup(), "Vtk.WriteFaceData", false);
        if(writeFaceVars)
        {
            auto faceVelocityVector = [](const auto& scvf, const auto& faceVars)
                                      {
                                          using VelocityVector = std::decay_t<decltype(scvf.unitOuterNormal())>;

                                          VelocityVector velocity(0.0);
                                          velocity[scvf.directionIndex()] = faceVars.velocitySelf();
                                          return velocity;
                                      };

            out.addFaceVariable(faceVelocityVector, "faceVelocity");

            auto faceNormalVelocity = [](const auto& faceVars)
                                      {
                                          return faceVars.velocitySelf();
                                      };

            out.addFaceVariable(faceNormalVelocity, "v");
        }
    }
};

} // end namespace Dumux

#endif
