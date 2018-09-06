// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIOFields
 */
#ifndef DUMUX_NAVIER_STOKES_IO_FIELDS_HH
#define DUMUX_NAVIER_STOKES_IO_FIELDS_HH

#include <dune/common/fvector.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/indices.hh>
#include <dune/istl/multitypeblockvector.hh> // TODO: needed? or is forward declare enough?

#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/io/name.hh>

namespace Dumux
{

/*!
 * \ingroup InputOutput
 * \ingroup NavierStokesModel
 * \brief helper function to determine the names of cell-centered primary variables of a model with staggered grid discretization
 * \note use this as input for the load solution function
 */
template<class ModelTraits, class IOFields, class FluidSystem, std::size_t cellCenterIdx = 0, class ...Args, std::size_t i>
std::function<std::string(int,int)> createStaggeredPVNameFunction(const Dune::MultiTypeBlockVector<Args...>&,
                                                                  Dune::index_constant<i> idx,
                                                                  const std::string& paramGroup = "")
{
    using CellCenterPriVarsType = typename std::tuple_element_t<cellCenterIdx, std::tuple<Args...>>::block_type;
    static constexpr auto offset = ModelTraits::numEq() - CellCenterPriVarsType::dimension;

    if (i == cellCenterIdx) // cell center
    {
        if (hasParamInGroup(paramGroup, "LoadSolution.CellCenterPriVarNames"))
        {
            const auto pvName = getParamFromGroup<std::vector<std::string>>(paramGroup, "LoadSolution.CellCenterPriVarNames");
            return [n = std::move(pvName)](int pvIdx, int state = 0){ return n[pvIdx]; };
        }
        else
        // add an offset to the pvIdx so that the velocities are skipped
        return [](int pvIdx, int state = 0){ return IOFields::template primaryVariableName<ModelTraits ,FluidSystem>(pvIdx + offset, state); };
    }
    else // face
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
        out.addVolumeVariable([](const auto& v){ return v.molarDensity(); }, IOName::molarDensity());
        out.addVolumeVariable([](const auto& v){ return v.density(); }, IOName::density());

        // add discretization-specific fields
        additionalOutput_(out, isStaggered<OutputModule>());
    }

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
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
    static void additionalOutput_(OutputModule& out)
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
