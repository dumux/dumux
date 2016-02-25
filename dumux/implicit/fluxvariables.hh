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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_IMPLICIT_FLUXVARIABLES_HH
#define DUMUX_IMPLICIT_FLUXVARIABLES_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables
 *        specializations are provided for combinations of diffusion processes
 */
template<class TypeTag, bool darcy, bool diffusion, bool energy>
class FluxVariables {};


// specialization for pure darcy flow
template<class TypeTag>
class FluxVariables<TypeTag, true, false, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using DarcyFluxVariables = typename GET_PROP_TYPE(TypeTag, DarcyFluxVariables);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        darcyFluxVars_.update(problem, scv);
    }

    const DarcyFluxVariables& darcyFluxVars() const
    {
        return darcyFluxVars_;
    }

    DarcyFluxVariables& darcyFluxVars()
    {
        return darcyFluxVars_;
    }

private:
    DarcyFluxVariables darcyFluxVars_;
};


// specialization for darcy flow with diffusion
template<class TypeTag>
class FluxVariables<TypeTag, true, true, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using DarcyFluxVariables = typename GET_PROP_TYPE(TypeTag, DarcyFluxVariables);
    using DiffusionFluxVariables = typename GET_PROP_TYPE(TypeTag, DiffusionFluxVariables);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, numPhases);
        numComponents = GET_PROP_VALUE(TypeTag, numComponents);
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        darcyFluxVars_.update(problem, scv);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusionFluxVariables_.update(problem, scv, phaseIdx, compIdx);
    }

    const DarcyFluxVariables& darcyFluxVars() const
    {
        return darcyFluxVars_;
    }

    DarcyFluxVariables& darcyFluxVars()
    {
        return darcyFluxVars_;
    }

    const DiffusionFluxVariables& diffusionFluxVars(const int phaseIdx, const int compIdx) const
    {
        return diffusionFluxVariables_[phaseIdx][compIdx];
    }

    DiffusionFluxVariables& diffusionFluxVars(const int phaseIdx, const int compIdx)
    {
        return diffusionFluxVariables_[phaseIdx][compIdx];
    }

private:
    DarcyFluxVariables darcyFluxVars_;
    std::array< std::array<DiffusionFluxVariables, numComponents>, numPhases> diffusionFluxVariables_;
};


// specialization for pure diffusion
template<class TypeTag>
class FluxVariables<TypeTag, false, true, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using DiffusionFluxVariables = typename GET_PROP_TYPE(TypeTag, DiffusionFluxVariables);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, numPhases);
        numComponents = GET_PROP_VALUE(TypeTag, numComponents);
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusionFluxVariables_.update(problem, scv, phaseIdx, compIdx);
    }

    const DiffusionFluxVariables& diffusionFluxVars(const int phaseIdx, const int compIdx) const
    {
        return diffusionFluxVariables_[phaseIdx][compIdx];
    }

    DiffusionFluxVariables& diffusionFluxVars(const int phaseIdx, const int compIdx)
    {
        return diffusionFluxVariables_[phaseIdx][compIdx];
    }

private:
    std::array< std::array<DiffusionFluxVariables, numComponents>, numPhases> diffusionFluxVariables_;
};


// specialization for non-isothermal darcy flow
template<class TypeTag>
class FluxVariables<TypeTag, true, false, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using DarcyFluxVariables = typename GET_PROP_TYPE(TypeTag, DarcyFluxVariables);
    using EnergyFluxVariables = typename GET_PROP_TYPE(TypeTag, EnergyFluxVariables);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        darcyFluxVars_.update(problem, scv);
        energyFluxVars_.update(problem, scv);
    }

    const DarcyFluxVariables& darcyFluxVars() const
    {
        return darcyFluxVars_;
    }

    DarcyFluxVariables& darcyFluxVars()
    {
        return darcyFluxVars_;
    }

    const EnergyFluxVariables& energyFluxVars() const
    {
        return energyFluxVars_;
    }

    EnergyFluxVariables& energyFluxVars()
    {
        return energyFluxVars_;
    }

private:
    DarcyFluxVariables darcyFluxVars_;
    EnergyFluxVariables energyFluxVars_;
};


// specialization for non-isothermal darcy flow with diffusion
template<class TypeTag>
class FluxVariables<TypeTag, true, true, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using DarcyFluxVariables = typename GET_PROP_TYPE(TypeTag, DarcyFluxVariables);
    using DiffusionFluxVariables = typename GET_PROP_TYPE(TypeTag, DiffusionFluxVariables);
    using EnergyFluxVariables = typename GET_PROP_TYPE(TypeTag, EnergyFluxVariables);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, numPhases);
        numComponents = GET_PROP_VALUE(TypeTag, numComponents);
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        darcyFluxVars_.update(problem, scv);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusionFluxVariables_.update(problem, scv, phaseIdx, compIdx);
        energyFluxVars_.update(problem, scv);
    }

    const DarcyFluxVariables& darcyFluxVars() const
    {
        return darcyFluxVars_;
    }

    DarcyFluxVariables& darcyFluxVars()
    {
        return darcyFluxVars_;
    }

    const EnergyFluxVariables& energyFluxVars() const
    {
        return energyFluxVars_;
    }

    EnergyFluxVariables& energyFluxVars()
    {
        return energyFluxVars_;
    }

    const DiffusionFluxVariables& diffusionFluxVars(const int phaseIdx, const int compIdx) const
    {
        return diffusionFluxVariables_[phaseIdx][compIdx];
    }

    DiffusionFluxVariables& diffusionFluxVars(const int phaseIdx, const int compIdx)
    {
        return diffusionFluxVariables_[phaseIdx][compIdx];
    }

private:
    DarcyFluxVariables darcyFluxVars_;
    std::array< std::array<DiffusionFluxVariables, numComponents>, numPhases> diffusionFluxVariables_;
    EnergyFluxVariables energyFluxVars_;
};


// specialization for non-isothermal pure diffusion
template<class TypeTag>
class FluxVariables<TypeTag, false, true, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using DiffusionFluxVariables = typename GET_PROP_TYPE(TypeTag, DiffusionFluxVariables);
    using EnergyFluxVariables = typename GET_PROP_TYPE(TypeTag, EnergyFluxVariables);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, numPhases);
        numComponents = GET_PROP_VALUE(TypeTag, numComponents);
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusionFluxVariables_.update(problem, scv, phaseIdx, compIdx);
        energyFluxVars_.update(problem, scv);
    }

    const EnergyFluxVariables& energyFluxVars() const
    {
        return energyFluxVars_;
    }

    EnergyFluxVariables& energyFluxVars()
    {
        return energyFluxVars_;
    }

    const DiffusionFluxVariables& diffusionFluxVars(const int phaseIdx, const int compIdx) const
    {
        return diffusionFluxVariables_[phaseIdx][compIdx];
    }

    DiffusionFluxVariables& diffusionFluxVars(const int phaseIdx, const int compIdx)
    {
        return diffusionFluxVariables_[phaseIdx][compIdx];
    }

private:
    std::array< std::array<DiffusionFluxVariables, numComponents>, numPhases> diffusionFluxVariables_;
    EnergyFluxVariables energyFluxVars_;
};


// specialization for pure heat transport
template<class TypeTag>
class FluxVariables<TypeTag, false, false, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using EnergyFluxVariables = typename GET_PROP_TYPE(TypeTag, EnergyFluxVariables);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        energyFluxVars_.update(problem, scv);
    }

    const EnergyFluxVariables& energyFluxVars() const
    {
        return energyFluxVars_;
    }

    EnergyFluxVariables& energyFluxVars()
    {
        return energyFluxVars_;
    }

private:
    EnergyFluxVariables energyFluxVars_;
};
} // end namespace

#endif
