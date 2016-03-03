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

namespace Properties
{
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
}

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables
 *        specializations are provided for combinations of diffusion_ processes
 */
template<class TypeTag, bool enableAdvection, bool enableMolecularDiffusion, bool enableEnergyBalance>
class FluxVariables {};


// specialization for pure advective flow (e.g. one-phase darcy equation)
template<class TypeTag>
class FluxVariables<TypeTag, true, false, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Stencil = std::set<IndexType>;
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        advection_.update(problem, scv);
    }

    const AdvectionType& advection() const
    {
        return advection_;
    }

    Stencil stencil() const
    {
        return advection().stencil();
    }

private:
    AdvectionType advection_;
};


// specialization for isothermal advection diffusion_ equations
template<class TypeTag>
class FluxVariables<TypeTag, true, true, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        advection_.update(problem, scv);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusion_.update(problem, scv, phaseIdx, compIdx);
    }

    const AdvectionType& advection() const
    {
        return advection_;
    }

    AdvectionType& advection()
    {
        return advection_;
    }

    const MolecularDiffusionType& diffusion(const int phaseIdx, const int compIdx) const
    {
        return diffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& diffusion(const int phaseIdx, const int compIdx)
    {
        return diffusion_[phaseIdx][compIdx];
    }

private:
    AdvectionType advection_;
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> diffusion_;
};


// specialization for pure diffusion_
template<class TypeTag>
class FluxVariables<TypeTag, false, true, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusion_.update(problem, scv, phaseIdx, compIdx);
    }

    const MolecularDiffusionType& diffusion(const int phaseIdx, const int compIdx) const
    {
        return diffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& diffusion(const int phaseIdx, const int compIdx)
    {
        return diffusion_[phaseIdx][compIdx];
    }

private:
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> diffusion_;
};


// specialization for non-isothermal advective flow (e.g. non-isothermal one-phase darcy equation)
template<class TypeTag>
class FluxVariables<TypeTag, true, false, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        advection_.update(problem, scv);
        heatConduction_.update(problem, scv);
    }

    const AdvectionType& advection() const
    {
        return advection_;
    }

    AdvectionType& advection()
    {
        return advection_;
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

private:
    AdvectionType advection_;
    HeatConductionType heatConduction_;
};


// specialization for non-isothermal advection diffusion_ equations
template<class TypeTag>
class FluxVariables<TypeTag, true, true, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        advection_.update(problem, scv);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusion_.update(problem, scv, phaseIdx, compIdx);
        heatConduction_.update(problem, scv);
    }

    const AdvectionType& advection() const
    {
        return advection_;
    }

    AdvectionType& advection()
    {
        return advection_;
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

    const MolecularDiffusionType& diffusion(const int phaseIdx, const int compIdx) const
    {
        return diffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& diffusion(const int phaseIdx, const int compIdx)
    {
        return diffusion_[phaseIdx][compIdx];
    }

private:
    AdvectionType advection_;
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> diffusion_;
    HeatConductionType heatConduction_;
};


// specialization for non-isothermal diffusion_
template<class TypeTag>
class FluxVariables<TypeTag, false, true, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusion_.update(problem, scv, phaseIdx, compIdx);
        heatConduction_.update(problem, scv);
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

    const MolecularDiffusionType& diffusion(const int phaseIdx, const int compIdx) const
    {
        return diffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& diffusion(const int phaseIdx, const int compIdx)
    {
        return diffusion_[phaseIdx][compIdx];
    }

private:
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> diffusion_;
    HeatConductionType heatConduction_;
};


// specialization for pure heat conduction (e.g. the heat equation)
template<class TypeTag>
class FluxVariables<TypeTag, false, false, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scv)
    {
        heatConduction_.update(problem, scv);
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

private:
    HeatConductionType heatConduction_;
};

} // end namespace

#endif
