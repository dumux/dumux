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
 *        specializations are provided for combinations of physical processes
 */
template<class TypeTag, bool enableAdvection, bool enableMolecularDiffusion, bool enableEnergyBalance>
class FluxVariables {};


// specialization for pure advective flow (e.g. 1p/2p/3p immiscible darcy flow)
template<class TypeTag>
class FluxVariables<TypeTag, true, false, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Stencil = std::vector<IndexType>;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);

public:
    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvf)
    {
        if (scvf.boundary())
        {
            if(!boundaryVolVars_)
                boundaryVolVars_ = Dune::Std::make_unique<VolumeVariables>();
            advection_.update(problem, element, scvf, boundaryVolVars_.get());
        }
        else
        {
            advection_.update(problem, element, scvf);
        }
    }

    void beginFluxComputation()
    {
        advection_.beginFluxComputation();
    }

    const AdvectionType& advection() const
    {
        return advection_;
    }

    const Stencil& stencil() const
    {
        return advection_.stencil();
    }

private:
    AdvectionType advection_;
    // boundary volume variables in case of Dirichlet boundaries
    std::unique_ptr<VolumeVariables> boundaryVolVars_;
};


// specialization for isothermal advection molecularDiffusion equations
template<class TypeTag>
class FluxVariables<TypeTag, true, true, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Stencil = std::vector<IndexType>;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvf)
    {
        if (scvf.boundary())
        {
            if(!boundaryVolVars_)
                boundaryVolVars_ = Dune::Std::make_unique<VolumeVariables>();

            advection_.update(problem, element, scvf, boundaryVolVars_.get());
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if (phaseIdx != compIdx)
                        molecularDiffusion(phaseIdx, compIdx).update(problem, element, scvf, phaseIdx, compIdx, boundaryVolVars_.get());
                }
        }
        else
        {
            advection_.update(problem, element, scvf);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if (phaseIdx != compIdx)
                        molecularDiffusion(phaseIdx, compIdx).update(problem, element, scvf, phaseIdx, compIdx);
                }
        }
    }

    void beginFluxComputation()
    {
        advection_.beginFluxComputation();
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                if (phaseIdx != compIdx)
                    molecularDiffusion(phaseIdx, compIdx).beginFluxComputation(true);
            }
        }
    }

    // TODO update transmissibilities?

    const AdvectionType& advection() const
    {
        return advection_;
    }

    AdvectionType& advection()
    {
        return advection_;
    }

    const MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx) const
    {
        if (compIdx < phaseIdx)
            return molecularDiffusion_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            return molecularDiffusion_[phaseIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion flux called for phaseIdx = compIdx");
    }

    MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx)
    {
        if (compIdx < phaseIdx)
            return molecularDiffusion_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            return molecularDiffusion_[phaseIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion flux called for phaseIdx = compIdx");
    }

    Stencil stencil() const
    {
        // std::set<IndexType> stencil(advection().stencil().begin(), advection().stencil().end());
        // TODO: insert all molecularDiffusion stencils of all components in all phases
        // stencil.insert();
        return advection().stencil();
    }

private:
    AdvectionType advection_;
    std::array< std::array<MolecularDiffusionType, numComponents-1>, numPhases> molecularDiffusion_;
    // boundary volume variables in case of Dirichlet boundaries
    std::unique_ptr<VolumeVariables> boundaryVolVars_;
};


// specialization for pure molecularDiffusion_
template<class TypeTag>
class FluxVariables<TypeTag, false, true, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                molecularDiffusion_.update(problem, scvf, phaseIdx, compIdx);
    }

    const MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx) const
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx)
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

private:
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> molecularDiffusion_;
};


// specialization for non-isothermal advective flow (e.g. non-isothermal one-phase darcy equation)
template<class TypeTag>
class FluxVariables<TypeTag, true, false, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        advection_.update(problem, scvf);
        heatConduction_.update(problem, scvf);
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


// specialization for non-isothermal advection molecularDiffusion_ equations
template<class TypeTag>
class FluxVariables<TypeTag, true, true, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
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
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        advection_.update(problem, scvf);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                molecularDiffusion_.update(problem, scvf, phaseIdx, compIdx);
        heatConduction_.update(problem, scvf);
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

    const MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx) const
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx)
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

private:
    AdvectionType advection_;
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> molecularDiffusion_;
    HeatConductionType heatConduction_;
};


// specialization for non-isothermal molecularDiffusion_
template<class TypeTag>
class FluxVariables<TypeTag, false, true, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                molecularDiffusion_.update(problem, scvf, phaseIdx, compIdx);
        heatConduction_.update(problem, scvf);
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

    const MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx) const
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx)
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

private:
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> molecularDiffusion_;
    HeatConductionType heatConduction_;
};


// specialization for pure heat conduction (e.g. the heat equation)
template<class TypeTag>
class FluxVariables<TypeTag, false, false, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        heatConduction_.update(problem, scvf);
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
