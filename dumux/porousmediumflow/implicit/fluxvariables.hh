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
#ifndef DUMUX_POROUSMEDIUMFLOW_IMPLICIT_FLUXVARIABLES_HH
#define DUMUX_POROUSMEDIUMFLOW_IMPLICIT_FLUXVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(EnableAdvection);
NEW_PROP_TAG(EnableMolecularDiffusion);
NEW_PROP_TAG(EnableEnergyBalance);
}

// forward declaration
template<class TypeTag, bool enableAdvection, bool enableMolecularDiffusion, bool enableEnergyBalance>
class PorousMediumFluxVariablesImpl;

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables class
 *        specializations are provided for combinations of physical processes
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag>
using PorousMediumFluxVariables = PorousMediumFluxVariablesImpl<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                                         GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                                         GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;


// specialization for pure advective flow (e.g. 1p/2p/3p immiscible darcy flow)
template<class TypeTag>
class PorousMediumFluxVariablesImpl<TypeTag, true, false, false> : public FluxVariablesBase<TypeTag>
{
    using ParentType = FluxVariablesBase<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

public:

    template<typename FunctionType>
    Scalar advectiveFlux(const int phaseIdx, const FunctionType& upwindTerm)
    {
        Scalar flux = AdvectionType::flux(this->problem(),
                                          this->element(),
                                          this->fvGeometry(),
                                          this->elemVolVars(),
                                          this->scvFace(),
                                          phaseIdx,
                                          this->elemFluxVarsCache());

        return this->applyUpwindScheme(upwindTerm, flux, phaseIdx);
    }
};


// specialization for isothermal advection molecularDiffusion equations
template<class TypeTag>
class PorousMediumFluxVariablesImpl<TypeTag, true, true, false> : public FluxVariablesBase<TypeTag>
{
    using ParentType = FluxVariablesBase<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
           numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    //! The constructor
    PorousMediumFluxVariablesImpl()
    {
        advFluxCached_.reset();
        advPreFlux_.fill(0.0);
    }

    template<typename FunctionType>
    Scalar advectiveFlux(const int phaseIdx, const FunctionType& upwindTerm)
    {
        if (!advFluxCached_[phaseIdx])
        {

            advPreFlux_[phaseIdx] = AdvectionType::flux(this->problem(),
                                                        this->element(),
                                                        this->fvGeometry(),
                                                        this->elemVolVars(),
                                                        this->scvFace(),
                                                        phaseIdx,
                                                        this->elemFluxVarsCache());
            advFluxCached_.set(phaseIdx, true);
        }

        return this->applyUpwindScheme(upwindTerm, advPreFlux_[phaseIdx], phaseIdx);
    }

    Dune::FieldVector<Scalar, numComponents> molecularDiffusionFlux(const int phaseIdx)
    {
        return MolecularDiffusionType::flux(this->problem(),
                                            this->element(),
                                            this->fvGeometry(),
                                            this->elemVolVars(),
                                            this->scvFace(),
                                            phaseIdx,
                                            this->elemFluxVarsCache());
    }

private:
    //! simple caching if advection flux is used twice with different upwind function
    std::bitset<numPhases> advFluxCached_;
    std::array<Scalar, numPhases> advPreFlux_;
};

// specialization for non-isothermal advective flow (e.g. non-isothermal one-phase darcy equation)
template<class TypeTag>
class PorousMediumFluxVariablesImpl<TypeTag, true, false, true> : public FluxVariablesBase<TypeTag>
{
    using ParentType = FluxVariablesBase<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

public:
    //! The constructor
    PorousMediumFluxVariablesImpl()
    {
        advFluxCached_.reset();
        advPreFlux_.fill(0.0);
    }

    template<typename FunctionType>
    Scalar advectiveFlux(const int phaseIdx, const FunctionType& upwindTerm)
    {
        if (!advFluxCached_[phaseIdx])
        {

            advPreFlux_[phaseIdx] = AdvectionType::flux(this->problem(),
                                                        this->element(),
                                                        this->fvGeometry(),
                                                        this->elemVolVars(),
                                                        this->scvFace(),
                                                        phaseIdx,
                                                        this->elemFluxVarsCache());
            advFluxCached_.set(phaseIdx, true);
        }

        return this->applyUpwindScheme(upwindTerm, advPreFlux_[phaseIdx], phaseIdx);
    }

    Scalar heatConductionFlux()
    {
        return HeatConductionType::flux(this->problem(),
                                        this->element(),
                                        this->fvGeometry(),
                                        this->elemVolVars(),
                                        this->scvFace(),
                                        this->elemFluxVarsCache());
    }

private:
    //! simple caching if advection flux is used twice with different upwind function
    std::bitset<numPhases> advFluxCached_;
    std::array<Scalar, numPhases> advPreFlux_;
};

// specialization for non-isothermal advection and difussion equations (e.g. non-isothermal three-phase three-component flow)
template<class TypeTag>
class PorousMediumFluxVariablesImpl<TypeTag, true, true, true> : public FluxVariablesBase<TypeTag>
{
    using ParentType = FluxVariablesBase<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
           numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    //! The constructor
    PorousMediumFluxVariablesImpl()
    {
        advFluxCached_.reset();
        advPreFlux_.fill(0.0);
    }

    template<typename FunctionType>
    Scalar advectiveFlux(const int phaseIdx, const FunctionType& upwindTerm)
    {
        if (!advFluxCached_[phaseIdx])
        {

            advPreFlux_[phaseIdx] = AdvectionType::flux(this->problem(),
                                                        this->element(),
                                                        this->fvGeometry(),
                                                        this->elemVolVars(),
                                                        this->scvFace(),
                                                        phaseIdx,
                                                        this->elemFluxVarsCache());
            advFluxCached_.set(phaseIdx, true);
        }

        return this->applyUpwindScheme(upwindTerm, advPreFlux_[phaseIdx], phaseIdx);
    }

    Dune::FieldVector<Scalar, numComponents> molecularDiffusionFlux(const int phaseIdx)
    {
        return MolecularDiffusionType::flux(this->problem(),
                                            this->element(),
                                            this->fvGeometry(),
                                            this->elemVolVars(),
                                            this->scvFace(),
                                            phaseIdx,
                                            this->elemFluxVarsCache());
    }

    Scalar heatConductionFlux()
    {
        return HeatConductionType::flux(this->problem(),
                                        this->element(),
                                        this->fvGeometry(),
                                        this->elemVolVars(),
                                        this->scvFace(),
                                        this->elemFluxVarsCache());
    }

private:
    //! simple caching if advection flux is used twice with different upwind function
    std::bitset<numPhases> advFluxCached_;
    std::array<Scalar, numPhases> advPreFlux_;
};

} // end namespace

#endif
