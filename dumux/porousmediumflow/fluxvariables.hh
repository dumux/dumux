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
 * \brief Base class for the flux variables in porous medium models
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_FLUXVARIABLES_HH
#define DUMUX_POROUSMEDIUMFLOW_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief The porous medium flux variables class that computes advective / convective,
 *        molecular diffusive and heat conduction fluxes.
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag>
class PorousMediumFluxVariables : public FluxVariablesBase<TypeTag>
{
    using ParentType = FluxVariablesBase<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
           numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    static constexpr bool enableAdvection = GET_PROP_VALUE(TypeTag, EnableAdvection);
    static constexpr bool enableMolecularDiffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);
    static constexpr bool enableEnergyBalance = GET_PROP_VALUE(TypeTag, EnableEnergyBalance);
    static constexpr bool enableThermalNonEquilibrium = GET_PROP_VALUE(TypeTag, EnableThermalNonEquilibrium);

public:

    //! The constructor
    PorousMediumFluxVariables()
    {
        advFluxIsCached_.reset();
        advFluxBeforeUpwinding_.fill(0.0);
    }

    template<typename FunctionType>
    Scalar advectiveFlux(const int phaseIdx, const FunctionType& upwindTerm) const
    {
        if (enableAdvection)
        {
            if (!advFluxIsCached_[phaseIdx])
            {

                advFluxBeforeUpwinding_[phaseIdx] = AdvectionType::flux(this->problem(),
                                                                        this->element(),
                                                                        this->fvGeometry(),
                                                                        this->elemVolVars(),
                                                                        this->scvFace(),
                                                                        phaseIdx,
                                                                        this->elemFluxVarsCache());
                advFluxIsCached_.set(phaseIdx, true);
            }

            return this->applyUpwindScheme(upwindTerm, advFluxBeforeUpwinding_[phaseIdx], phaseIdx);
        }
        else
        {
            return 0.0;
        }
    }

    Dune::FieldVector<Scalar, numComponents> molecularDiffusionFlux(const int phaseIdx) const
    {
        if (enableMolecularDiffusion)
        {
            return MolecularDiffusionType::flux(this->problem(),
                                                this->element(),
                                                this->fvGeometry(),
                                                this->elemVolVars(),
                                                this->scvFace(),
                                                phaseIdx,
                                                this->elemFluxVarsCache());
        }
        else
        {
            return Dune::FieldVector<Scalar, numComponents>(0.0);
        }
    }

    template <bool enable = !enableThermalNonEquilibrium, typename std::enable_if_t<enable, int> = 0>
    Scalar heatConductionFlux() const
    {
        if (enableEnergyBalance)
        {
            return HeatConductionType::flux(this->problem(),
                                            this->element(),
                                            this->fvGeometry(),
                                            this->elemVolVars(),
                                            this->scvFace(),
                                            this->elemFluxVarsCache());
        }
        else
        {
            return 0.0;
        }
    }

    template <bool enable = enableThermalNonEquilibrium, typename std::enable_if_t<enable, int> = 0>
    Scalar heatConductionFlux(const int phaseIdx) const
    {
            return HeatConductionType::flux(this->problem(),
                                            this->element(),
                                            this->fvGeometry(),
                                            this->elemVolVars(),
                                            this->scvFace(),
                                            phaseIdx,
                                            this->elemFluxVarsCache());

    }

private:
    //! simple caching if advection flux is used twice with different upwind function
    mutable std::bitset<numPhases> advFluxIsCached_;
    mutable std::array<Scalar, numPhases> advFluxBeforeUpwinding_;
};

} // end namespace Dumux

#endif
