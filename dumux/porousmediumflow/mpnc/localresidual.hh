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
 * \ingroup MPNCModel
 * \brief MpNc specific details needed to approximately calculate the local
 *        defect in the fully implicit scheme.
 *
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_HH

#include <dune/istl/bvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup MPNCModel
 * \brief MpNc specific details needed to approximately calculate the local
 *        defect in the fully implicit scheme.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the
 * MpNc flow.
 */
template<class TypeTag>
class MPNCLocalResidual : public CompositionalLocalResidual<TypeTag>
{
    using ParentType = CompositionalLocalResidual<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementResidualVector = Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, NumEqVector)>;
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {phase0NcpIdx = Indices::phase0NcpIdx};

public:
    using ParentType::ParentType;


    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero for instationary problems.
     * \param problem The problem to solve

     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevVolVars The volume averaged variables for all
     *                    sub-control volumes of the element at the previous
     *                    time level
     * \param curVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all
     *                vertices of the element
     */
    ElementResidualVector eval(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& prevElemVolVars,
                               const ElementVolumeVariables& curElemVolVars,
                               const ElementBoundaryTypes &bcTypes,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());
        residual = 0.0;

        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            this->asImp().evalStorage(residual, problem, element, fvGeometry, prevElemVolVars, curElemVolVars, scv);
            this->asImp().evalSource(residual, problem, element, fvGeometry, curElemVolVars, scv);
             //here we need to set the constraints of the mpnc model into the residual
                const auto localScvIdx = scv.indexInElement();
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    residual[localScvIdx][phase0NcpIdx + phaseIdx] = curElemVolVars[scv].phaseNcp(phaseIdx);
                }
        }

        for (auto&& scvf : scvfs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            this->asImp().evalFlux(residual, problem, element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache, scvf);
        }

        return residual;
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero for stationary problem.
     * \param problem The problem to solve

     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param curVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all
     *                vertices of the element
     */
    ElementResidualVector eval(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& curElemVolVars,
                               const ElementBoundaryTypes &bcTypes,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());
        residual = 0.0;

        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            this->asImp().evalSource(residual, problem, element, fvGeometry, curElemVolVars, scv);

               //here we need to set the constraints of the mpnc model into the residual
                const auto localScvIdx = scv.indexInElement();
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    residual[localScvIdx][phase0NcpIdx + phaseIdx] = curElemVolVars[scv].phaseNcp(phaseIdx);
                }
        }

        for (auto&& scvf : scvfs(fvGeometry))
        {
            //! foward to the local residual specialized for the discretization methods
            this->asImp().evalFlux(residual, problem, element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache, scvf);
        }

        return residual;
    }
};

} // end namespace

#endif
