// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief Calculates the element-wise residual for the box scheme
 */
#ifndef DUMUX_FACECENTERED_LOCAL_RESIDUAL_HH
#define DUMUX_FACECENTERED_LOCAL_RESIDUAL_HH

#include <dune/geometry/type.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/assembly/fvlocalresidual.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief The element-wise residual for the box scheme
 * \tparam TypeTag the TypeTag
 */
template<class TypeTag>
class FaceCenteredLocalResidual : public FVLocalResidual<TypeTag>
{
    using ParentType = FVLocalResidual<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Extrusion = Extrusion_t<GridGeometry>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    //! evaluate flux residuals for one sub control volume face and add to residual
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const ElementFluxVariablesCache& elemFluxVarsCache,
                  const SubControlVolumeFace& scvf) const
    {
        const auto flux = evalFlux(problem, element, fvGeometry, elemVolVars, elemBcTypes, elemFluxVarsCache, scvf);
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        residual[insideScv.localDofIndex()] += flux;
    }

    //! evaluate flux residuals for one sub control volume face
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementBoundaryTypes& elemBcTypes,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        if (elemBcTypes.hasNeumann())
            return this->asImp().computeFluxWithNeumannBoundaries(problem, element, fvGeometry, elemVolVars, elemBcTypes, elemFluxVarsCache, scvf);


        return this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
    }

    using ParentType::evalStorage;

    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVolVars The volume averaged variables for all
     *                        sub-control volumes of the element at the previous time level
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     * \param scv The sub control volume the storage term is integrated over
     */
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const SubControlVolume& scv) const
    {
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];

        // mass balance within the element. this is the
        // \f$\frac{m}{\partial t}\f$ term if using implicit or explicit
        // euler as time discretization.
        //

        //! Compute storage with the model specific storage residual
        NumEqVector prevStorage = this->asImp().computeStorage(problem, scv, prevVolVars, true/*isPreviousStorage*/);
        NumEqVector storage = this->asImp().computeStorage(problem, scv, curVolVars, false/*isPreviousStorage*/);

        prevStorage *= prevVolVars.extrusionFactor();
        storage *= curVolVars.extrusionFactor();

        storage -= prevStorage;
        storage *= Extrusion::volume(scv);
        storage /= this->timeLoop().timeStepSize();

        residual[scv.localDofIndex()] += storage;
    }
};

} // end namespace Dumux

#endif
