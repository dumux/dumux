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
 * \ingroup MultiDomainModel
 * \brief The Darcy local residual to compute the interface velocity
 *
 * For coupled Stokes/Darcy systems, the free flow velocity can be computed
 * directly out of the mass flux from the porous medium into the free flow domain.
 * If this option is chosen as a boundary condition in the Stokes sub problem,
 * the Darcy local residual needs the method 'computeIntegralFluxAcrossBoundary'.
 */
#ifndef DUMUX_DARCY_COUPLING_LOCAL_RESIDUAL_HH
#define DUMUX_DARCY_COUPLING_LOCAL_RESIDUAL_HH

#include <dumux/porousmediumflow/immiscible/localresidual.hh> // TODO 1p
//#include <dumux/porousmediumflow/compositional/localresidual.hh> // TODO 2p2c

namespace Dumux
{
namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(StokesProblemTypeTag);
NEW_PROP_TAG(DarcyProblemTypeTag);
}
/*!
 * \ingroup MultiDomainModel
 * \brief TODO
 *
 * This class is also used for the non-isothermal model, which means
 * that it uses static polymorphism.
 */

template<class TypeTag>
class DarcyCouplingLocalResidual : public ImmiscibleLocalResidual<TypeTag> // TODO change name
//class DarcyCouplingLocalResidual : public CompositionalLocalResidual<TypeTag> // TODO change name
{
protected:
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry) ;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

public:
     /*!
     * \brief Evaluates the integral flux leaving or entering the domain at a given element
     *        by partially evaluating the local residual
     *
     * \param element The DUNE Codim<0> entity for which the residual ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     *
     * \param prevVolVars The volume averaged variables for all sub-control volumes of the element at the previous
     *                   time level
     * \param curVolVars The volume averaged variables for all sub-control volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all vertices of the element
     *
     * \param returnValue Specifies whether the integral value shall be computed and returned.
     *                    If set false, the user has to do the actual evaluation of the actual flux.
     */
    PrimaryVariables computeIntegralFluxAcrossBoundary(const Element &element,
                                                       const FVElementGeometry& fvGeometry,
                                                       const ElementVolumeVariables& prevElemVolVars,
                                                       const ElementVolumeVariables& curElemVolVars,
                                                       const ElementBoundaryTypes &bcTypes,
                                                       const ElementFluxVariablesCache& elemFluxVarsCache,
                                                       const bool returnValue = true)
    {
        // resize the vectors for all terms
        const auto numScv = fvGeometry.numScv();
        this->residual_.resize(numScv);
        this->storageTerm_.resize(numScv);

        this->residual_ = 0.0;
        this->storageTerm_ = 0.0;

        // eval volume terms and fluxes, but NOT boundary conditions
        this->evalVolumeTerms_(element, fvGeometry, prevElemVolVars, curElemVolVars, bcTypes);
        this->evalFluxes_(element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache);

        PrimaryVariables flux(0.0);

        if(returnValue)
        {
            // only consider those scvs for the calculation of the flux
            for(auto&& scv : scvs(fvGeometry))
            {
                if(this->problem().model().onBoundary(scv))
                    flux += this->residual_[scv.indexInElement()];
            }
        }
        return flux;
    }
};

}




#endif // DUMUX_DARCY_COUPLING_LOCAL_RESIDUAL_HH
