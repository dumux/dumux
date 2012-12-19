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
 *
 * \brief This file contains the data which is required to calculate
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_BOX_DARCY_FLUX_VARIABLES_HH
#define DUMUX_BOX_DARCY_FLUX_VARIABLES_HH

#include <dumux/implicit/common/implicitdarcyfluxvariables.hh>

namespace Dumux
{
    
/*!
 * \ingroup BoxModel
 * \ingroup BoxFluxVariables
 * \brief Evaluates the normal component of the Darcy velocity 
 * on a (sub)control volume face.
 */
template <class TypeTag>
class BoxDarcyFluxVariables : public ImplicitDarcyFluxVariables<TypeTag>
{
    typedef ImplicitDarcyFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    /*
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param faceIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    DUNE_DEPRECATED_MSG("Use ImplicitDarcyFluxVariables from "
                         "dumux/implicit/common/implicitdarcyfluxvariables.hh.")
    BoxDarcyFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int faceIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
    : ParentType(problem, element, fvGeometry, faceIdx, elemVolVars, onBoundary)
    {}
};

} // end namespace

#endif
