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
 *        all fluxes of fluid phases over a face of a finite volume,
 *        according to the Forchheimer-relation between velocity and pressure.
 *
 */
#ifndef DUMUX_BOX_FORCHHEIMER_FLUX_VARIABLES_HH
#define DUMUX_BOX_FORCHHEIMER_FLUX_VARIABLES_HH

#include <dumux/implicit/common/implicitforchheimerfluxvariables.hh>

namespace Dumux
{
    
/*!
 * \ingroup BoxModel
 * \ingroup BoxFluxVariables
 * \brief Evaluates the normal component of the Forchheimer velocity
 *        on a (sub)control volume face.
 *
 *        The commonly used Darcy relation looses it's validity for \f$ Re < 1\f$.
 *        If one encounters flow velocities in porous media above this Reynolds number,
 *        the Forchheimer relation can be used. Like the Darcy relation, it relates
 *        the gradient in potential to velocity.
 *        However, this relation is not linear (as in the Darcy case) any more.
 *
 *        Therefore, a Newton scheme is implemented in this class in order to calculate a
 *        velocity from the current set of variables. This velocity can subsequently be used
 *        e.g. by the localresidual.
 *
 *        For Reynolds numbers above \f$ 500\f$ the (Standard) forchheimer relation also
 *        looses it's validity.
 *
 *        The Forchheimer equation looks as follows:
 *         \f[ \nabla \left({p_\alpha + \rho_\alpha g z} \right)=  - \frac{\mu_\alpha}{k_{r \alpha} K} \underline{v_\alpha} - \frac{c_F}{\eta_{\alpha r} \sqrt{K}} \rho_\alpha |\underline{v_\alpha}| \underline{v_\alpha} \f]
 *
 *        For the formulation that is actually used in this class, see the documentation of the function calculating the residual.
 *
 *        - This algorithm does not find a solution if the fluid is incompressible and the
 *          initial pressure distribution is uniform.
 *        - This algorithm assumes the relative passabilities to have the same functional
 *          form as the relative permeabilities.
 */
template <class TypeTag>
class BoxForchheimerFluxVariables
    : public ImplicitForchheimerFluxVariables<TypeTag>
{
    typedef ImplicitForchheimerFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    
public:
    /*!
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
    DUNE_DEPRECATED_MSG("Use ImplicitForchheimerFluxVariables from "
                        "dumux/implicit/common/implicitforchheimerfluxvariables.hh.")
    BoxForchheimerFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const unsigned int faceIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
    :   ParentType(problem, element, fvGeometry, faceIdx, elemVolVars, onBoundary)
    {}
};

} // end namespace

#endif
