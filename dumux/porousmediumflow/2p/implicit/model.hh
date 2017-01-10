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
 * \brief Adaption of the fully implicit scheme to the two-phase flow model.
 */

#ifndef DUMUX_TWOP_MODEL_HH
#define DUMUX_TWOP_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPModel
 * \brief Adaption of the fully implicit scheme to the two-phase flow model.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPFormulation::pwsn</tt> or <tt>TwoPFormulation::pnsw</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */
template<class TypeTag >
class TwoPModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseModel);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem& problem)
    {
        ParentType::init(problem);

        // register standardized vtk output fields
        auto& vtkOutputModule = problem.vtkOutputModule();
        vtkOutputModule.addSecondaryVariable("Sw", [](const VolumeVariables& v){ return v.saturation(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("Sn", [](const VolumeVariables& v){ return v.saturation(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pw", [](const VolumeVariables& v){ return v.pressure(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pn", [](const VolumeVariables& v){ return v.pressure(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pc", [](const VolumeVariables& v){ return v.capillaryPressure(); });
        vtkOutputModule.addSecondaryVariable("rhoW", [](const VolumeVariables& v){ return v.density(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("rhoN", [](const VolumeVariables& v){ return v.density(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("mobW", [](const VolumeVariables& v){ return v.mobility(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("mobN", [](const VolumeVariables& v){ return v.mobility(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("temperature", [](const VolumeVariables& v){ return v.temperature(); });
        vtkOutputModule.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });
    }
};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
