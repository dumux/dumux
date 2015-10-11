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
 * \brief The implicit non-isothermal model.
 */
#ifndef DUMUX_NI_MODEL_HH
#define DUMUX_NI_MODEL_HH

#include <dune/common/version.hh>

#include "niproperties.hh"

namespace Dumux {
/*!
 * \ingroup NIModel
 * \brief The implicit non-isothermal model.
 *
 * This model implements a generic energy balance for single and multi-phase
 * transport problems. Currently the non-isothermal model can be used on top of
 * the 1p2c, 2p, 2p2c and 3p3c models. Comparison to simple analytical solutions
 * for pure convective and conductive problems are found in the 1p2c test. Also refer
 * to this test for details on how to activate the non-isothermal model.
 *
 * For the energy balance, local thermal equilibrium is assumed. This
 * results in one energy conservation equation for the porous solid
 * matrix and the fluids:
 \f{align*}{
 \phi \frac{\partial \sum_\alpha \varrho_\alpha u_\alpha S_\alpha}{\partial t}
 & +
 \left( 1 - \phi \right) \frac{\partial (\varrho_s c_s T)}{\partial t}
 -
 \sum_\alpha \text{div}
 \left\{
 \varrho_\alpha h_\alpha
 \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left( \textbf{grad}\,p_\alpha - \varrho_\alpha \mbox{\bf g} \right)
 \right\} \\
    & - \text{div} \left(\lambda_{pm} \textbf{grad} \, T \right)
    - q^h = 0.
 \f}
 * where \f$h_\alpha\f$ is the specific enthalpy of a fluid phase
 * \f$\alpha\f$ and \f$u_\alpha = h_\alpha -
 * p_\alpha/\varrho_\alpha\f$ is the specific internal energy of the
 * phase.
 *
 */

template<class TypeTag>
class NIModel : public GET_PROP_TYPE(TypeTag, IsothermalModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, IsothermalModel) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        dim = GridView::dimension
    };
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        // call the ParentType output function
        ParentType::addOutputVtkFields(sol, writer);

        if (GET_PROP_VALUE(TypeTag, NiOutputLevel) == 0)
            return;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();

        // create required scalar fields
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        ScalarField &temperature = *writer.allocateManagedBuffer(numDofs);

        ElementIterator eIt = this->gridView().template begin<0>();
        ElementIterator eEndIt = this->gridView().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            if(eIt->partitionType() == Dune::InteriorEntity)
            {
                FVElementGeometry fvGeometry;
                fvGeometry.update(this->gridView_(), *eIt);

                ElementVolumeVariables elemVolVars;
                elemVolVars.update(this->problem_(),
                                   *eIt,
                                   fvGeometry,
                                   false /* oldSol? */);

                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int dofIdxGlobal = this->dofMapper().subIndex(*eIt, scvIdx, dofCodim);
                    temperature[dofIdxGlobal] = elemVolVars[scvIdx].temperature();
                }
            }
        }

        writer.attachDofData(temperature, "temperature", isBox);
    }
};

}
#endif
