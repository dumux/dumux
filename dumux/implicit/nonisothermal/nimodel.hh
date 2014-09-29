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
 * \brief the implicit non-isothermal model
 */
#ifndef DUMUX_NI_MODEL_HH
#define DUMUX_NI_MODEL_HH

#include "niproperties.hh"

namespace Dumux {
/*!
 * \ingroup NIModel
 *  *
 * This model implements a non-isothermal flow for single and multi-phase
 * transport problems
 *
 * \todo The equations have to be added here
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
                    int globalIdx = this->dofMapper().map(*eIt, scvIdx, dofCodim);
                    temperature[globalIdx] = elemVolVars[scvIdx].temperature();
                }
            }
        }

        writer.attachDofData(temperature, "temperature", isBox);
    }
};

}
#endif
