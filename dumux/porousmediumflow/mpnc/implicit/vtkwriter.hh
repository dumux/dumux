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
 * \brief Writes the VTK output files for a solution of the MpNc model.
 */

#ifndef DUMUX_MPNC_VTK_WRITER_HH
#define DUMUX_MPNC_VTK_WRITER_HH

#include "properties.hh"

#include <dumux/io/vtkmultiwriter.hh>

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \brief Writes the VTK output files for a
 * solution of the MpNc model.
 */
template<class TypeTag>
class MPNCVtkWriter
{
    using MPNCVtkCommonModule = typename GET_PROP_TYPE(TypeTag, MPNCVtkCommonModule);
    using MPNCVtkMassModule = typename GET_PROP_TYPE(TypeTag, MPNCVtkMassModule);
    using MPNCVtkEnergyModule = typename GET_PROP_TYPE(TypeTag, MPNCVtkEnergyModule);
    using MPNCVtkCustomModule = typename GET_PROP_TYPE(TypeTag, MPNCVtkCustomModule);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

public:
    MPNCVtkWriter(const Problem &problem)
        : problem_(problem)
        , commonWriter_(problem)
        , massWriter_(problem)
        , energyWriter_(problem)
        , customWriter_(problem)
    {
    }

    /*!
     * \brief Add  the current solution to the VtkMultiWriter.
     */
    template <class MultiWriter>
    void addCurrentSolution(MultiWriter &writer)
    {
        // tell sub-writers to allocate their buffers
        commonWriter_.allocBuffers(writer);
        massWriter_.allocBuffers(writer);
        energyWriter_.allocBuffers(writer);
        customWriter_.allocBuffers(writer);

        // iterate over grid
        FVElementGeometry fvGeometry;
        ElementVolumeVariables elemVolVars;
        ElementBoundaryTypes elemBcTypes;

        for (const auto& element : elements(problem_.gridView(), Dune::Partitions::interior))
        {
            fvGeometry.update(problem_.gridView(), element);
            elemBcTypes.update(problem_, element);
            this->problem_.model().setHints(element, elemVolVars);
            elemVolVars.update(problem_,
                               element,
                               fvGeometry,
                               false);
            this->problem_.model().updateCurHints(element, elemVolVars);

            // tell the sub-writers to do what ever they need to with
            // their internal buffers when a given element is seen.
            commonWriter_.processElement(element,
                                         fvGeometry,
                                         elemVolVars,
                                         elemBcTypes);
            massWriter_.processElement(element,
                                       fvGeometry,
                                       elemVolVars,
                                       elemBcTypes);
            energyWriter_.processElement(element,
                                         fvGeometry,
                                         elemVolVars,
                                         elemBcTypes);
            customWriter_.processElement(element,
                                         fvGeometry,
                                         elemVolVars,
                                         elemBcTypes);
        }

        // write everything to the output file
        commonWriter_.commitBuffers(writer);
        massWriter_.commitBuffers(writer);
        energyWriter_.commitBuffers(writer);
        customWriter_.commitBuffers(writer);
    }

private:
    const Problem &problem_;

    MPNCVtkCommonModule commonWriter_;
    MPNCVtkMassModule massWriter_;
    MPNCVtkEnergyModule energyWriter_;
    MPNCVtkCustomModule customWriter_;
};

}

#endif
