// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_MPNC_VTK_WRITER_HH
#define DUMUX_MPNC_VTK_WRITER_HH

#include "MpNcproperties.hh"

#include <dumux/io/vtkmultiwriter.hh>

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \brief Writes the VTK output files for a
 * solution of the Mp-Nc model.
 */
template<class TypeTag>
class MPNCVtkWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCVtkCommonModule)) MPNCVtkCommonModule;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCVtkMassModule)) MPNCVtkMassModule;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCVtkEnergyModule)) MPNCVtkEnergyModule;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCVtkCustomModule)) MPNCVtkCustomModule;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum { dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

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
        FVElementGeometry fvElemGeom;
        ElementVolumeVariables elemVolVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator elemIt = problem_.gridView().template begin<0>();
        ElementIterator elemEndIt = problem_.gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            fvElemGeom.update(problem_.gridView(), *elemIt);
            elemBcTypes.update(problem_, *elemIt, fvElemGeom);
            this->problem_.model().setHints(*elemIt, elemVolVars);
            elemVolVars.update(problem_,
                               *elemIt,
                               fvElemGeom,
                               false);
            this->problem_.model().updateCurHints(*elemIt, elemVolVars);

            // tell the sub-writers to do what ever they need to with
            // their internal buffers when a given element is seen.
            commonWriter_.processElement(*elemIt,
                                         fvElemGeom,
                                         elemVolVars,
                                         elemBcTypes);
            massWriter_.processElement(*elemIt,
                                       fvElemGeom,
                                       elemVolVars,
                                       elemBcTypes);
            energyWriter_.processElement(*elemIt,
                                         fvElemGeom,
                                         elemVolVars,
                                         elemBcTypes);
            customWriter_.processElement(*elemIt,
                                         fvElemGeom,
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
