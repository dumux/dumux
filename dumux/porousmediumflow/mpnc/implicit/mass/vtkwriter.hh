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
 * \brief VTK writer module for the mass related quantities of the
 *        MpNc model.
 */
#ifndef DUMUX_MPNC_VTK_WRITER_MASS_HH
#define DUMUX_MPNC_VTK_WRITER_MASS_HH

#include "../vtkwritermodule.hh"

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 *
 * \brief VTK writer module for the mass related quantities of the
 *        MpNc model.
 *
 * This is the specialization for the case _without_ kinetic mass
 * transfer between phases.
 */
template<class TypeTag, bool enableKinetic /* = false */>
class MPNCVtkWriterMass : public MPNCVtkWriterModule<TypeTag>
{
    static_assert(!enableKinetic,
                  "No kinetic mass transfer module included, "
                  "but kinetic mass transfer enabled.");

    typedef MPNCVtkWriterModule<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename ParentType::ComponentVector ComponentVector;
    bool fugacityOutput_;

public:
    MPNCVtkWriterMass(const Problem &problem)
        : ParentType(problem)
    {
        fugacityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddFugacities);
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    template <class MultiWriter>
    void allocBuffers(MultiWriter &writer)
    {
        if (fugacityOutput_) this->resizeComponentBuffer_(fugacity_, isBox);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     *
     *        \param element The finite element
     *        \param fvGeometry The finite-volume geometry in the fully implicit scheme
     *        \param elemVolVars The volume variables of the current element
     *        \param elemBcTypes The types of the boundary conditions for all vertices of the element
     */
    void processElement(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementBoundaryTypes &elemBcTypes)
    {
        for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx) {
            const unsigned int dofIdxGlobal = this->problem_.model().dofMapper().subIndex(element, scvIdx, dofCodim);
            const VolumeVariables &volVars = elemVolVars[scvIdx];

            if (fugacityOutput_) {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    fugacity_[compIdx][dofIdxGlobal] = volVars.fugacity(compIdx);
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    template <class MultiWriter>
    void commitBuffers(MultiWriter &writer)
    {
        if (fugacityOutput_)
            this->commitComponentBuffer_(writer, "f_%s", fugacity_, isBox);
    }

private:
    ComponentVector fugacity_;
};

}

#endif
