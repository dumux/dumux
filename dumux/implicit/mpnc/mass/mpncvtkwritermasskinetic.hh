/*****************************************************************************
 *   Copyright (C) 2011 by Philipp Nuske                                     *
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
/*!
 * \file
 *
 * \brief VTK writer module for the mass related quantities of the
 *        MpNc model (kinetic mass transfer case).
 */
#ifndef DUMUX_MPNC_VTK_WRITER_MASS_KINETIC_HH
#define DUMUX_MPNC_VTK_WRITER_MASS_KINETIC_HH

#include <dumux/implicit/mpnc/mpncvtkwritermodule.hh>
#include <dumux/implicit/mpnc/mass/mpncvtkwritermass.hh>

namespace Dumux
{

/*!
 * \ingroup MPNCModel
 *
 * \brief VTK writer module for the mass related quantities of the
 *        MpNc model.
 *
 * This is the specialization for the case with kinetic mass transfer.
 */
template<class TypeTag>
class MPNCVtkWriterMass<TypeTag, /* enableKinetic = */ true>
    : public MPNCVtkWriterModule<TypeTag>
{
    typedef MPNCVtkWriterModule<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { dim = GridView::dimension };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx};
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    enum { xEquilOutput   = GET_PROP_VALUE(TypeTag, VtkAddxEquil) };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename ParentType::ScalarVector ScalarVector;
    typedef typename ParentType::PhaseVector PhaseVector;
    typedef typename ParentType::ComponentVector ComponentVector;
    typedef typename ParentType::PhaseComponentMatrix PhaseComponentMatrix;

    bool deltaPOutput_;
public:
    MPNCVtkWriterMass(const Problem &problem)
        : ParentType(problem)
    {
        deltaPOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddDeltaP);
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    template <class MultiWriter>
    void allocBuffers(MultiWriter &writer)
    {
        this->resizePhaseComponentBuffer_(moleFracEquil_);
        if (deltaPOutput_) this->resizeScalarBuffer_(deltaP_, isBox);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const Element &elem,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementBoundaryTypes &elemBcTypes)
    {
        int numLocalVertices = elem.geometry().corners();
        for (int localVertexIdx = 0; localVertexIdx< numLocalVertices; ++localVertexIdx) {
            const unsigned int globalIdx = this->problem_.vertexMapper().map(elem, localVertexIdx, dim);
            const VolumeVariables &volVars = elemVolVars[localVertexIdx];

            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                for (unsigned  int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                    moleFracEquil_[phaseIdx][compIdx][globalIdx]          = volVars.xEquil(phaseIdx, compIdx);
                }
            }
            if (deltaPOutput_) {
                deltaP_[globalIdx] = volVars.fluidState().pressure(nPhaseIdx) - 100000.;
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    template <class MultiWriter>
    void commitBuffers(MultiWriter &writer)
    {
        if(xEquilOutput){
            this->commitPhaseComponentBuffer_(writer, "xEquil_%s^%s", moleFracEquil_);
        }
        if (deltaPOutput_)
            this->commitScalarBuffer_(writer, "pnMinus1e5", deltaP_, isBox);
    }

private:
    PhaseComponentMatrix moleFracEquil_;
    ScalarVector deltaP_;
};

}

#endif
