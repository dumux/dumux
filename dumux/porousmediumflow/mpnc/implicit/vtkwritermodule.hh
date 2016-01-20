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
 * \brief A VTK writer module which adheres to the required API but
 *        does nothing.
 */

#ifndef DUMUX_MPNC_VTK_BASE_WRITER_HH
#define DUMUX_MPNC_VTK_BASE_WRITER_HH

#include <array>
#include <cstdio>

#include <dune/istl/bvector.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 *
 * \brief A VTK writer module which adheres to the required API but
 *        does nothing.
 *
 * This class also provides some convenience methods for buffer
 * management.
 */
template<class TypeTag>
class MPNCVtkWriterModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { dim             = GridView::dimension };

public:
    typedef std::vector<Dune::FieldVector<Scalar, 1> > ScalarVector;
    typedef std::array<ScalarVector, numPhases> PhaseVector;
    typedef std::array<ScalarVector, numComponents> ComponentVector;
    typedef std::array<ComponentVector,  numPhases> PhaseComponentMatrix;

    MPNCVtkWriterModule(const Problem &problem)
        : problem_(problem)
    {
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    template <class MultiWriter>
    void allocBuffers(MultiWriter &writer)
    {
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemCurVolVars,
                        const ElementBoundaryTypes &elemBcTypes)
    {
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    template <class MultiWriter>
    void commitBuffers(MultiWriter &writer)
    {
    }

protected:
    /*!
     * \brief Allocate the space for a buffer storing a scalar quantity
     */
    void resizeScalarBuffer_(ScalarVector &buffer,
                             bool vertexCentered = true)
    {
        Scalar n; // numVertices for vertexCentereed, numVolumes for volume centered
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        buffer.resize(n);
        std::fill(buffer.begin(), buffer.end(), 0.0);
    }

    /*!
     * \brief Allocate the space for a buffer storing a phase-specific
     *        quantity
     */
    void resizePhaseBuffer_(PhaseVector &buffer,
                            bool vertexCentered = true)
    {
        Scalar n; // numVertices for vertexCentereed, numVolumes for volume centered
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            buffer[phaseIdx].resize(n);
            std::fill(buffer[phaseIdx].begin(), buffer[phaseIdx].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a component
     *        specific quantity
     */
    void resizeComponentBuffer_(ComponentVector &buffer,
                                bool vertexCentered = true)
    {
        Scalar n;// numVertices for vertexCentereed, numVolumes for volume centered
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            buffer[compIdx].resize(n);
            std::fill(buffer[compIdx].begin(), buffer[compIdx].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a phase and
     *        component specific buffer
     */
    void resizePhaseComponentBuffer_(PhaseComponentMatrix &buffer,
                                     bool vertexCentered = true)
    {
        Scalar n;// numVertices for vertexCentereed, numVolumes for volume centered
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                buffer[phaseIdx][compIdx].resize(n);
                std::fill(buffer[phaseIdx][compIdx].begin(), buffer[phaseIdx][compIdx].end(), 0.0);
            }
        }
    }

    /*!
     * \brief Add a phase-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitScalarBuffer_(MultiWriter &writer,
                             const char *name,
                             ScalarVector &buffer,
                             bool vertexCentered = true)
    {
        if (vertexCentered)
            writer.attachVertexData(buffer, name, 1);
        else
            writer.attachCellData(buffer, name, 1);
    }

    /*!
     * \brief Add a phase-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitPhaseBuffer_(MultiWriter &writer,
                            const char *pattern,
                            PhaseVector &buffer,
                            bool vertexCentered = true)
    {
        char name[512];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            snprintf(name, 512, pattern, FluidSystem::phaseName(phaseIdx));

            if (vertexCentered)
                writer.attachVertexData(buffer[phaseIdx], name, 1);
            else
                writer.attachCellData(buffer[phaseIdx], name, 1);
        }
    }

    /*!
     * \brief Add a component-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitComponentBuffer_(MultiWriter &writer,
                                const char *pattern,
                                ComponentVector &buffer,
                                bool vertexCentered = true)
    {
        char name[512];
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            snprintf(name, 512, pattern, FluidSystem::componentName(compIdx));

            if (vertexCentered)
                writer.attachVertexData(buffer[compIdx], name, 1);
            else
                writer.attachCellData(buffer[compIdx], name, 1);
        }
    }

    /*!
     * \brief Add a phase and component specific quantities to the output.
     */
    template <class MultiWriter>
    void commitPhaseComponentBuffer_(MultiWriter &writer,
                                     const char *pattern,
                                     PhaseComponentMatrix &buffer,
                                     bool vertexCentered = true)
    {
        char name[512];
        for (int phaseIdx= 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                snprintf(name, 512, pattern,
                         FluidSystem::phaseName(phaseIdx),
                         FluidSystem::componentName(compIdx));

                if (vertexCentered)
                    writer.attachVertexData(buffer[phaseIdx][compIdx], name, 1);
                else
                    writer.attachCellData(buffer[phaseIdx][compIdx], name, 1);
            }
        }
    }

    const Problem &problem_;
};

}

#endif
