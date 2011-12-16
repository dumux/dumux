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
#ifndef DUMUX_MPNC_VTK_BASE_WRITER_HH
#define DUMUX_MPNC_VTK_BASE_WRITER_HH

#include "MpNcproperties.hh"

#include <dumux/io/vtkmultiwriter.hh>
#include <dune/istl/bvector.hh>

#include <tr1/array>

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
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;

    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim             = GridView::dimension };

public:
    typedef std::vector<Dune::FieldVector<Scalar, 1> > ScalarBuffer;
    typedef std::tr1::array<ScalarBuffer, numPhases> PhaseBuffer;
    typedef std::tr1::array<ScalarBuffer, numComponents> ComponentBuffer;
    typedef std::tr1::array<ComponentBuffer,  numPhases> PhaseComponentBuffer;

//    typedef Dune::FieldVector<Scalar, dim> VelocityVector;
//    typedef Dune::BlockVector<VelocityVector> VelocityField;
//    typedef std::tr1::array<VelocityField, numPhases> PhaseVelocityField;

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
    void processElement(const Element &elem,
                        const FVElementGeometry &fvElemGeom,
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
    void resizeScalarBuffer_(ScalarBuffer &buffer,
                             bool vertexCentered = true)
    {
        Scalar n;
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
    void resizePhaseBuffer_(PhaseBuffer &buffer,
                            bool vertexCentered = true)
    {
        Scalar n;
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int i = 0; i < numPhases; ++i) {
            buffer[i].resize(n);
            std::fill(buffer[i].begin(), buffer[i].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a component
     *        specific quantity
     */
    void resizeComponentBuffer_(ComponentBuffer &buffer,
                                bool vertexCentered = true)
    {
        Scalar n;
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int i = 0; i < numComponents; ++i) {
            buffer[i].resize(n);
            std::fill(buffer[i].begin(), buffer[i].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a phase and
     *        component specific buffer
     */
    void resizePhaseComponentBuffer_(PhaseComponentBuffer &buffer,
                                     bool vertexCentered = true)
    {
        Scalar n;
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int i = 0; i < numPhases; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                buffer[i][j].resize(n);
                std::fill(buffer[i][j].begin(), buffer[i][j].end(), 0.0);
            }
        }
    }

    /*!
     * \brief Add a phase-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitScalarBuffer_(MultiWriter &writer,
                             const char *name,
                             ScalarBuffer &buffer,
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
                            PhaseBuffer &buffer,
                            bool vertexCentered = true)
    {
        for (int i = 0; i < numPhases; ++i) {
            std::string phaseName = FluidSystem::phaseName(i);
            std::string name;
            name = (boost::format(pattern)%phaseName).str();

            if (vertexCentered)
                writer.attachVertexData(buffer[i], name.c_str(), 1);
            else
                writer.attachCellData(buffer[i], name.c_str(), 1);
        }
    }

    /*!
     * \brief Add a component-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitComponentBuffer_(MultiWriter &writer,
                                const char *pattern,
                                ComponentBuffer &buffer,
                                bool vertexCentered = true)
    {
        for (int i = 0; i < numComponents; ++i) {
            std::string compName = FluidSystem::componentName(i);
            std::string name;
            name = (boost::format(pattern)%compName).str();

            if (vertexCentered)
                writer.attachVertexData(buffer[i], name.c_str(), 1);
            else
                writer.attachCellData(buffer[i], name.c_str(), 1);
        }
    }

    /*!
     * \brief Add a phase and component specific quantities to the output.
     */
    template <class MultiWriter>
    void commitPhaseComponentBuffer_(MultiWriter &writer,
                                     const char *pattern,
                                     PhaseComponentBuffer &buffer,
                                     bool vertexCentered = true)
    {
        for (int i= 0; i < numPhases; ++i) {
            std::string phaseName = FluidSystem::phaseName(i);
            for (int j = 0; j < numComponents; ++j) {
                std::string compName = FluidSystem::componentName(j);
                std::string name;
                name = (boost::format(pattern)%phaseName%compName).str();

                if (vertexCentered)
                    writer.attachVertexData(buffer[i][j], name.c_str(), 1);
                else
                    writer.attachCellData(buffer[i][j], name.c_str(), 1);
            }
        }
    }

    const Problem &problem_;
};

}

#endif
