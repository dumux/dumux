// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Parallel
 * \brief Provides data handles for parallel communication which operate on vertices
 * \note This is useful for schemes with degrees of freedom on vertices (box scheme)
 */
#ifndef DUMUX_VERTEX_HANDLES_HH
#define DUMUX_VERTEX_HANDLES_HH

#warning "This header is deprecated and will be removed after release 3.2. Use parallel/vectorcommdatahandle.hh"

#include "vectorcommdatahandle.hh"


namespace Dumux {

    template<class FieldType, class Container, class VertexMapper>
    class [[deprecated("Use VectorCommDataHandleSum<VertexMapper, Container, entityCodim> instead.")]] VertexHandleSum : public VectorCommDataHandleSum<VertexMapper, Container, 0/*dummy*/>
    {
        using ParentType = VectorCommDataHandleSum<VertexMapper, Container, 0/*dummy*/>;
    public:
        VertexHandleSum(Container &container,
                        const VertexMapper &mapper)
            : ParentType(mapper, container)
        { };

        bool contains(int dim, int codim) const
        {
            // only communicate vertices
            return codim == dim;
        }
    };

    template<class FieldType, class Container, class VertexMapper>
    class [[deprecated("Use VectorCommDataHandleMax<VertexMapper, Container, entityCodim> instead.")]] VertexHandleMax : public VectorCommDataHandleMax<VertexMapper, Container, 0/*dummy*/>
    {
        using ParentType = VectorCommDataHandleMax<VertexMapper, Container, 0/*dummy*/>;
    public:
        VertexHandleMax(Container &container,
                        const VertexMapper &mapper)
            : ParentType(mapper, container)
        { };

        bool contains(int dim, int codim) const
        {
            // only communicate vertices
            return codim == dim;
        }
    };

    template<class FieldType, class Container, class VertexMapper>
    class [[deprecated("Use VectorCommDataHandleMin<VertexMapper, Container, entityCodim> instead.")]] VertexHandleMin : public VectorCommDataHandleMin<VertexMapper, Container, 0/*dummy*/>
    {
        using ParentType = VectorCommDataHandleMin<VertexMapper, Container, 0/*dummy*/>;
    public:
        VertexHandleMin(Container &container,
                        const VertexMapper &mapper)
            : ParentType(mapper, container)
        { };

        bool contains(int dim, int codim) const
        {
            // only communicate vertices
            return codim == dim;
        }
    };
} // end namespace Dumux

#endif
