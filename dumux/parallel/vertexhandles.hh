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
 * \ingroup Assembly
 * \brief Provides data handles for parallel communication which operate on vertices
 * \note This is useful for schemes with degrees of freedom on vertices (box scheme)
 */
#ifndef DUMUX_VERTEX_HANDLES_HH
#define DUMUX_VERTEX_HANDLES_HH

#include <dune/grid/common/datahandleif.hh>


namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief Data handle for parallel communication which sums up all
 *        values are attached to vertices
 */
template <class FieldType, class Container, class VertexMapper>
class VertexHandleSum
    : public Dune::CommDataHandleIF< VertexHandleSum<FieldType, Container, VertexMapper>,
                                     FieldType >
{
public:
    VertexHandleSum(Container &container,
                    const VertexMapper &mapper)
        : mapper_(mapper)
        , container_(container)
    { };

    bool contains(int dim, int codim) const
    {
        // only communicate vertices
        return codim == dim;
    }

    bool fixedsize(int dim, int codim) const
    {
        // for each vertex we communicate a single field vector which
        // has a fixed size
        return true;
    }

    template<class EntityType>
    size_t size (const EntityType &e) const
    {
        // communicate a field type per entity
        return 1;
    }

    template<class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp &buff, const EntityType &e) const
    {
        int vIdx = mapper_.index(e);
        buff.write(container_[vIdx]);
    }

    template<class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        int vIdx = mapper_.index(e);

        FieldType tmp;
        buff.read(tmp);
        container_[vIdx] += tmp;
    }

private:
    const VertexMapper &mapper_;
    Container &container_;
};

/*!
 * \brief Data handle for parallel communication which takes the
 *        maximum of all values that are attached to vertices
 */
template <class FieldType, class Container, class VertexMapper>
class VertexHandleMax
    : public Dune::CommDataHandleIF< VertexHandleMax<FieldType, Container, VertexMapper>,
                                     FieldType >
{
public:
    VertexHandleMax(Container &container,
                    const VertexMapper &mapper)
        : mapper_(mapper)
        , container_(container)
    { };

    bool contains(int dim, int codim) const
    {
        // only communicate vertices
        return codim == dim;
    }

    bool fixedsize(int dim, int codim) const
    {
        // for each vertex we communicate a single field vector which
        // has a fixed size
        return true;
    }

    template<class EntityType>
    size_t size (const EntityType &e) const
    {
        // communicate a field type per entity
        return 1;
    }

    template<class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp &buff, const EntityType &e) const
    {
        int vIdx = mapper_.index(e);
        buff.write(container_[vIdx]);
    }

    template<class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        int vIdx = mapper_.index(e);

        FieldType tmp;
        buff.read(tmp);
        using std::max;
        container_[vIdx] = max(container_[vIdx], tmp);
    }

private:
    const VertexMapper &mapper_;
    Container &container_;
};


/*!
 * \brief Provides data handle for parallel communication which takes
 *        the minimum of all values that are attached to vertices
 */
template <class FieldType, class Container, class VertexMapper>
class VertexHandleMin
    : public Dune::CommDataHandleIF< VertexHandleMin<FieldType, Container, VertexMapper>,
                                     FieldType >
{
public:
    VertexHandleMin(Container &container,
                    const VertexMapper &mapper)
        : mapper_(mapper)
        , container_(container)
    { };

    bool contains(int dim, int codim) const
    {
        // only communicate vertices
        return codim == dim;
    }

    bool fixedsize(int dim, int codim) const
    {
        // for each vertex we communicate a single field vector which
        // has a fixed size
        return true;
    }

    template<class EntityType>
    size_t size (const EntityType &e) const
    {
        // communicate a field type per entity
        return 1;
    }

    template<class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp &buff, const EntityType &e) const
    {
        int vIdx = mapper_.index(e);
        buff.write(container_[vIdx]);
    }

    template<class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        int vIdx = mapper_.index(e);
        FieldType tmp;
        buff.read(tmp);
        using std::min;
        container_[vIdx] = min(container_[vIdx], tmp);
    }

private:
    const VertexMapper &mapper_;
    Container &container_;
};

} // end namespace Dumux

#endif
