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
 * \ingroup Assembly
 * \brief Detects which entries in the Jacobian have to be recomputed
 */
#ifndef DUMUX_PARTIAL_REASSEMBLER_HH
#define DUMUX_PARTIAL_REASSEMBLER_HH

#include <algorithm>
#include <vector>

#include <dune/grid/common/gridenums.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/io/format.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/discretization/method.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#include "entitycolor.hh"

namespace Dumux {

class DefaultPartialReassembler
{
public:
    template<typename... Args>
    DefaultPartialReassembler(Args&&... args)
    { DUNE_THROW(Dune::InvalidStateException, "DefaultPartialReassembler should be never constructed!"); }

    template<typename... Args>
    void report(Args&&... args) {}

    template<typename... Args>
    void resetJacobian(Args&&... args) const {}

    template<typename... Args>
    void computeColors(Args&&... args) {}

    template<typename... Args>
    void resetColors(Args&&... args) {}

    EntityColor dofColor(size_t idx) const
    { return EntityColor::red; }

    EntityColor elementColor(size_t idx) const
    { return EntityColor::red; }

    EntityColor vertexColor(size_t idx) const
    { return EntityColor::red; }
};

//! the partial reassembler engine specialized for discretization methods
template<class Assembler, DiscretizationMethod discMethod>
class PartialReassemblerEngine
{
public:
    PartialReassemblerEngine(const Assembler&)
    { DUNE_THROW(Dune::NotImplemented, "PartialReassembler for this discretization method!"); }

    EntityColor elementColor(size_t idx) const
    { return EntityColor::red; }

    EntityColor dofColor(size_t idx) const
    { return EntityColor::red; }

    template<typename... Args>
    std::size_t computeColors(Args&&... args) { return 0; }

    template<typename... Args>
    void resetJacobian(Args&&... args) const {}

    template<typename... Args>
    void resetColors(Args&&... args) {}
};

/*!
 * \ingroup Assembly
 * \brief The partial reassembler engine specialized for the box method
 */
template<class Assembler>
class PartialReassemblerEngine<Assembler, DiscretizationMethod::box>
{
    using Scalar = typename Assembler::Scalar;
    using GridGeometry = typename Assembler::GridGeometry;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using VertexMapper = typename GridGeometry::VertexMapper;
    static constexpr int dim = GridGeometry::GridView::dimension;

public:
    PartialReassemblerEngine(const Assembler& assembler)
    : elementColor_(assembler.gridGeometry().elementMapper().size(), EntityColor::red)
    , vertexColor_(assembler.gridGeometry().vertexMapper().size(), EntityColor::red)
    {}

    // returns number of green elements
    std::size_t computeColors(const Assembler& assembler,
                              const std::vector<Scalar>& distanceFromLastLinearization,
                              Scalar threshold)
    {
        const auto& gridGeometry = assembler.gridGeometry();
        const auto& gridView = gridGeometry.gridView();
        const auto& elementMapper = gridGeometry.elementMapper();
        const auto& vertexMapper = gridGeometry.vertexMapper();

        // set all vertices to green
        vertexColor_.assign(vertexColor_.size(), EntityColor::green);

        // mark the red vertices
        for (unsigned int i = 0; i < vertexColor_.size(); ++i)
        {
            using std::max;
            if (distanceFromLastLinearization[i] > threshold)
                // mark vertex as red if discrepancy is larger than
                // the relative tolerance
                vertexColor_[i] = EntityColor::red;
        }

        // Mark all red elements
        for (const auto& element : elements(gridView))
        {
            // find out whether the current element features a red vertex
            bool isRed = false;

            int numVertices = element.subEntities(dim);

            for (int i = 0; i < numVertices; ++i) {
                int globalI = vertexMapper.subIndex(element, i, dim);

                if (vertexColor_[globalI] == EntityColor::red) {
                    isRed = true;
                    break;
                }
            }

            int eIdx = elementMapper.index(element);
            // if a vertex is red, the element color is also red, otherwise green
            if (isRed)
                elementColor_[eIdx] = EntityColor::red;
            else
                elementColor_[eIdx] = EntityColor::green;
        }

        // mark orange vertices
        for (const auto& element : elements(gridView))
        {
            int eIdx = elementMapper.index(element);

            // only red elements tint vertices yellow
            if (elementColor_[eIdx] == EntityColor::red)
            {
                int numVertices = element.subEntities(dim);

                for (int i = 0; i < numVertices; ++i) {
                    int globalI = vertexMapper.subIndex(element, i, dim);

                    // red vertices don't become orange
                    if (vertexColor_[globalI] != EntityColor::red)
                        vertexColor_[globalI] = EntityColor::orange;
                }
            }
        }

        // at this point we communicate the yellow vertices to the
        // neighboring processes because a neigbor process may not see
        // the red vertex for yellow border vertices
        VectorCommDataHandleMin<VertexMapper, std::vector<EntityColor>, dim>
            minHandle(vertexMapper, vertexColor_);
        gridView.communicate(minHandle,
                             Dune::InteriorBorder_InteriorBorder_Interface,
                             Dune::ForwardCommunication);

        // mark yellow elements
        for (const auto& element : elements(gridView))
        {
            int eIdx = elementMapper.index(element);

            // only treat non-red elements
            if (elementColor_[eIdx] != EntityColor::red)
            {
                // check whether the element features a orange vertex
                bool isOrange = false;
                int numVertices = element.subEntities(dim);

                for (int i = 0; i < numVertices; ++i) {
                    int globalI = vertexMapper.subIndex(element, i, dim);

                    if (vertexColor_[globalI] == EntityColor::orange) {
                        isOrange = true;
                        break;
                    }
                }

                if (isOrange)
                    elementColor_[eIdx] = EntityColor::yellow;
            }
        }

        // change orange vertices to yellow ones if it has at least
        // one green element as a neighbor
        for (const auto& element : elements(gridView))
        {
            int eIdx = elementMapper.index(element);

            // only green elements are considered
            if (elementColor_[eIdx] == EntityColor::green)
            {
                int numVertices = element.subEntities(dim);

                for (int i = 0; i < numVertices; ++i) {
                    int globalI = vertexMapper.subIndex(element, i, dim);

                    // if a vertex is orange, recolor it to yellow
                    if (vertexColor_[globalI] == EntityColor::orange)
                        vertexColor_[globalI] = EntityColor::yellow;
                }
            }
        }

        // demote the border orange vertices
        VectorCommDataHandleMax<VertexMapper, std::vector<EntityColor>, dim>
            maxHandle(vertexMapper, vertexColor_);
        gridView.communicate(maxHandle,
                             Dune::InteriorBorder_InteriorBorder_Interface,
                             Dune::ForwardCommunication);

        // promote the remaining orange vertices to red
        for (unsigned int i=0; i < vertexColor_.size(); ++i) {
            // if a vertex is green or yellow don't do anything
            if (vertexColor_[i] == EntityColor::green || vertexColor_[i] == EntityColor::yellow)
                continue;

            // set the vertex to red
            vertexColor_[i] = EntityColor::red;
        }

        // count green elements
        return std::count_if(elementColor_.begin(), elementColor_.end(),
                             [](EntityColor c){ return c == EntityColor::green; });
    }

    void resetJacobian(Assembler& assembler) const
    {
        auto& jacobian = assembler.jacobian();

        // loop over all dofs
        for (unsigned int rowIdx = 0; rowIdx < jacobian.N(); ++rowIdx)
        {
            // reset all entries corrosponding to a non-green vertex
            if (vertexColor_[rowIdx] != EntityColor::green)
            {
                // set all matrix entries in the row to 0
                auto colIt = jacobian[rowIdx].begin();
                const auto& colEndIt = jacobian[rowIdx].end();
                for (; colIt != colEndIt; ++colIt) {
                    *colIt = 0.0;
                }
            }
        }
    }

    void resetColors()
    {
        elementColor_.assign(elementColor_.size(), EntityColor::red);
        vertexColor_.assign(vertexColor_.size(), EntityColor::red);
    }

    EntityColor elementColor(size_t idx) const
    { return elementColor_[idx]; }

    EntityColor vertexColor(size_t idx) const
    { return vertexColor_[idx]; }

    EntityColor dofColor(size_t idx) const
    { return vertexColor_[idx]; }

private:
    //! entity colors for partial reassembly
    std::vector<EntityColor> elementColor_;
    std::vector<EntityColor> vertexColor_;
};

/*!
 * \ingroup Assembly
 * \brief The partial reassembler engine specialized for the cellcentered TPFA method
 */
template<class Assembler>
class PartialReassemblerEngine<Assembler, DiscretizationMethod::cctpfa>
{
    using Scalar = typename Assembler::Scalar;
    using GridGeometry = typename Assembler::GridGeometry;
    using JacobianMatrix = typename Assembler::JacobianMatrix;

public:
    PartialReassemblerEngine(const Assembler& assembler)
    : elementColor_(assembler.gridGeometry().elementMapper().size(), EntityColor::red)
    {}

    // returns number of green elements
    std::size_t computeColors(const Assembler& assembler,
                              const std::vector<Scalar>& distanceFromLastLinearization,
                              Scalar threshold)
    {
        const auto& gridGeometry = assembler.gridGeometry();
        const auto& gridView = gridGeometry.gridView();
        const auto& elementMapper = gridGeometry.elementMapper();

        // mark the red elements
        for (const auto& element : elements(gridView))
        {
            int eIdx = elementMapper.index(element);

            if (distanceFromLastLinearization[eIdx] > threshold)
            {
                // mark element as red if discrepancy is larger than
                // the relative tolerance
                elementColor_[eIdx] = EntityColor::red;
            }
            else
            {
                elementColor_[eIdx] = EntityColor::green;
            }
        }

        // mark the neighbors also red
        const auto& connectivityMap = gridGeometry.connectivityMap();
        for (unsigned eIdx = 0; eIdx < elementColor_.size(); ++eIdx)
        {
            if (elementColor_[eIdx] == EntityColor::red)
                continue; // element is red already!

            if (distanceFromLastLinearization[eIdx] > threshold)
            {
                for (const auto& connectedDof : connectivityMap[eIdx])
                    elementColor_[connectedDof.globalJ] = EntityColor::red;
            }
        }

        // count green elements
        return std::count_if(elementColor_.begin(), elementColor_.end(),
                             [](EntityColor c){return c == EntityColor::green;});

    }

    void resetJacobian(Assembler& assembler) const
    {
        auto& jacobian = assembler.jacobian();
        const auto& connectivityMap = assembler.gridGeometry().connectivityMap();

        // loop over all dofs
        for (unsigned int colIdx = 0; colIdx < jacobian.M(); ++colIdx)
        {
            // reset all entries corresponding to a non-green element
            if (elementColor_[colIdx] != EntityColor::green)
            {
                // set all matrix entries in the column to 0
                jacobian[colIdx][colIdx] = 0;
                for (const auto& dataJ : connectivityMap[colIdx])
                    jacobian[dataJ.globalJ][colIdx] = 0;
            }
        }
    }

    void resetColors()
    {
        elementColor_.assign(elementColor_.size(), EntityColor::red);
    }

    EntityColor elementColor(size_t idx) const
    { return elementColor_[idx]; }

    EntityColor dofColor(size_t idx) const
    { return elementColor_[idx]; }

private:
    //! entity colors for partial reassembly
    std::vector<EntityColor> elementColor_;
};

/*!
 * \ingroup Assembly
 * \brief The partial reassembler engine specialized for the cellcentered MPFA method
 */
template<class Assembler>
class PartialReassemblerEngine<Assembler, DiscretizationMethod::ccmpfa>
: public PartialReassemblerEngine<Assembler, DiscretizationMethod::cctpfa>
{
    using ParentType = PartialReassemblerEngine<Assembler, DiscretizationMethod::cctpfa>;
public:
    using ParentType::ParentType;
};

//! helper struct to determine whether the an engine class has vertex colors
struct hasVertexColor
{
    template<class Engine>
    auto operator()(Engine&& e) -> decltype(e.vertexColor(0)) {}
};

/*!
 * \ingroup Assembly
 * \brief detects which entries in the Jacobian have to be recomputed
 * \tparam TypeTag The TypeTag
 */
template<class Assembler>
class PartialReassembler
{
    using Scalar = typename Assembler::Scalar;
    using GridGeometry = typename Assembler::GridGeometry;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using VertexMapper = typename GridGeometry::VertexMapper;

    static constexpr DiscretizationMethod discMethod = GridGeometry::discMethod;
    using Engine = PartialReassemblerEngine<Assembler, discMethod>;

public:

    /*!
     * \brief constructor
     * \param assembler the assembler
     */
    PartialReassembler(const Assembler& assembler)
    : engine_(assembler)
    , greenElems_(0)
    {
        const auto& gridGeometry = assembler.gridGeometry();
        totalElems_ = gridGeometry.elementMapper().size();
        totalElems_ = gridGeometry.gridView().comm().sum(totalElems_);
    }

    /*!
     * \brief Determine the colors of entities for partial reassembly.
     *
     * The following approach is used:
     *
     * - Set all elements to 'green'
     * - Mark all elements as 'red' which exhibit an relative error above
     *   the tolerance
     * - Mark all neighbors of 'red' elements also 'red'
     *
     * \param assembler the assembler
     * \param distanceFromLastLinearization The distance from the last linearization
     * \param threshold Reassemble only if the distance from the last
     *        linearization is above this value.
     */
    void computeColors(const Assembler& assembler,
                       const std::vector<Scalar>& distanceFromLastLinearization,
                       Scalar threshold)
    {
        greenElems_ = engine_.computeColors(assembler, distanceFromLastLinearization, threshold);
    }

    void resetColors()
    {
        engine_.resetColors();
    }

    void resetJacobian(Assembler& assembler) const
    {
        engine_.resetJacobian(assembler);
    }

    /*!
     * \brief called by the assembler after successful assembly
     */
    template <class Communication>
    void report(const Communication& comm, std::ostream& outStream)
    {
        if (comm.size() > 1)
            greenElems_ = comm.sum(greenElems_);

        const auto reassembledElems = totalElems_ - greenElems_;
        const auto width = Fmt::formatted_size("{}", totalElems_);
        outStream << Fmt::format(", reassembled {:{}} ({:3}%) elements",
                                 reassembledElems, width, 100*reassembledElems/totalElems_);
    }

    EntityColor elementColor(size_t idx) const
    { return engine_.elementColor(idx); }

    EntityColor dofColor(size_t idx) const
    { return engine_.dofColor(idx); }

    template<bool enable = decltype(isValid(hasVertexColor()).template check<Engine>())::value,
             typename std::enable_if_t<enable, int> = 0>
    EntityColor vertexColor(size_t idx) const
    { return engine_.vertexColor(idx); }

private:
    Engine engine_;
    size_t totalElems_;
    size_t greenElems_;
};

} // namespace Dumux

#endif
