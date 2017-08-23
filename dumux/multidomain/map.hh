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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
 /*!
 * \file
 * \ingroup BoundaryCoupling
 * \brief The zwo container types used for storing the information of the two
 *        subdomains
 */

#ifndef DUMUX_COUPLING_STOKES_DARCY_MAP_HH
#define DUMUX_COUPLING_STOKES_DARCY_MAP_HH

#include <dumux/common/propertysystem.hh>

namespace Dumux
{

/*!
 * \ingroup BoundaryCoupling
 * \brief The zwo container types used for storing the information of the two
 *        subdomains
 */
template<typename TypeTag>
class DarcyToStokesMap
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

public:

    void setDarcyElementIndex(const unsigned int idx)
    {
        darcyElementIdx_ = idx;
    }

    void setDarcyScvIdx(const unsigned int idx)
    {
        darcyScvIdx_ = idx;
    }

    void addStokesElementIndex(const unsigned int idx)
    {
        stokesElementIdx_.push_back(idx);
    }

    void addStokesCCDofIndex(const unsigned int idx)
    {
        stokesCCDofIdx_.push_back(idx);
    }

    void addstokesFaceDofIndex(const unsigned int idx)
    {
        stokesFaceDofIdx_.push_back(idx);
    }

    void setCouplingArea(const Scalar area)
    {
        couplingArea_ = area;
    }

    auto& stokesElementIndices() const
    {
        return stokesElementIdx_;
    }

    auto& stokesCCDofIndices() const
    {
        return stokesCCDofIdx_;
    }

    auto& stokesFaceDofIndices() const
    {
        return stokesFaceDofIdx_;
    }

    auto couplingArea() const
    {
        return couplingArea_;
    }

    auto darcyElementIndex() const
    {
        return darcyElementIdx_;
    }

    auto darcyScvIdx() const
    {
        return darcyScvIdx_;
    }

private:
        unsigned int darcyElementIdx_;
        unsigned int darcyScvIdx_;
        std::vector<unsigned int> stokesElementIdx_;
        std::vector<unsigned int> stokesCCDofIdx_;
        std::vector<unsigned int> stokesFaceDofIdx_;
        Scalar couplingArea_;

};


template<typename TypeTag>
class StokesToDarcyMap
{

};


// // helper function to determine the coupled entities
// template<typename GlobalPosition, typename StokesProblem>
// inline std::vector<unsigned int> getCoupledStokesElements(const GlobalPosition& vertexPos, const StokesProblem& problem)
// {
//     const auto &darcyPos =  vertex.geometry().center();
//     const auto &stokesTree = problem.boundingBoxTree();
//     return stokesTree.computeEntityCollisions(darcyPos);
// }

} // end namespace

#endif
