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
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
#ifndef DUMUX_MORTAR_PROJECTOR_HH
#define DUMUX_MORTAR_PROJECTOR_HH

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/discretization/projection/projector.hh>

namespace Dumux {

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template<class MortarSolution>
class MortarProjectorBase
{
public:
    /*!
     * \brief TODO doc me.
     */
    virtual MortarSolution projectMortarToSubDomain(const MortarSolution& x) const = 0;

    /*!
     * \brief TODO doc me.
     */
    virtual MortarSolution projectSubDomainToMortar(const MortarSolution& x) const = 0;
};

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template<class MortarSolution>
class MortarProjector : public MortarProjectorBase<MortarSolution>
{
    using Scalar = typename MortarSolution::field_type;
    using Projector = Dumux::Projector<Scalar>;

public:
    /*!
     * \brief TODO doc me.
     */
    MortarProjector(Projector&& toSubDomain,
                    Projector&& toMortar)
    : subDomainToMortarProjector_(std::move(toMortar))
    , mortarToSubDomainProjector_(std::move(toSubDomain))
    {}

    /*!
     * \brief TODO doc me.
     */
    virtual MortarSolution projectMortarToSubDomain(const MortarSolution& x) const
    { return mortarToSubDomainProjector_.project(x); }

    /*!
     * \brief TODO doc me.
     */
    virtual MortarSolution projectSubDomainToMortar(const MortarSolution& x) const
    { return subDomainToMortarProjector_.project(x); }

private:
    Projector subDomainToMortarProjector_;
    Projector mortarToSubDomainProjector_;
};

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template<class MortarSolution>
class TransposedMortarProjector : public MortarProjectorBase<MortarSolution>
{
    using Scalar = typename MortarSolution::field_type;
    using Projector = Dumux::Projector<Scalar>;

public:
    /*!
     * \brief TODO doc me.
     */
    template<class Matrix>
    TransposedMortarProjector(Matrix&& massMatrix, Matrix&& projMatrix)
    : B_(projMatrix)
    , mortarToSubDomainProjector_(std::move(massMatrix), std::move(projMatrix))
    {}

    /*!
     * \brief TODO doc me.
     */
    MortarSolution projectMortarToSubDomain(const MortarSolution& x) const override
    { return mortarToSubDomainProjector_.project(x); }

    /*!
     * \brief TODO doc me.
     */
    MortarSolution projectSubDomainToMortar(const MortarSolution& x) const override
    {
        MortarSolution result;
        result.resize(B_.M());

        B_.mtv(x, result);
        return result;
    }

private:
    typename Projector::Matrix B_;
    Projector mortarToSubDomainProjector_;
};

} // end namespace Dumux

#endif
