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
 * \brief Local operator base class for multidomain problems
 */

#ifndef DUMUX_MULTIDOMAIN_LOCAL_OPERATOR_HH
#define DUMUX_MULTIDOMAIN_LOCAL_OPERATOR_HH

#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

#include <dumux/implicit/box/boxproperties.hh>

namespace Dumux {

namespace PDELab {

/*!
 * \ingroup MultidomainModel
 * \brief Local operator base class for multidomain problems
 */
template<class TypeTag>
class MultiDomainLocalOperator
 : public Dune::PDELab::FullVolumePattern,
   public Dune::PDELab::LocalOperatorDefaultFlags
{
    // copying the local operator for PDELab is not a good idea
    MultiDomainLocalOperator(const MultiDomainLocalOperator &);

    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename Grid::Traits::template Codim<0>::EntityPointer EntityPointer;

    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};

public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };

    //! \brief The constructor
    MultiDomainLocalOperator(Model &model)
    : model_(model)
    {}

    /*!
     * \brief Volume integral depending on test and ansatz functions
     *
     * \tparam EG Element geometry
     * \tparam LFSU Local function space for ansatz functions
     * \tparam X Coefficient vector
     * \tparam LFSV Local function space for test functions
     * \tparam R Residual vector
     *
     * \param eg Element geometry
     * \param lfsu Local functions space for ansatz functions
     * \param x Coefficient vector
     * \param lfsv Local function space for test functions
     * \param r Residual vector
     */
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
            const LFSV& lfsv, R& r) const
    {
        typedef typename LFSU::Traits::SizeType size_type;

        const EntityPointer elementPtr = model_.gridView().grid().subDomainEntityPointer(eg.entity());
        model_.localResidual().eval(*elementPtr);

        int numVertices = x.size()/numEq;
        for (size_type comp = 0; comp < r.size(); comp++)
            r.accumulate(lfsv, comp, model_.localResidual().residual(comp%numVertices)[comp/numVertices]);
    }

    /*!
     * \brief Jacobian of volume term
     *
     * \tparam EG Element geometry
     * \tparam LFSU Local function space for ansatz functions
     * \tparam X Coefficient vector
     * \tparam LFSV Local function space for test functions
     * \tparam M Matrix
     *
     * \param eg Element geometry
     * \param lfsu Local functions space for ansatz functions
     * \param x Coefficient vector
     * \param lfsv Local function space for test functions
     * \param mat Matrix
     *
     */
    template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume (const EG& eg,
            const LFSU& lfsu,
            const X& x,
            const LFSV& lfsv,
            M& mat) const
    {
        typedef typename LFSU::Traits::SizeType size_typeU;
        typedef typename LFSV::Traits::SizeType size_typeV;

        const EntityPointer elementPtr = model_.gridView().grid().subDomainEntityPointer(eg.entity());
        model_.localJacobian().assemble(*elementPtr);

        int numVertices = x.size()/numEq;
        for (size_typeV j=0; j<lfsv.size(); j++) {
            for (size_typeU i=0; i<lfsu.size(); i++) {
                mat.accumulate(lfsv, i, lfsu, j, (model_.localJacobian().mat(i%numVertices,j%numVertices))[i/numVertices][j/numVertices]);
            }
        }
    }

private:
    Model& model_;
};

} // namespace PDELab
} // namespace Dumux

#endif // DUMUX_MULTIDOMAIN_LOCAL_OPERATOR_HH
