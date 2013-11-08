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
#ifndef DUMUX_MULTIDOMAIN_BOX_LOCAL_OPERATOR_HH
#define DUMUX_MULTIDOMAIN_BOX_LOCAL_OPERATOR_HH

#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

#include <dumux/implicit/box/boxproperties.hh>

namespace Dumux {

namespace PDELab {

template<class TypeTag>
class MultiDomainBoxLocalOperator
:
public Dune::PDELab::FullVolumePattern,
public Dune::PDELab::LocalOperatorDefaultFlags
{
	// copying the local operator for PDELab is not a good idea
	MultiDomainBoxLocalOperator(const MultiDomainBoxLocalOperator &);

	typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
	typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
	typedef typename Grid::Traits::template Codim<0>::EntityPointer EntityPointer;

	enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};

public:
	// pattern assembly flags
	enum { doPatternVolume = true };

	// residual assembly flags
	enum { doAlphaVolume = true };

	MultiDomainBoxLocalOperator(Model &model)
	: model_(model)
	{}

	/*!
	 * \brief Volume integral depending on test and ansatz functions
	 *
	 * \param eg docme
     * \param lfsu docme
     * \param x docme
	 * \param lfsv docme
	 * \param r docme
	 *
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
	 * \param eg docme
     * \param lfsu docme, is basis
     * \param x docme
	 * \param lfsv docme, is test
	 * \param mat docme
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

#endif
