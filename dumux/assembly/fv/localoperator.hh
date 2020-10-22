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
 * \brief An element-wise local operator for finite-volume schemes.
 */
#ifndef DUMUX_FV_LOCAL_OPERATOR_HH
#define DUMUX_FV_LOCAL_OPERATOR_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Assembly
 * \brief The element-wise local operator for finite volume schemes.
 *        This allows for element-wise evaluation of individual terms
 *        of the equations to be solved.
 * \tparam OP The model-specific operators
 */
template<class OP>
class FVLocalOperator
{
    using LC = typename OP::LocalContext;
    using ElementVariables = typename LC::ElementVariables;
    using GridVariables = typename ElementVariables::GridVariables;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using FVElementGeometry = typename LC::ElementGridGeometry;
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int numEq = NumEqVectorTraits<PrimaryVariables>::numEq();
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;

public:
    //! export the expected local context type
    using LocalContext = LC;

    //! export the underlying operators
    using Operators = OP;

    //! vector type storing the operator values on all dofs of an element
    //! TODO: Use ReservedBlockVector
    using ElementResidualVector = Dune::BlockVector<NumEqVector>;

    //! Constructor from a local context
    explicit FVLocalOperator(const LocalContext& context)
    : context_(context)
    {}

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the terms of the local residual that do not appear in
     *        time derivatives. These are the sources and the fluxes.
     */
    ElementResidualVector evalFluxesAndSources() const
    {
        auto result = getEmptyResidual();
        const auto& problem = context_.elementVariables().gridVariables().gridVolVars().problem();
        const auto& fvGeometry = context_.elementGridGeometry();

        // source term
        for (const auto& scv : scvs(fvGeometry))
            result[scv.localDofIndex()] -= Operators::source(problem, context_, scv);

        // flux term
        for (const auto& scvf : scvfs(fvGeometry))
            addFlux_(result, scvf);

        return result;
    }

    /*!
     * \brief Compute the storage term, i.e. the term appearing in the time derivative.
     */
    ElementResidualVector evalStorage() const
    {
        const auto& problem = context_.elementVariables().gridVariables().gridVolVars().problem();
        const auto& fvGeometry = context_.elementGridGeometry();

        // TODO: Until now, FVLocalResidual defined storage as the entire
        //       time derivative. Now it means the term above the time derivative.
        //       We should think about the correct naming here...
        // TODO: Should storage() NOT multiply with volume?? That was different until
        //       now but flux() also returns the flux multiplied with area so this should
        //       be more consistent
        auto result = getEmptyResidual();
        for (const auto& scv : scvs(fvGeometry))
            result[scv.localDofIndex()] += Operators::storage(problem, context_, scv);

        return result;
    }

    ElementResidualVector getEmptyResidual() const
    {
        ElementResidualVector res(context_.elementGridGeometry().numScv());
        res = 0.0;
        return res;
    }

    // \}

    /*!
     * \name Interfaces for analytic Jacobian computation
     */
    // \{

    //! \todo TODO: Add interfaces. Or, should this be here at all!?

    //\}

    // \}

protected:

    //! compute and add the flux across the given face to the container (cc schemes)
    template<bool b = isBox, std::enable_if_t<!b, int> = 0>
    void addFlux_(ElementResidualVector& r, const SubControlVolumeFace& scvf) const
    {
        const auto& evv = context_.elementVariables().elemVolVars();
        const auto& problem = evv.gridVolVars().problem();

        // TODO: Modify problem interfaces to receive context
        const auto& fvGeometry = context_.elementGridGeometry();
        const auto& element = fvGeometry.element();
        const auto& efvc = context_.elementVariables().elemFluxVarsCache();

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto localDofIdx = insideScv.localDofIndex();

        if (!scvf.boundary())
            r[localDofIdx] += Operators::flux(problem, context_, scvf);
        else
        {
            const auto& bcTypes = problem.boundaryTypes(element, scvf);

            // Dirichlet boundaries
            if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
                r[localDofIdx] += Operators::flux(problem, context_, scvf);

            // Neumann and Robin ("solution dependent Neumann") boundary conditions
            else if (bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
            {
                auto neumannFluxes = problem.neumann(element, fvGeometry, evv, efvc, scvf);

                // multiply neumann fluxes with the area and the extrusion factor
                neumannFluxes *= Extrusion::area(scvf)*evv[insideScv].extrusionFactor();
                r[localDofIdx] += neumannFluxes;
            }

            else
                DUNE_THROW(Dune::NotImplemented, "Mixed boundary conditions for cell-centered schemes. " <<
                                                 "Use pure boundary conditions by converting Dirichlet BCs to Robin BCs");
        }
    }

    //! compute and add the flux across the given face to the container (box scheme)
    template<bool b = isBox, std::enable_if_t<b, int> = 0>
    void addFlux_(ElementResidualVector& r, const SubControlVolumeFace& scvf) const
    {
        const auto& evv = context_.elementVariables().elemVolVars();
        const auto& problem = evv.gridVolVars().problem();

        // TODO: Modify problem interfaces to receive context
        const auto& fvGeometry = context_.elementGridGeometry();
        const auto& element = fvGeometry.element();
        const auto& efvc = context_.elementVariables().elemFluxVarsCache();

        // inner faces
        if (!scvf.boundary())
        {
            const auto flux = Operators::flux(problem, context_, scvf);
            r[fvGeometry.scv(scvf.insideScvIdx()).localDofIndex()] += flux;

            if (scvf.numOutsideScvs() > 0)
                r[fvGeometry.scv(scvf.outsideScvIdx()).localDofIndex()] -= flux;
        }

        // boundary faces
        else
        {
            const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& bcTypes = problem.boundaryTypes(element, scv);

            // Treat Neumann and Robin ("solution dependent Neumann") boundary conditions.
            if (bcTypes.hasNeumann())
            {
                const auto neumannFluxes = problem.neumann(element, fvGeometry, evv, efvc, scvf);
                const auto area = Extrusion::area(scvf)*evv[scv].extrusionFactor();

                // only add fluxes to equations for which Neumann is set
                for (int eqIdx = 0; eqIdx < NumEqVector::size(); ++eqIdx)
                    if (bcTypes.isNeumann(eqIdx))
                        r[scv.localDofIndex()][eqIdx] += neumannFluxes[eqIdx]*area;
            }
        }
    }

private:
    const LocalContext& context_;
};

} // end namespace Dumux::Experimental

#endif
