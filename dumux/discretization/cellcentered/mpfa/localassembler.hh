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
 * \ingroup CCMpfaDiscretization
 * \brief Defines the general interface of classes used for the assembly
 *        of the local systems of equations involved in the transmissibility
 *        computaion in mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_LOCAL_ASSEMBLER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_LOCAL_ASSEMBLER_HH

#include <dune/common/exceptions.hh>

#include <dumux/discretization/cellcentered/mpfa/methods.hh>

namespace Dumux
{
//! Forward declaration of the implementation
template< class P, class EG, class EV, MpfaMethods M > class InteractionVolumeAssemblerImpl;

//! Alias to select the right implementation.
template< class P, class EG, class EV, MpfaMethods M >
using InteractionVolumeAssembler = InteractionVolumeAssemblerImpl< P, EG, EV, M >;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Defines the general interface of the local assembler
 *        classes for the assembly of the interaction volume-local
 *        transmissibility matrix. Specializations have to be provided
 *        for the available interaction volume implementations. these
 *        should derive from this base clases.
 *
 * \tparam P The problem type
 * \tparam EG The element finite volume geometry
 * \tparam EV The element volume variables type
 */
template< class P, class EG, class EV >
class InteractionVolumeAssemblerBase
{
    using Problem = P;
    using FVElementGeometry = EG;
    using ElementVolumeVariables = EV;

  public:
    /*!
     * \brief The constructor.
     *        Sets pointers to the objects required for a subsequent call to assemble().
     *
     * \param problem The problem to be solved (boundary/initial conditions etc.)
     * \param fvGeometry The local view on the finite volume grid geometry
     * \param elemVolVars The local view on the primary/secondary variables
     */
    InteractionVolumeAssemblerBase(const Problem& problem,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars)
    {
        problemPtr_ = &problem;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
    }

    // return functions to the local views & problem
    const Problem& problem() const { return *problemPtr_; }
    const FVElementGeometry& fvGeometry() const { return *fvGeometryPtr_; }
    const ElementVolumeVariables& elemVolVars() const { return *elemVolVarsPtr_; }

    /*!
     * \brief Assembles the matrices involved in the flux
     *        expressions and the local system of equations
     *        within an mpfa interaction volume.
     *
     * \tparam IV The interaction volume type implementation
     * \tparam TensorFunc Lambda to obtain the tensor w.r.t.
     *                    which the local system is to be solved
     *
     * \param handle The data handle in which the matrices are stored
     * \param iv The interaction volume
     * \param getT Lambda to evaluate the scv-wise tensors
     */
    template< class DataHandle, class IV, class TensorFunc >
    void assembleMatrices(DataHandle& handle, IV& iv, const TensorFunc& getT)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assembleMatrices() function");
    }

    /*!
     * \brief Assembles the vector of primary (cell) unknowns and (maybe)
     *        Dirichlet boundary conditions within an interaction volume.
     *
     * \tparam IV The interaction volume type implementation
     * \tparam GetU Lambda to obtain the cell unknowns from grid indices
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The mpfa-o interaction volume
     * \param getU Lambda to obtain the desired cell/Dirichlet value from vol vars
     */
    template< class DataHandle, class IV, class GetU >
    void assembleU(DataHandle& handle, const IV& iv, const GetU& getU)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assemble() function for the cell/Dirichlet unknowns");
    }

    /*!
     * \brief Assembles the gravitational flux contributions on the scvfs within an
     *        interaction volume.
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The mpfa-o interaction volume
     * \param getRho Lambda to obtain the density from volume variables
     */
    template< class DataHandle, class IV, class GetRho >
    void assembleGravity(DataHandle& handle, const IV& iv, const GetRho& getRho)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assembleGravity() function");
    }

  private:
    // pointers to the data required for assembly
    const Problem* problemPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;
};

} // end namespace Dumux

//! include all specializations for different mpfa schemes
#include <dumux/discretization/cellcentered/mpfa/omethod/localassembler.hh>

#endif
