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

namespace Dumux
{
//! Forward declaration of the implementation
template< class InteractionVolume > class InteractionVolumeAssemblerImpl;

//! Alias to select the right implementation.
template< class InteractionVolume >
using InteractionVolumeAssembler = InteractionVolumeAssemblerImpl< InteractionVolume >;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Defines the general interface of the local assembler
 *        classes for the assembly of the interaction volume-local
 *        transmissibility matrix. Specializations have to be provided
 *        for the available interaction volume implementations. these
 *        should derive from this base clases.
 *
 * \tparam Interaction Volume The interaction volume implementation
 *                            used by the mpfa scheme
 */
template<class InteractionVolume>
class InteractionVolumeAssemblerBase
{
    using Problem = typename InteractionVolume::Problem;
    using FVElementGeometry = typename InteractionVolume::FVElementGeometry;
    using ElementVolumeVariables = typename InteractionVolume::ElementVolumeVariables;

    using Traits = typename InteractionVolume::Traits;
    using Matrix = typename Traits::Matrix;
    using Vector = typename Traits::Vector;

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
     * \brief General interface of a function assembling the
     *        interaction volume-local transmissibility matrix.
     *
     * \tparam GetTensorFunction Lambda to obtain the tensor w.r.t.
     *                           which the local system is to be solved
     *
     * \param T The transmissibility matrix to be assembled
     * \param iv The interaction volume
     * \param getTensor Lambda to evaluate the scv-wise tensors
     */
    template< class GetTensorFunction >
    void assemble(Matrix& T, InteractionVolume& iv, const GetTensorFunction& getTensor)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assemble() function");
    }

    /*!
     * \brief General interface of a function assembling the interaction
     *        volume-local transmissibilities matrix for surface grids. The
     *        transmissibility associated with "outside" faces are stored
     *        in a separate container.
     *
     * \param outsideTij tij on "outside" faces to be assembled
     * \param T The transmissibility matrix tij to be assembled
     * \param iv The mpfa-o interaction volume
     * \param getTensor Lambda to evaluate the scv-wise tensors
     */
    template< class OutsideTijContainer, class GetTensorFunction >
    void assemble(OutsideTijContainer& outsideTij, Matrix& T, InteractionVolume& iv, const GetTensorFunction& getTensor)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assemble() function to be used on surface grids");
    }

    /*!
     * \brief General interface of a function assembling the interaction
     *        volume-local transmissibility matrix in the case that gravity
     *        is to be considered in the local system of equations.
     *
     * \tparam GravityContainer The type of container used to store the
     *                          gravitational acceleration per scvf & phase
     * \tparam GetTensorFunction Lambda to obtain the tensor w.r.t.
     *                           which the local system is to be solved
     *
     * \param T The transmissibility matrix to be assembled
     * \param g Container to assemble gravity per scvf & phase
     * \param CA Matrix to store matrix product C*A^-1
     * \param iv The interaction volume
     * \param getTensor Lambda to evaluate the scv-wise tensors
     */
    template< class GravityContainer, class GetTensorFunction >
    void assembleWithGravity(Matrix& T,
                             GravityContainer& g,
                             Matrix& CA,
                             InteractionVolume& iv,
                             const GetTensorFunction& getTensor)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assembleWithGravity() function");
    }

    /*!
     * \brief General interface of a function assembling the interaction
     *        volume-local transmissibility matrix in the case that gravity
     *        is to be considered in the local system of equations. This
     *        specialization is to be used on surface grids, where the
     *        gravitational flux contributions on "outside" faces are stored
     *        in a separate container.
     *
     * \param outsideTij tij on "outside" faces to be assembled
     * \param T The transmissibility matrix to be assembled
     * \param outsideG Container to assemble gravity on "outside" faces
     * \param g Container to assemble gravity per scvf & phase
     * \param CA Matrix to store matrix product C*A^-1
     * \param A Matrix to store the inverse A^-1
     * \param iv The interaction volume
     * \param getTensor Lambda to evaluate the scv-wise tensors
     */
    template< class GravityContainer,
              class OutsideGravityContainer,
              class OutsideTijContainer,
              class GetTensorFunction >
    void assembleWithGravity(OutsideTijContainer& outsideTij,
                             Matrix& T,
                             OutsideGravityContainer& outsideG,
                             GravityContainer& g,
                             Matrix& CA,
                             Matrix& A,
                             InteractionVolume& iv,
                             const GetTensorFunction& getTensor)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assembleWithGravity() function to be used on surface grids");
    }

    /*!
     * \brief General interface for the assembly of the vector of
     *        primary (cell) unknowns and (maybe) Dirichlet boundary
     *        conditions within an interaction volume.
     *
     * \tparam GetU Lambda to obtain the cell unknowns from grid indices
     *
     * \param u The vector to be filled with the cell unknowns
     * \param iv The interaction volume
     * \param getU Lambda to obtain the desired cell value from grid indices
     */
    template< class GetU >
    void assemble(Vector& u, const InteractionVolume& iv, const GetU& getU)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assemble() function for the cell unknowns");
    }

    /*!
     * \brief General interface for the assembly of the gravitational
     *        flux contributions on the scvfs within an interaction volume.
     *
     * \note  For each face, the gravity term in the form of \f$\rho \mathbf{n K g}\f$ is
     *        evaluated. Thus, make sure to only call this with a lambda that returns the
     *        hydraulic conductivity.
     *
     * \param g Container to assemble gravity per scvf & phase
     * \param iv The interaction volume
     * \param CA Projection matrix transforming the gravity terms in the local system of
     *        equations to the entire set of faces within the interaction volume
     * \param getTensor Lambda to evaluate scv-wise hydraulic conductivities
     */
    template< class GravityContainer, class GetTensorFunction >
    void assembleGravity(GravityContainer& g,
                         const InteractionVolume& iv,
                         const Matrix& CA,
                         const GetTensorFunction& getTensor)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assembleGravity() function");
    }

    /*!
     * \brief General interface for the assembly of the gravitational
     *        flux contributions on the scvfs within an interaction volume.
     *        This specialization is to be used on surface grids, where the gravitational
     *        flux contributions on "outside" faces are stored in a separate container.
     *
     * \note  For each face, the gravity term in the form of \f$\rho \mathbf{n K g}\f$ is
     *        evaluated. Thus, make sure to only call this with a lambda that returns the
     *        hydraulic conductivity.
     *
     * \param g Container to store gravity per scvf & phase
     * \param outsideG Container to store gravity per "outside" scvf & phase
     * \param iv The mpfa-o interaction volume
     * \param CA Projection matrix transforming the gravity terms in the local system of
     *        equations to the entire set of faces within the interaction volume
     * \param A Matrix needed for the "reconstruction" of face unknowns as a function of gravity
     * \param getTensor Lambda to evaluate scv-wise hydraulic conductivities
     */
    template< class GravityContainer,
              class OutsideGravityContainer,
              class GetTensorFunction >
    void assembleGravity(GravityContainer& g,
                         OutsideGravityContainer& outsideG,
                         const InteractionVolume& iv,
                         const Matrix& CA,
                         const Matrix& A,
                         const GetTensorFunction& getTensor)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assembleGravity() function to be used on surface grids");
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
