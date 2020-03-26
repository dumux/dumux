// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \ingroup SequentialTwoPModel
 * \brief Provides methods for transmissibility calculation 2-d.
 */
#ifndef DUMUX_FVMPFAL2D_TRANSMISSIBILITYCALCULATOR_HH
#define DUMUX_FVMPFAL2D_TRANSMISSIBILITYCALCULATOR_HH

// dumux environment
#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/linteractionvolume.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Provides methods for transmissibility calculation in 2-d.
 *
 *  The transmissibilities are calculated using the MPFA L-method.
 *
 *  Aavatsmark et al. A compact multipoint flux calculation method with improved robustness.
 *  Numerical Methods for Partial Differential Equations 24. 2008
 */
template<class TypeTag>
class FvMpfaL2dTransmissibilityCalculator
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using DimVector = Dune::FieldVector<Scalar, dim>;

    using InteractionVolume = FVMPFALInteractionVolume<TypeTag>;

public:
    using TransmissibilityType = Dune::FieldMatrix<Scalar, dim, 2*dim - dim + 1>;//!< Type of the transmissibility matrix

    //! return values for the transmissibility functions
    enum
    {
        leftTriangle = -1,//!< Left L-shape
        noTransmissibility = 0,//!< No transmissibility calculated
        rightTriangle = 1 //!< Right L-shape
    };

    //! Calculates tranmissibility matrix
    int calculateTransmissibility(
            TransmissibilityType& transmissibility,
            InteractionVolume& interactionVolume,
            std::vector<DimVector >& lambda,
            int idx1, int idx2, int idx3, int idx4);

    //! Calculates tranmissibility matrix of left L-shape
    int calculateLeftHNTransmissibility(TransmissibilityType& transmissibility,
    InteractionVolume& interactionVolume,
    std::vector<DimVector >& lambda,
    int idx1, int idx2, int idx3);

    //! Calculates tranmissibility matrix of right L-shape
    int calculateRightHNTransmissibility(TransmissibilityType& transmissibility,
    InteractionVolume& interactionVolume,
    std::vector<DimVector >& lambda,
    int idx1, int idx2, int idx3);

    /*!
     * \brief Constructs a FvMpfaL2dTransmissibilityCalculator object
     *
     * \param problem A problem class object
     */
    FvMpfaL2dTransmissibilityCalculator(Problem& problem) :
        problem_(problem), R_(0)
    {
        if (dim != 2)
        {
            DUNE_THROW(Dune::NotImplemented, "Dimension not supported!");
        }

        // evaluate matrix R
        if (dim == 2)
        {
            R_[0][1] = 1;
            R_[1][0] = -1;
        }
    }

private:
    Problem& problem_;
    DimMatrix R_;
};

/*!
 * \ingroup SequentialTwoPModel
 * \brief Calculates tranmissibility matrix
 *
 * Calculates tranmissibility matrix of an L-shape for a certain flux face.
 * Automatically selects one of the two possible L-shape (left, or right).
 *
 * \param transmissibility Matrix for the resulting transmissibility
 * \param interactionVolume The interaction volume object (includes geometric information)
 * \param lambda Mobilities of cells 1-4
 * \param idx1 Index of cell 1 of the L-stencil
 * \param idx2 Index of cell 2 of the L-stencil
 * \param idx3 Index of cell 3 of the L-stencil
 * \param idx4 Index of cell 4 of the L-stencil
 */
template<class TypeTag>
int FvMpfaL2dTransmissibilityCalculator<TypeTag>::calculateTransmissibility(
        TransmissibilityType& transmissibility,
        InteractionVolume& interactionVolume,
        std::vector<DimVector >& lambda,
        int idx1, int idx2, int idx3, int idx4)
{
    auto element1 = interactionVolume.getSubVolumeElement(idx1);
    auto element2 = interactionVolume.getSubVolumeElement(idx2);
    auto element3 = interactionVolume.getSubVolumeElement(idx3);
    auto element4 = interactionVolume.getSubVolumeElement(idx4);

    if (element3 == element4 && element1.level() != element2.level())
    {
        return noTransmissibility;
    }

    // get global coordinate of cell centers
    const GlobalPosition& globalPos1 = element1.geometry().center();
    const GlobalPosition& globalPos2 = element2.geometry().center();
    const GlobalPosition& globalPos3 = element3.geometry().center();
    const GlobalPosition& globalPos4 = element4.geometry().center();

    const GlobalPosition& globalPosCenter = interactionVolume.getCenterPosition();

    const DimMatrix& K1 = problem_.spatialParams().intrinsicPermeability(element1);
    const DimMatrix& K2 = problem_.spatialParams().intrinsicPermeability(element2);
    const DimMatrix& K3 = problem_.spatialParams().intrinsicPermeability(element3);
    const DimMatrix& K4 = problem_.spatialParams().intrinsicPermeability(element4);

    const GlobalPosition& globalPosFace12 = interactionVolume.getFacePosition(idx1, 0);
    const GlobalPosition& globalPosFace23 = interactionVolume.getFacePosition(idx2, 0);

    // compute normal vectors nu1-nu7 in triangle R for first half edge
    DimVector nu1R1(0);
    R_.mv(globalPosFace12 - globalPos2, nu1R1);

    DimVector nu2R1(0);
    R_.mv(globalPos2 - globalPosFace23, nu2R1);

    DimVector nu3R1(0);
    R_.mv(globalPosFace23 - globalPos3, nu3R1);

    DimVector nu4R1(0);
    R_.mv(globalPos3 - globalPosCenter, nu4R1);

    DimVector nu5R1(0);
    R_.mv(globalPosCenter - globalPos1, nu5R1);

    DimVector nu6R1(0);
    R_.mv(globalPos1 - globalPosFace12, nu6R1);

    DimVector nu7R1(0);
    R_.mv(globalPosCenter - globalPos2, nu7R1);

    // compute T, i.e., the area of quadrilateral made by normal vectors 'nu'
    DimVector Rnu2R1(0);
    R_.mv(nu2R1, Rnu2R1);
    Scalar T1R1 = nu1R1 * Rnu2R1;

    DimVector Rnu4R1(0);
    R_.mv(nu4R1, Rnu4R1);
    Scalar T2R1 = nu3R1 * Rnu4R1;

    DimVector Rnu6R1(0);
    R_.mv(nu6R1, Rnu6R1);
    Scalar T3R1 = nu5R1 * Rnu6R1;

    // compute components needed for flux calculation, denoted as 'omega' and 'chi'
    DimVector K2nu1R1(0);
    K2.mv(nu1R1, K2nu1R1);
    DimVector K2nu2R1(0);
    K2.mv(nu2R1, K2nu2R1);
    DimVector K4nu3R1(0);
    K3.mv(nu3R1, K4nu3R1);
    DimVector K4nu4R1(0);
    K3.mv(nu4R1, K4nu4R1);
    DimVector K1nu5R1(0);
    K1.mv(nu5R1, K1nu5R1);
    DimVector K1nu6R1(0);
    K1.mv(nu6R1, K1nu6R1);

    DimVector Rnu1R1(0);
    R_.mv(nu1R1, Rnu1R1);

    DimVector &outerNormaln1R1 = interactionVolume.getNormal(idx2, 0);
    DimVector &outerNormaln2 = interactionVolume.getNormal(idx1, 0);

    Scalar omega111R1 = lambda[idx2][0] * (outerNormaln1R1 * K2nu1R1) * interactionVolume.getFaceArea(idx2, 0) / T1R1;
    Scalar omega112R1 = lambda[idx2][0] * (outerNormaln1R1 * K2nu2R1) * interactionVolume.getFaceArea(idx2, 0) / T1R1;
    Scalar omega211R1 = lambda[idx2][1] * (outerNormaln2 * K2nu1R1) * interactionVolume.getFaceArea(idx2, 1) / T1R1;
    Scalar omega212R1 = lambda[idx2][1] * (outerNormaln2 * K2nu2R1) * interactionVolume.getFaceArea(idx2, 1) / T1R1;
    Scalar omega123R1 = lambda[idx3][1] * (outerNormaln1R1 * K4nu3R1) * interactionVolume.getFaceArea(idx3, 1) / T2R1;
    Scalar omega124R1 = lambda[idx3][1] * (outerNormaln1R1 * K4nu4R1) * interactionVolume.getFaceArea(idx3, 1) / T2R1;
    Scalar omega235R1 = lambda[idx1][0] * (outerNormaln2 * K1nu5R1) * interactionVolume.getFaceArea(idx1, 0) / T3R1;
    Scalar omega236R1 = lambda[idx1][0] * (outerNormaln2 * K1nu6R1) * interactionVolume.getFaceArea(idx1, 0) / T3R1;
    Scalar chi711R1 = (nu7R1 * Rnu1R1) / T1R1;
    Scalar chi712R1 = (nu7R1 * Rnu2R1) / T1R1;

    // compute transmissibility matrix TR1 = CA^{-1}B+D
    DimMatrix C(0), A(0);
    Dune::FieldMatrix<Scalar, dim, 2 * dim - dim + 1> D(0), B(0);

    // evaluate matrix C, D, A, B
    C[0][0] = -omega111R1;
    C[0][1] = -omega112R1;
    C[1][0] = -omega211R1;
    C[1][1] = -omega212R1;

    D[0][0] = omega111R1 + omega112R1;
    D[1][0] = omega211R1 + omega212R1;

    A[0][0] = omega111R1 - omega124R1 - omega123R1 * chi711R1;
    A[0][1] = omega112R1 - omega123R1 * chi712R1;
    A[1][0] = omega211R1 - omega236R1 * chi711R1;
    A[1][1] = omega212R1 - omega235R1 - omega236R1 * chi712R1;

    B[0][0] = omega111R1 + omega112R1 + omega123R1 * (1.0 - chi711R1 - chi712R1);
    B[0][1] = -omega123R1 - omega124R1;
    B[1][0] = omega211R1 + omega212R1 + omega236R1 * (1.0 - chi711R1 - chi712R1);
    B[1][2] = -omega235R1 - omega236R1;

    // compute TR1
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));
    Dune::FieldMatrix<Scalar, dim, 2 * dim - dim + 1> TR1(D);

    // 2.use triangle L to compute the transmissibility of half edge
    const GlobalPosition& globalPosFace14 = interactionVolume.getFacePosition(idx1, 1);

    // compute normal vectors nu1-nu7 in triangle L for first half edge
    DimVector nu1L1(0);
    R_.mv(globalPosFace12 - globalPos1, nu1L1);

    DimVector nu2L1(0);
    R_.mv(globalPos1 - globalPosFace14, nu2L1);

    DimVector nu3L1(0);
    R_.mv(globalPosFace14 - globalPos4, nu3L1);

    DimVector nu4L1(0);
    R_.mv(globalPos4 - globalPosCenter, nu4L1);

    DimVector nu5L1(0);
    R_.mv(globalPosCenter - globalPos2, nu5L1);

    DimVector nu6L1(0);
    R_.mv(globalPos2 - globalPosFace12, nu6L1);

    DimVector nu7L1(0);
    R_.mv(globalPosCenter - globalPos1, nu7L1);

    // compute T, i.e., the area of quadrilateral made by normal vectors 'nu'
    DimVector Rnu2L1(0);
    R_.mv(nu2L1, Rnu2L1);
    Scalar T1L1 = nu1L1 * Rnu2L1;

    DimVector Rnu4L1(0);
    R_.mv(nu4L1, Rnu4L1);
    Scalar T2L1 = nu3L1 * Rnu4L1;

    DimVector Rnu6L1(0);
    R_.mv(nu6L1, Rnu6L1);
    Scalar T3L1 = nu5L1 * Rnu6L1;

    // compute components needed for flux calculation, denoted as 'omega' and 'chi'
    DimVector K1nu1L1(0);
    K1.mv(nu1L1, K1nu1L1);
    DimVector K1nu2L1(0);
    K1.mv(nu2L1, K1nu2L1);
    DimVector K3nu3L1(0);
    K4.mv(nu3L1, K3nu3L1);
    DimVector K3nu4L1(0);
    K4.mv(nu4L1, K3nu4L1);
    DimVector K2nu5L1(0);
    K2.mv(nu5L1, K2nu5L1);
    DimVector K2nu6L1(0);
    K2.mv(nu6L1, K2nu6L1);

    DimVector Rnu1L1(0);
    R_.mv(nu1L1, Rnu1L1);

    DimVector &outerNormaln1L1 = interactionVolume.getNormal(idx1, 1);

    Scalar omega111L1 = lambda[idx1][1] * (outerNormaln1L1 * K1nu1L1) * interactionVolume.getFaceArea(idx1, 1) / T1L1;
    Scalar omega112L1 = lambda[idx1][1] * (outerNormaln1L1 * K1nu2L1) * interactionVolume.getFaceArea(idx1, 1) / T1L1;
    Scalar omega211L1 = lambda[idx1][0] * (outerNormaln2 * K1nu1L1) * interactionVolume.getFaceArea(idx1, 0) / T1L1;
    Scalar omega212L1 = lambda[idx1][0] * (outerNormaln2 * K1nu2L1) * interactionVolume.getFaceArea(idx1, 0) / T1L1;
    Scalar omega123L1 = lambda[idx4][0] * (outerNormaln1L1 * K3nu3L1) * interactionVolume.getFaceArea(idx4, 0) / T2L1;
    Scalar omega124L1 = lambda[idx4][0] * (outerNormaln1L1 * K3nu4L1) * interactionVolume.getFaceArea(idx4, 0) / T2L1;
    Scalar omega235L1 = lambda[idx2][1] * (outerNormaln2 * K2nu5L1) * interactionVolume.getFaceArea(idx2, 1) / T3L1;
    Scalar omega236L1 = lambda[idx2][1] * (outerNormaln2 * K2nu6L1) * interactionVolume.getFaceArea(idx2, 1) / T3L1;
    Scalar chi711L1 = (nu7L1 * Rnu1L1) / T1L1;
    Scalar chi712L1 = (nu7L1 * Rnu2L1) / T1L1;

    // compute transmissibility matrix TL1 = CA^{-1}B+D
    C = 0;
    A = 0;
    D = 0;
    B = 0;

    // evaluate matrix C, D, A, B
    C[0][0] = -omega111L1;
    C[0][1] = -omega112L1;
    C[1][0] = -omega211L1;
    C[1][1] = -omega212L1;

    D[0][0] = omega111L1 + omega112L1;
    D[1][0] = omega211L1 + omega212L1;

    A[0][0] = omega111L1 - omega124L1 - omega123L1 * chi711L1;
    A[0][1] = omega112L1 - omega123L1 * chi712L1;
    A[1][0] = omega211L1 - omega236L1 * chi711L1;
    A[1][1] = omega212L1 - omega235L1 - omega236L1 * chi712L1;

    B[0][0] = omega111L1 + omega112L1 + omega123L1 * (1.0 - chi711L1 - chi712L1);
    B[0][1] = -omega123L1 - omega124L1;
    B[1][0] = omega211L1 + omega212L1 + omega236L1 * (1.0 - chi711L1 - chi712L1);
    B[1][2] = -omega235L1 - omega236L1;

    // compute TL1
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));
    Dune::FieldMatrix<Scalar, dim, 2 * dim - dim + 1> TL1(D);

    //selection criterion
    using std::abs;
    Scalar sR = abs(TR1[1][2] - TR1[1][0]);
    Scalar sL = abs(TL1[1][0] - TL1[1][2]);

    // 3.decide which triangle (which transmissibility coefficients) to use
    if (sR <= sL)
    {
        transmissibility = TR1;
        return rightTriangle;
    }
    else
    {
        transmissibility = TL1;
        return leftTriangle;
    }
}

/*!
 * \ingroup SequentialTwoPModel
 * \brief Calculates tranmissibility matrix
 *
 * Calculates tranmissibility matrix of an L-shape for a certain flux face.
 * Calculates only the transmissibility of the left L-shape (needed at hanging nodes HN).
 *
 * \param transmissibilityLeft Matrix for the resulting transmissibility
 * \param interactionVolume The interaction volume object (includes geometric information)
 * \param lambda Mobilities of cells 1-3
 * \param idx1 Index of cell 1 of the L-stencil
 * \param idx2 Index of cell 2 of the L-stencil
 * \param idx3 Index of cell 3 of the L-stencil
 */
template<class TypeTag>
int FvMpfaL2dTransmissibilityCalculator<TypeTag>::calculateLeftHNTransmissibility(
        TransmissibilityType& transmissibilityLeft,
        InteractionVolume& interactionVolume,
        std::vector<DimVector >& lambda,
        int idx1, int idx2, int idx3)
{
    auto element1 = interactionVolume.getSubVolumeElement(idx1);
    auto element2 = interactionVolume.getSubVolumeElement(idx2);
    auto element3 = interactionVolume.getSubVolumeElement(idx3);

    if (element1.level() != element3.level())
    {
        return noTransmissibility;
    }

    // get global coordinate of cell centers
    const GlobalPosition& globalPos1 = element1.geometry().center();
    const GlobalPosition& globalPos2 = element2.geometry().center();
    const GlobalPosition& globalPos3 = element3.geometry().center();

    const GlobalPosition& globalPosCenter = interactionVolume.getCenterPosition();

    const DimMatrix& K1 = problem_.spatialParams().intrinsicPermeability(element1);
    const DimMatrix& K2 = problem_.spatialParams().intrinsicPermeability(element2);
    const DimMatrix& K3 = problem_.spatialParams().intrinsicPermeability(element3);

    const GlobalPosition& globalPosFace12 = interactionVolume.getFacePosition(idx1, 0);
    DimVector &outerNormaln2 = interactionVolume.getNormal(idx1, 0);

    // compute transmissibility matrix TR1 = CA^{-1}B+D
    DimMatrix C(0), A(0);
    TransmissibilityType D(0), B(0);

    // 2.use triangle L to compute the transmissibility of half edge
    const GlobalPosition& globalPosFace14 = interactionVolume.getFacePosition(idx1, 1);

    // compute normal vectors nu1-nu7 in triangle L for first half edge
    DimVector nu1L1(0);
    R_.mv(globalPosFace12 - globalPos1, nu1L1);

    DimVector nu2L1(0);
    R_.mv(globalPos1 - globalPosFace14, nu2L1);

    DimVector nu3L1(0);
    R_.mv(globalPosFace14 - globalPos3, nu3L1);

    DimVector nu4L1(0);
    R_.mv(globalPos3 - globalPosCenter, nu4L1);

    DimVector nu5L1(0);
    R_.mv(globalPosCenter - globalPos2, nu5L1);

    DimVector nu6L1(0);
    R_.mv(globalPos2 - globalPosFace12, nu6L1);

    DimVector nu7L1(0);
    R_.mv(globalPosCenter - globalPos1, nu7L1);

    // compute T, i.e., the area of quadrilateral made by normal vectors 'nu'
    DimVector Rnu2L1(0);
    R_.mv(nu2L1, Rnu2L1);
    Scalar T1L1 = nu1L1 * Rnu2L1;

    DimVector Rnu4L1(0);
    R_.mv(nu4L1, Rnu4L1);
    Scalar T2L1 = nu3L1 * Rnu4L1;

    DimVector Rnu6L1(0);
    R_.mv(nu6L1, Rnu6L1);
    Scalar T3L1 = nu5L1 * Rnu6L1;

    // compute components needed for flux calculation, denoted as 'omega' and 'chi'
    DimVector K1nu1L1(0);
    K1.mv(nu1L1, K1nu1L1);
    DimVector K1nu2L1(0);
    K1.mv(nu2L1, K1nu2L1);
    DimVector K3nu3L1(0);
    K3.mv(nu3L1, K3nu3L1);
    DimVector K3nu4L1(0);
    K3.mv(nu4L1, K3nu4L1);
    DimVector K2nu5L1(0);
    K2.mv(nu5L1, K2nu5L1);
    DimVector K2nu6L1(0);
    K2.mv(nu6L1, K2nu6L1);

    DimVector Rnu1L1(0);
    R_.mv(nu1L1, Rnu1L1);

    DimVector &outerNormaln1L1 = interactionVolume.getNormal(idx1, 1);

    Scalar omega111L1 = lambda[idx1][1] * (outerNormaln1L1 * K1nu1L1) * interactionVolume.getFaceArea(idx1, 1) / T1L1;
    Scalar omega112L1 = lambda[idx1][1] * (outerNormaln1L1 * K1nu2L1) * interactionVolume.getFaceArea(idx1, 1) / T1L1;
    Scalar omega211L1 = lambda[idx1][0] * (outerNormaln2 * K1nu1L1) * interactionVolume.getFaceArea(idx1, 0) / T1L1;
    Scalar omega212L1 = lambda[idx1][0] * (outerNormaln2 * K1nu2L1) * interactionVolume.getFaceArea(idx1, 0) / T1L1;
    Scalar omega123L1 = lambda[idx3][0] * (outerNormaln1L1 * K3nu3L1) * interactionVolume.getFaceArea(idx3, 0) / T2L1;
    Scalar omega124L1 = lambda[idx3][0] * (outerNormaln1L1 * K3nu4L1) * interactionVolume.getFaceArea(idx3, 0) / T2L1;
    Scalar omega235L1 = lambda[idx2][1] * (outerNormaln2 * K2nu5L1) * interactionVolume.getFaceArea(idx2, 1) / T3L1;
    Scalar omega236L1 = lambda[idx2][1] * (outerNormaln2 * K2nu6L1) * interactionVolume.getFaceArea(idx2, 1) / T3L1;
    Scalar chi711L1 = (nu7L1 * Rnu1L1) / T1L1;
    Scalar chi712L1 = (nu7L1 * Rnu2L1) / T1L1;

    // compute transmissibility matrix TL1 = CA^{-1}B+D
    // evaluate matrix C, D, A, B
    C[0][0] = -omega111L1;
    C[0][1] = -omega112L1;
    C[1][0] = -omega211L1;
    C[1][1] = -omega212L1;

    D[0][0] = omega111L1 + omega112L1;
    D[1][0] = omega211L1 + omega212L1;

    A[0][0] = omega111L1 - omega124L1 - omega123L1 * chi711L1;
    A[0][1] = omega112L1 - omega123L1 * chi712L1;
    A[1][0] = omega211L1 - omega236L1 * chi711L1;
    A[1][1] = omega212L1 - omega235L1 - omega236L1 * chi712L1;

    B[0][0] = omega111L1 + omega112L1 + omega123L1 * (1.0 - chi711L1 - chi712L1);
    B[0][1] = -omega123L1 - omega124L1;
    B[1][0] = omega211L1 + omega212L1 + omega236L1 * (1.0 - chi711L1 - chi712L1);
    B[1][2] = -omega235L1 - omega236L1;

    // compute TL1
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));

    // 3.decide which triangle (which transmissibility coefficients) to use
    transmissibilityLeft = D;
    return leftTriangle;
}

/*!
 * \ingroup SequentialTwoPModel
 * \brief Calculates tranmissibility matrix
 *
 * Calculates tranmissibility matrix of an L-shape for a certain flux face.
 * Calculates only the transmissibility of the right L-shape (needed at hanging nodes HN).
 *
 * \param transmissibilityRight Matrix for the resulting transmissibility
 * \param interactionVolume The interaction volume object (includes geometric information)
 * \param lambda Mobilities of cells 1-3
 * \param idx1 Index of cell 1 of the L-stencil
 * \param idx2 Index of cell 2 of the L-stencil
 * \param idx3 Index of cell 3 of the L-stencil
 */
template<class TypeTag>
int FvMpfaL2dTransmissibilityCalculator<TypeTag>::calculateRightHNTransmissibility(
        TransmissibilityType& transmissibilityRight,
        InteractionVolume& interactionVolume,
        std::vector<DimVector >& lambda,
        int idx1, int idx2, int idx3)
{
    auto element1 = interactionVolume.getSubVolumeElement(idx1);
    auto element2 = interactionVolume.getSubVolumeElement(idx2);
    auto element3 = interactionVolume.getSubVolumeElement(idx3);

    if (element2.level() != element3.level())
    {
        return noTransmissibility;
    }

    // get global coordinate of cell centers
    const GlobalPosition& globalPos1 = element1.geometry().center();
    const GlobalPosition& globalPos2 = element2.geometry().center();
    const GlobalPosition& globalPos3 = element3.geometry().center();

    const GlobalPosition& globalPosCenter = interactionVolume.getCenterPosition();

    const DimMatrix& K1 = problem_.spatialParams().intrinsicPermeability(element1);
    const DimMatrix& K2 = problem_.spatialParams().intrinsicPermeability(element2);
    const DimMatrix& K3 = problem_.spatialParams().intrinsicPermeability(element3);

    const GlobalPosition& globalPosFace12 = interactionVolume.getFacePosition(idx1, 0);
    const GlobalPosition& globalPosFace23 = interactionVolume.getFacePosition(idx2, 0);

    // compute normal vectors nu1-nu7 in triangle R for first half edge
    DimVector nu1R1(0);
    R_.mv(globalPosFace12 - globalPos2, nu1R1);

    DimVector nu2R1(0);
    R_.mv(globalPos2 - globalPosFace23, nu2R1);

    DimVector nu3R1(0);
    R_.mv(globalPosFace23 - globalPos3, nu3R1);

    DimVector nu4R1(0);
    R_.mv(globalPos3 - globalPosCenter, nu4R1);

    DimVector nu5R1(0);
    R_.mv(globalPosCenter - globalPos1, nu5R1);

    DimVector nu6R1(0);
    R_.mv(globalPos1 - globalPosFace12, nu6R1);

    DimVector nu7R1(0);
    R_.mv(globalPosCenter - globalPos2, nu7R1);

    // compute T, i.e., the area of quadrilateral made by normal vectors 'nu'
    DimVector Rnu2R1(0);
    R_.mv(nu2R1, Rnu2R1);
    Scalar T1R1 = nu1R1 * Rnu2R1;

    DimVector Rnu4R1(0);
    R_.mv(nu4R1, Rnu4R1);
    Scalar T2R1 = nu3R1 * Rnu4R1;

    DimVector Rnu6R1(0);
    R_.mv(nu6R1, Rnu6R1);
    Scalar T3R1 = nu5R1 * Rnu6R1;

    // compute components needed for flux calculation, denoted as 'omega' and 'chi'
    DimVector K2nu1R1(0);
    K2.mv(nu1R1, K2nu1R1);
    DimVector K2nu2R1(0);
    K2.mv(nu2R1, K2nu2R1);
    DimVector K3nu3R1(0);
    K3.mv(nu3R1, K3nu3R1);
    DimVector K3nu4R1(0);
    K3.mv(nu4R1, K3nu4R1);
    DimVector K1nu5R1(0);
    K1.mv(nu5R1, K1nu5R1);
    DimVector K1nu6R1(0);
    K1.mv(nu6R1, K1nu6R1);

    DimVector Rnu1R1(0);
    R_.mv(nu1R1, Rnu1R1);

    DimVector &outerNormaln1R1 = interactionVolume.getNormal(idx2, 0);
    DimVector &outerNormaln2 = interactionVolume.getNormal(idx1, 0);

    Scalar omega111R1 = lambda[idx2][0] * (outerNormaln1R1 * K2nu1R1) * interactionVolume.getFaceArea(idx2, 0) / T1R1;
    Scalar omega112R1 = lambda[idx2][0] * (outerNormaln1R1 * K2nu2R1) * interactionVolume.getFaceArea(idx2, 0) / T1R1;
    Scalar omega211R1 = lambda[idx2][1] * (outerNormaln2 * K2nu1R1) * interactionVolume.getFaceArea(idx2, 1) / T1R1;
    Scalar omega212R1 = lambda[idx2][1] * (outerNormaln2 * K2nu2R1) * interactionVolume.getFaceArea(idx2, 1) / T1R1;
    Scalar omega123R1 = lambda[idx3][1] * (outerNormaln1R1 * K3nu3R1) * interactionVolume.getFaceArea(idx3, 1) / T2R1;
    Scalar omega124R1 = lambda[idx3][1] * (outerNormaln1R1 * K3nu4R1) * interactionVolume.getFaceArea(idx3, 1) / T2R1;
    Scalar omega235R1 = lambda[idx1][0] * (outerNormaln2 * K1nu5R1) * interactionVolume.getFaceArea(idx1, 0) / T3R1;
    Scalar omega236R1 = lambda[idx1][0] * (outerNormaln2 * K1nu6R1) * interactionVolume.getFaceArea(idx1, 0) / T3R1;
    Scalar chi711R1 = (nu7R1 * Rnu1R1) / T1R1;
    Scalar chi712R1 = (nu7R1 * Rnu2R1) / T1R1;

    // compute transmissibility matrix TR1 = CA^{-1}B+D
    DimMatrix C(0), A(0);
    TransmissibilityType D(0), B(0);

    // evaluate matrix C, D, A, B
    C[0][0] = -omega111R1;
    C[0][1] = -omega112R1;
    C[1][0] = -omega211R1;
    C[1][1] = -omega212R1;

    D[0][0] = omega111R1 + omega112R1;
    D[1][0] = omega211R1 + omega212R1;

    A[0][0] = omega111R1 - omega124R1 - omega123R1 * chi711R1;
    A[0][1] = omega112R1 - omega123R1 * chi712R1;
    A[1][0] = omega211R1 - omega236R1 * chi711R1;
    A[1][1] = omega212R1 - omega235R1 - omega236R1 * chi712R1;

    B[0][0] = omega111R1 + omega112R1 + omega123R1 * (1.0 - chi711R1 - chi712R1);
    B[0][1] = -omega123R1 - omega124R1;
    B[1][0] = omega211R1 + omega212R1 + omega236R1 * (1.0 - chi711R1 - chi712R1);
    B[1][2] = -omega235R1 - omega236R1;

    // compute TR1
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));

    transmissibilityRight = D;
    return rightTriangle;
}
} // end namespace Dumux
#endif
