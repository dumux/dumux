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
 * \brief Provides methods for transmissibility calculation in 3-d.
 */
#ifndef DUMUX_FVMPFAL3D_TRANSMISSIBILITYCALCULATOR_HH
#define DUMUX_FVMPFAL3D_TRANSMISSIBILITYCALCULATOR_HH

// dumux environment
#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Provides methods for transmissibility calculation in 3-d.
 *
 *  The transmissibilities are calculated using the MPFA L-method.
 *
 * Aavatsmark et al. 2010 \cite aavatsmark2010 <BR>
 * Wolff et al. 2013 \cite wolff2013b <BR>
 *
 *  Various parameters can be defined via an input parameter file or the property system:
 *
 * MPFAEnableSimpleLStencil - enables the two centered flux stencils
 * MPFAEnableComplexLStencil - enables the two non-centered flux stencils
 * MPFAEnableTPFA - enables the use of TPFA if neighboring cells are of the same grid level
 * MPFATransmissibilityCriterion - 0: Criterion of [1] (default), 1: Criterion of [2]
 *
 */
template<class TypeTag>
class FvMpfaL3dTransmissibilityCalculator
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    enum TransCriterion
        {
            sDiff = 0, sSum = 1
        };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using DimVector = Dune::FieldVector<Scalar, dim>;

    using InteractionVolume = GetPropType<TypeTag, Properties::MPFAInteractionVolume>;


public:
    using TransmissibilityType = Dune::FieldMatrix<Scalar, dim, 2*dim - dim + 1>;//!< Type of the transmissibility matrix

    int chooseTransmissibility(TransmissibilityType& transmissibilityOne,
                               TransmissibilityType& transmissibilityTwo, int lTypeOne, int lTypeTwo);

    int transmissibility(Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                         InteractionVolume& interactionVolume,
                         std::vector<DimVector >& lambda,
                         int idx1, int idx2, int idx3, int idx4, int idx5, int idx6);

    int transmissibility(Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                         InteractionVolume& interactionVolume,
                         std::vector<DimVector >& lambda,
                         int idx1, int idx2, int idx3, int idx4, int idx5, int idx6,
                         Dune::FieldVector<bool, 4> &useCases);

    int transmissibilityTPFA(Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                             InteractionVolume& interactionVolume,
                             std::vector<DimVector >& lambda,
                             int idx1, int idx2);

    int transmissibilityCaseOne(
                                Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                InteractionVolume& interactionVolume,
                                std::vector<DimVector >& lambda,
                                int idx1, int idx2, int idx3, int idx5);

    int transmissibilityCaseTwo(
                                Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                InteractionVolume& interactionVolume,
                                std::vector<DimVector >& lambda,
                                int idx1, int idx2, int idx4, int idx6);
    int transmissibilityCaseThree(
                                  Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                  InteractionVolume& interactionVolume,
                                  std::vector<DimVector >& lambda,
                                  int idx1, int idx2, int idx4, int idx5);

    int transmissibilityCaseFour(
                                 Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                 InteractionVolume& interactionVolume,
                                 std::vector<DimVector >& lambda,
                                 int idx1, int idx2, int idx3, int idx6);

    /*!
     * \brief Constructs a FvMpfaL3dTransmissibilityCalculator object
     *
     * \param problem A problem class object
     */
    FvMpfaL3dTransmissibilityCalculator(Problem& problem) :
        problem_(problem)
    {
        if (dim != 3)
        {
            DUNE_THROW(Dune::NotImplemented, "Dimension not supported!");
        }

        enableSimpleLStencil_ = getParam<bool>("MPFA.EnableSimpleLStencil", true);
        enableComplexLStencil_ = getParam<bool>("MPFA.EnableComplexLStencil", true);
        enableTPFA_= getParam<bool>("MPFA.EnableTPFA", false);

        if (!enableSimpleLStencil_ && !enableComplexLStencil_)
        {
            std::cout<<"MPFA disabled: Use TPFA!\n";
            enableTPFA_ = true;
        }

        transChoiceThreshold_ = getParam<Scalar>("MPFA.TransmissibilityCriterionThreshold", 1e-8);

        transCriterion_ = getParam<int>("MPFA.TransmissibilityCriterion", 0);
        if (transCriterion_ == sDiff)
        {
            std::cout << "MPFAL 3D: Using standard transmissibility criterion\n";
        }
        else if (transCriterion_ == sSum)
        {
            std::cout << "MPFAL 3D: Using accumulative transmissibility criterion\n";
        }
        else
        {
            transCriterion_ = sDiff;
            std::cout << "MPFAL 3D: Defined transmissiblity criterion not implemented! Using standard criterion!\n";
        }
    }

private:
    Problem& problem_;
    bool enableSimpleLStencil_;
    bool enableComplexLStencil_;
    bool enableTPFA_;
    Scalar transChoiceThreshold_;
    int transCriterion_;
};

/*!
 * \brief Compares two transmissibility matrices according to a L-selection criterion
 *
 * Compares two transmissibility matrices according to the L-selection criterion which is chosen via the parameter/property
 * MPFATransmissibilityCriterion (Criterion of [1], 1: Criterion of [2]) and returns the number of the preferable L-shape (1-4).
 *
 * \param transmissibilityOne A first transmissibility matrix
 * \param transmissibilityTwo A second transmissibility matrix
 * \param lTypeOne Type of the L-shape of the first matrix (1-4)
 * \param lTypeTwo Type of the L-shape of the second matrix (1-4)
 *
 * \return The type of the preferable L-shape (1-4)
 */
template<class TypeTag>
int FvMpfaL3dTransmissibilityCalculator<TypeTag>::chooseTransmissibility(TransmissibilityType& transmissibilityOne,
                                                                         TransmissibilityType& transmissibilityTwo,
                                                                         int lTypeOne, int lTypeTwo)
{
    if (transCriterion_ == sDiff)
    {
        using std::abs;
        Scalar sOne = abs(transmissibilityOne[0][0] - transmissibilityOne[0][1]);
        Scalar sTwo = abs(transmissibilityTwo[0][0] - transmissibilityTwo[0][1]);

        //Decide whether to take case1 or case2
        if (sOne < sTwo - transChoiceThreshold_)
        {
            return lTypeOne;
        }
        else
        {
            return lTypeTwo;
        }
    }
    else if (transCriterion_ == sSum)
    {

        Scalar tSumOne = 0;
        Scalar tSumTwo = 0;
        using std::abs;

        if (lTypeOne == 1)
            tSumOne = abs(transmissibilityOne[0][0] + transmissibilityOne[0][2] + transmissibilityOne[0][3]);
        else if (lTypeOne == 2)
            tSumOne = abs(transmissibilityOne[0][1] + transmissibilityOne[0][2] + transmissibilityOne[0][3]);
        else if (lTypeOne == 3)
            tSumOne = abs(transmissibilityOne[0][0] + transmissibilityOne[0][3]);
        else if (lTypeOne == 4)
            tSumOne = abs(transmissibilityOne[0][0] + transmissibilityOne[0][2]);
        else
            DUNE_THROW(Dune::NotImplemented,"Transmissibility type not implemented");

        if (lTypeTwo == 1)
            tSumTwo = abs(transmissibilityTwo[0][0] + transmissibilityTwo[0][2] + transmissibilityTwo[0][3]);
        else if (lTypeTwo == 2)
            tSumTwo = abs(transmissibilityTwo[0][1] + transmissibilityTwo[0][2] + transmissibilityTwo[0][3]);
        else if (lTypeTwo == 3)
            tSumTwo = abs(transmissibilityTwo[0][0] + transmissibilityTwo[0][3]);
        else if (lTypeTwo == 4)
            tSumTwo = abs(transmissibilityTwo[0][0] + transmissibilityTwo[0][2]);
        else
            DUNE_THROW(Dune::NotImplemented,"Transmissibility type not implemented");


        if (tSumOne > tSumTwo + transChoiceThreshold_)
        {
            return lTypeOne;
        }
        else
        {
            return lTypeTwo;
        }
    }
    else
    {
        DUNE_THROW(Dune::NotImplemented,"Transmissibility criterion not implemented");
    }

}

/*!
 * \brief Calculates the transmissibility matrix for the flux face between the cells with the local index idx1 and idx2
 *
 *  \param transmissibility Reference to the local transmissibility matrix
 *  \param interactionVolume An interaction volume object
 *  \param lambda A vector containing the mobilities of the interaction volume cells in the order of the local interaction volume index
 *  \param idx1 Local index of cell 1 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx2 Local index of cell 2 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx3 Local index of cell 3 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx4 Local index of cell 4 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx5 Local index of cell 5 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx6 Local index of cell 6 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *
 *  \return The type of the chosen L-shape (1-4)
 */
template<class TypeTag>
int FvMpfaL3dTransmissibilityCalculator<TypeTag>::transmissibility(TransmissibilityType& transmissibility,
                                                                   InteractionVolume& interactionVolume,
                                                                   std::vector<DimVector>& lambda,
                                                                   int idx1, int idx2, int idx3, int idx4,
                                                                   int idx5, int idx6)
{
    int level1 = interactionVolume.getSubVolumeElement(idx1).level();
    int level2 = interactionVolume.getSubVolumeElement(idx2).level();

    if (enableTPFA_ && level1 == level2)
    {
        return transmissibilityTPFA(transmissibility, interactionVolume, lambda, idx1, idx2);
    }
    else if (enableTPFA_ && !enableSimpleLStencil_ && !enableComplexLStencil_)
    {
        return transmissibilityTPFA(transmissibility, interactionVolume, lambda, idx1, idx2);
    }

    if (enableSimpleLStencil_ && enableComplexLStencil_)
    {
        TransmissibilityType transTemp(0);

        // decide which L-stencil (which transmissibility coefficients) to use
        // calculate transmissibility differences
        int transTypeTemp = transmissibilityCaseOne(transTemp, interactionVolume, lambda, idx1, idx2, idx3, idx5);
        int transType = transmissibilityCaseTwo(transmissibility, interactionVolume, lambda, idx1, idx2, idx4, idx6);

        transType = chooseTransmissibility(transTemp, transmissibility, 1, 2);

        if (transType == 1)
            transmissibility = transTemp;

        TransmissibilityType transCaseThree(0);

        int DUNE_UNUSED transTypeCaseThree = transmissibilityCaseThree(transCaseThree, interactionVolume,
                                                                               lambda, idx1, idx2, idx4, idx5);
        transTypeTemp = transmissibilityCaseFour(transTemp, interactionVolume, lambda, idx1, idx2, idx3, idx6);

        transTypeTemp = chooseTransmissibility(transCaseThree, transTemp, 3, 4);

        if (transTypeTemp == 3)
            transTemp = transCaseThree;

        transType = chooseTransmissibility(transmissibility, transTemp, transType, transTypeTemp);

        if (transType == transTypeTemp)
            transmissibility = transTemp;

        return transType;
    }
    else if (enableSimpleLStencil_)
    {
        TransmissibilityType transTemp(0);

        int DUNE_UNUSED transTypeTemp = transmissibilityCaseOne(transTemp, interactionVolume,
                                                                                     lambda, idx1, idx2, idx3, idx5);
        int transType = transmissibilityCaseTwo(transmissibility, interactionVolume, lambda, idx1, idx2, idx4, idx6);

        transType = chooseTransmissibility(transTemp, transmissibility, 1, 2);

        if (transType == 1)
            transmissibility = transTemp;

        return transType;
    }
    else if (enableComplexLStencil_)
    {
        TransmissibilityType transTemp(0);

        int DUNE_UNUSED transTypeTemp = transmissibilityCaseThree(transTemp, interactionVolume, lambda, idx1, idx2, idx4, idx5);
        int transType = transmissibilityCaseFour(transmissibility, interactionVolume, lambda, idx1, idx2, idx3, idx6);

        transType = chooseTransmissibility(transTemp, transmissibility, 3, 4);

        if (transType == 3)
            transmissibility = transTemp;

        return transType;
    }
    else
    {
        DUNE_THROW(Dune::NotImplemented,"Stencil type not implemented");
    }
}

/*!
 * \brief Calculates the transmissibility matrix for the flux face between the cells with the local index idx1 and idx2
 *
 *  \param transmissibility Reference to the local transmissibility matrix
 *  \param interactionVolume An interaction volume object
 *  \param lambda A vector containing the mobilities of the interaction volume cells in the order of the local interaction volume index
 *  \param idx1 Local index of cell 1 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx2 Local index of cell 2 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx3 Local index of cell 3 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx4 Local index of cell 4 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx5 Local index of cell 5 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx6 Local index of cell 6 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param useCases Vector of enabling/disabling single L-shape types (Can be necessary for certain configurations in the case of hanging nodes)
 *
 *  \return The type of the chosen L-shape (1-4)
 */
template<class TypeTag>
int FvMpfaL3dTransmissibilityCalculator<TypeTag>::transmissibility(TransmissibilityType& transmissibility,
                                                                   InteractionVolume& interactionVolume,
                                                                   std::vector<DimVector>& lambda, int idx1,
                                                                   int idx2, int idx3, int idx4, int idx5,
                                                                   int idx6, Dune::FieldVector<bool, 4> &useCases)
{
    int level1 = interactionVolume.getSubVolumeElement(idx1).level();
    int level2 = interactionVolume.getSubVolumeElement(idx2).level();

    if (enableTPFA_ && level1 == level2)
    {
        return transmissibilityTPFA(transmissibility, interactionVolume, lambda, idx1, idx2);
    }
    else if (enableTPFA_ && !enableSimpleLStencil_ && !enableComplexLStencil_)
    {
        return transmissibilityTPFA(transmissibility, interactionVolume, lambda, idx1, idx2);
    }
    else if (enableTPFA_ && !enableSimpleLStencil_ && !useCases[2] && !useCases[3])
    {
        return transmissibilityTPFA(transmissibility, interactionVolume, lambda, idx1, idx2);
    }
    else if (enableTPFA_ && !enableComplexLStencil_ && !useCases[0] && !useCases[1])
    {
        return transmissibilityTPFA(transmissibility, interactionVolume, lambda, idx1, idx2);
    }

    if (enableSimpleLStencil_ && enableComplexLStencil_)
    {
        TransmissibilityType transTemp(0);
        int transTypeTemp = 0;
        int transType = 0;

        // decide which L-stencil (which transmissibility coefficients) to use
        // calculate transmissibility differences

        if (useCases[0])
            transTypeTemp = transmissibilityCaseOne(transTemp, interactionVolume, lambda, idx1, idx2, idx3, idx5);

        if (useCases[1])
            transType = transmissibilityCaseTwo(transmissibility, interactionVolume, lambda, idx1, idx2, idx4, idx6);

        if (useCases[0] && useCases[1])
        {
            transType = chooseTransmissibility(transTemp, transmissibility, 1, 2);

            if (transType == 1)
                transmissibility = transTemp;
        }
        else if (useCases[0])
        {
            transmissibility = transTemp;
            transType = transTypeTemp;
        }

        TransmissibilityType transCaseThree(0);
        int transTypeCaseThree = 0;

        if (useCases[2])
            transTypeCaseThree = transmissibilityCaseThree(transCaseThree, interactionVolume, lambda, idx1, idx2, idx4, idx5);

        if (useCases[3])
            transTypeTemp = transmissibilityCaseFour(transTemp, interactionVolume, lambda, idx1, idx2, idx3, idx6);

        if (useCases[2] && useCases[3])
        {
            transTypeTemp = chooseTransmissibility(transCaseThree, transTemp, 3, 4);

            if (transTypeTemp == 3)
                transTemp = transCaseThree;
        }
        else if (useCases[2])
        {
            transTemp = transCaseThree;
            transTypeTemp = transTypeCaseThree;
        }

        if ((useCases[0] || useCases[1]) && (useCases[2] || useCases[3]))
        {
            transType = chooseTransmissibility(transmissibility, transTemp, transType, transTypeTemp);

            if (transType == transTypeTemp)
                transmissibility = transTemp;
        }
        else if (useCases[2] || useCases[3])
        {
            transmissibility = transTemp;
            transType = transTypeTemp;
        }

        return transType;
    }
    else if (enableSimpleLStencil_)
    {
        if (useCases[0] && useCases[1])
        {
            TransmissibilityType transTemp(0);

            int DUNE_UNUSED transTypeTemp = transmissibilityCaseOne(transTemp, interactionVolume,
                                                                    lambda, idx1, idx2, idx3, idx5);
            int transType = transmissibilityCaseTwo(transmissibility, interactionVolume,
                                                                    lambda, idx1, idx2, idx4, idx6);

            transType = chooseTransmissibility(transTemp, transmissibility, 1, 2);

            if (transType == 1)
                transmissibility = transTemp;

            return transType;
        }
        else if (useCases[0])
        {
            return transmissibilityCaseOne(transmissibility, interactionVolume, lambda, idx1, idx2, idx3, idx5);
        }
        else if (useCases[1])
        {
            return transmissibilityCaseTwo(transmissibility, interactionVolume, lambda, idx1, idx2, idx4, idx6);
        }
        else
        {
            DUNE_THROW(Dune::RangeError,"Only simple L-stencils allowed but no simple stencil chosen!");
        }
    }
    else if (enableComplexLStencil_)
    {
        if (useCases[2] && useCases[3])
        {
            TransmissibilityType transTemp(0);

            int DUNE_UNUSED transTypeTemp = transmissibilityCaseThree(transTemp, interactionVolume,
                                                                         lambda, idx1, idx2, idx4, idx5);
            int transType = transmissibilityCaseFour(transmissibility, interactionVolume,
                                                                         lambda, idx1, idx2, idx3, idx6);

            transType = chooseTransmissibility(transTemp, transmissibility, 3, 4);

            if (transType == 3)
                transmissibility = transTemp;

            return transType;
        }
        else if (useCases[2])
        {
            return transmissibilityCaseThree(transmissibility, interactionVolume, lambda, idx1, idx2, idx4, idx5);
        }
        else if (useCases[3])
        {
            return transmissibilityCaseFour(transmissibility, interactionVolume, lambda, idx1, idx2, idx3, idx6);
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented,"Transmissibility type not implemented");
        }
    }
    else
    {
        DUNE_THROW(Dune::NotImplemented,"Stencil type not implemented");
    }
}

/*!
 * \brief Calculates a TPFA transmissibility matrix for the flux face between the cells with the local index idx1 and idx2
 *
 *  \param transmissibility Reference to the local transmissibility matrix
 *  \param interactionVolume An interaction volume object
 *  \param lambda A vector containing the mobilities of the interaction volume cells in the order of the local interaction volume index
 *  \param idx1 Local index of cell 1 for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx2 Local index of cell 2 for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *
 */
template<class TypeTag>
int FvMpfaL3dTransmissibilityCalculator<TypeTag>::transmissibilityTPFA(
                                                         Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                                         InteractionVolume& interactionVolume,
                                                         std::vector<DimVector >& lambda, int idx1, int idx2)
{
    auto element1 = interactionVolume.getSubVolumeElement(idx1);
    auto element2 = interactionVolume.getSubVolumeElement(idx2);

    // get global coordinate of cell centers
    const GlobalPosition& globalPos1 = element1.geometry().center();
    const GlobalPosition& globalPos2 = element2.geometry().center();

    GlobalPosition distVec = globalPos1 - globalPos2;

    Scalar dist = distVec.two_norm();

    DimVector outerNormal(0);
    Scalar lambda12 =0;
    Scalar lambda21 =0;
    Scalar faceArea = 0;

    if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 2) || (idx1 == 2 && idx2 == 0))
    {

        outerNormal = interactionVolume.getNormal(idx1, 0);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][1];

        faceArea = interactionVolume.getFaceArea(idx1,0);
    }
    else if ((idx1 == 5 && idx2 == 4) || (idx1 == 4 && idx2 == 6) || (idx1 == 6 && idx2 == 7) || (idx1 == 7 && idx2 == 5))
    {
        outerNormal = interactionVolume.getNormal(idx1, 2);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][1];

        faceArea = interactionVolume.getFaceArea(idx1,2);
    }
    else if ((idx1 == 4 && idx2 == 0) || (idx1 == 7 && idx2 == 3))
    {
        outerNormal = interactionVolume.getNormal(idx1, 0);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][2];

        faceArea = interactionVolume.getFaceArea(idx1,0);
    }
    else
    {
        outerNormal = interactionVolume.getNormal(idx1, 2);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][0];

        faceArea = interactionVolume.getFaceArea(idx1,2);
    }

    const DimMatrix& K1 = problem_.spatialParams().intrinsicPermeability(element1);
    const DimMatrix& K2 = problem_.spatialParams().intrinsicPermeability(element2);

    // compute vectorized permeabilities
    DimVector permeability(0);
    K1.mv(outerNormal, permeability);

    Scalar perm1 = permeability * outerNormal;

    permeability = 0;
    K2.mv(outerNormal, permeability);

    Scalar perm2 = permeability * outerNormal;

    Scalar meanMob = 0;
    if ((lambda12 + lambda21) > 0.0)
        meanMob = 2*lambda12 * perm1 * lambda21 * perm2/(lambda12 * perm1 + lambda21 * perm2);

    Scalar entry = meanMob / dist * faceArea;

    transmissibility = 0;
    transmissibility[0][0] = entry;
    transmissibility[0][1] = -entry;

    return 1;
}

/*!
 * \brief Calculates the transmissibility matrix for the flux face between the cells with the local index idx1 and idx2 using the L-shape type 1
 *
 * For more details on L-shape type 1 see:
 *
 * Wolff, 2013 \cite wolff2013a <BR>
 * Wolff et al. 2013 \cite wolff2013b <BR>
 *
 *  \param transmissibility Reference to the local transmissibility matrix
 *  \param interactionVolume An interaction volume object
 *  \param lambda A vector containing the mobilities of the interaction volume cells in the order of the local interaction volume index
 *  \param idx1 Local index of cell 1 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx2 Local index of cell 2 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx3 Local index of cell 3 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa
 *  \param idx5 Local index of cell 5 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *
 *  \return 1 (L-shape 1)
 */
template<class TypeTag>
int FvMpfaL3dTransmissibilityCalculator<TypeTag>::transmissibilityCaseOne(
                                                         Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                                         InteractionVolume& interactionVolume,
                                                         std::vector<DimVector >& lambda, int idx1, int idx2, int idx3, int idx5)
{
    auto element1 = interactionVolume.getSubVolumeElement(idx1);
    auto element2 = interactionVolume.getSubVolumeElement(idx2);
    auto element3 = interactionVolume.getSubVolumeElement(idx3);
    auto element5 = interactionVolume.getSubVolumeElement(idx5);

    // get global coordinate of cell centers
    const GlobalPosition& globalPos1 = element1.geometry().center();
    const GlobalPosition& globalPos2 = element2.geometry().center();
    const GlobalPosition& globalPos3 = element3.geometry().center();
    const GlobalPosition& globalPos5 = element5.geometry().center();

    const GlobalPosition& globalPosCenter = interactionVolume.getCenterPosition();

    const DimMatrix& K1 = problem_.spatialParams().intrinsicPermeability(element1);
    const DimMatrix& K2 = problem_.spatialParams().intrinsicPermeability(element2);
    const DimMatrix& K3 = problem_.spatialParams().intrinsicPermeability(element3);
    const DimMatrix& K5 = problem_.spatialParams().intrinsicPermeability(element5);

    // 1.use centered L-stencil 1 to compute the transmissibility of 1/4 face
    GlobalPosition edgeCoord1213(0);
    GlobalPosition edgeCoord1215(0);
    GlobalPosition edgeCoord1315(0);

    GlobalPosition globalPosFace12(0);
    GlobalPosition globalPosFace13(0);
    GlobalPosition globalPosFace15(0);

    DimVector outerNormaln1(0);
    DimVector outerNormaln2(0);
    DimVector outerNormaln3(0);

    Scalar lambda12 = 0;
    Scalar lambda21 = 0;
    Scalar lambda13 = 0;
    Scalar lambda31 = 0;
    Scalar lambda15 = 0;
    Scalar lambda51 = 0;
    Scalar faceArea12 = 0;
    Scalar faceArea21 = 0;
    Scalar faceArea13 = 0;
    Scalar faceArea31 = 0;
    Scalar faceArea15 = 0;
    Scalar faceArea51 = 0;

    if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 2) || (idx1 == 2 && idx2 == 0))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace13 = interactionVolume.getFacePosition(idx1,1);
        globalPosFace15 = interactionVolume.getFacePosition(idx1,2);

        edgeCoord1213 = interactionVolume.getEdgePosition(idx1, 0, 1);
        edgeCoord1215 = interactionVolume.getEdgePosition(idx1, 0, 0);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 1, 1);

        outerNormaln1 = interactionVolume.getNormal(idx1, 0);
        outerNormaln2 = interactionVolume.getNormal(idx1, 1);
        outerNormaln3 = interactionVolume.getNormal(idx1, 2);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][1];
        lambda13 = lambda[idx1][1];
        lambda31 = lambda[idx3][0];
        lambda15 = lambda[idx1][2];
        lambda51 = lambda[idx5][0];

        faceArea12 = interactionVolume.getFaceArea(idx1,0);
        faceArea21 = interactionVolume.getFaceArea(idx2,1);
        faceArea13 = interactionVolume.getFaceArea(idx1,1);
        faceArea31 = interactionVolume.getFaceArea(idx3,0);
        faceArea15 = interactionVolume.getFaceArea(idx1,2);
        faceArea51 = interactionVolume.getFaceArea(idx5,0);
    }
    else if ((idx1 == 5 && idx2 == 4) || (idx1 == 4 && idx2 == 6) || (idx1 == 6 && idx2 == 7) || (idx1 == 7 && idx2 == 5))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace13 = interactionVolume.getFacePosition(idx1,1);
        globalPosFace15 = interactionVolume.getFacePosition(idx1,0);

        edgeCoord1213 = interactionVolume.getEdgePosition(idx1, 2, 1);
        edgeCoord1215 = interactionVolume.getEdgePosition(idx1, 2, 0);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 1, 1);

        outerNormaln1 = interactionVolume.getNormal(idx1, 2);
        outerNormaln2 = interactionVolume.getNormal(idx1, 1);
        outerNormaln3 = interactionVolume.getNormal(idx1, 0);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][1];
        lambda13 = lambda[idx1][1];
        lambda31 = lambda[idx3][2];
        lambda15 = lambda[idx1][0];
        lambda51 = lambda[idx5][2];

        faceArea12 = interactionVolume.getFaceArea(idx1,2);
        faceArea21 = interactionVolume.getFaceArea(idx2,1);
        faceArea13 = interactionVolume.getFaceArea(idx1,1);
        faceArea31 = interactionVolume.getFaceArea(idx3,2);
        faceArea15 = interactionVolume.getFaceArea(idx1,0);
        faceArea51 = interactionVolume.getFaceArea(idx5,2);
    }
    else if ((idx1 == 4 && idx2 == 0) || (idx1 == 7 && idx2 == 3))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace13 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace15 = interactionVolume.getFacePosition(idx1,1);

        edgeCoord1213 = interactionVolume.getEdgePosition(idx1, 0, 1);
        edgeCoord1215 = interactionVolume.getEdgePosition(idx1, 0, 0);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 2, 1);

        outerNormaln1 = interactionVolume.getNormal(idx1, 0);
        outerNormaln2 = interactionVolume.getNormal(idx1, 2);
        outerNormaln3 = interactionVolume.getNormal(idx1, 1);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][2];
        lambda13 = lambda[idx1][2];
        lambda31 = lambda[idx3][1];
        lambda15 = lambda[idx1][1];
        lambda51 = lambda[idx5][2];

        faceArea12 = interactionVolume.getFaceArea(idx1,0);
        faceArea21 = interactionVolume.getFaceArea(idx2,2);
        faceArea13 = interactionVolume.getFaceArea(idx1,2);
        faceArea31 = interactionVolume.getFaceArea(idx3,1);
        faceArea15 = interactionVolume.getFaceArea(idx1,1);
        faceArea51 = interactionVolume.getFaceArea(idx5,2);
    }
    else
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace13 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace15 = interactionVolume.getFacePosition(idx1,1);

        edgeCoord1213 = interactionVolume.getEdgePosition(idx1, 2, 1);
        edgeCoord1215 = interactionVolume.getEdgePosition(idx1, 2, 0);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 0, 1);

        outerNormaln1 = interactionVolume.getNormal(idx1, 2);
        outerNormaln2 = interactionVolume.getNormal(idx1, 0);
        outerNormaln3 = interactionVolume.getNormal(idx1, 1);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][0];
        lambda13 = lambda[idx1][0];
        lambda31 = lambda[idx3][1];
        lambda15 = lambda[idx1][1];
        lambda51 = lambda[idx5][0];

        faceArea12 = interactionVolume.getFaceArea(idx1,2);
        faceArea21 = interactionVolume.getFaceArea(idx2,0);
        faceArea13 = interactionVolume.getFaceArea(idx1,0);
        faceArea31 = interactionVolume.getFaceArea(idx3,1);
        faceArea15 = interactionVolume.getFaceArea(idx1,1);
        faceArea51 = interactionVolume.getFaceArea(idx5,0);
    }

    // compute normal vectors for case 1 (idx1, idx2, idx3, idx5)
    DimVector crossProductVector1(globalPosFace13-globalPos1);
    DimVector crossProductVector2(globalPosFace15-globalPos1);
    DimVector nu11 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace15-globalPos1;
    crossProductVector2 = globalPosFace12-globalPos1;
    DimVector nu12 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace12-globalPos1;
    crossProductVector2 = globalPosFace13-globalPos1;
    DimVector nu13 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = globalPosFace12-globalPos2;
    crossProductVector2 = edgeCoord1215-globalPos2;
    DimVector nu21 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1215-globalPos2;
    crossProductVector2 = edgeCoord1213-globalPos2;
    DimVector nu22 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1213-globalPos2;
    crossProductVector2 = globalPosFace12-globalPos2;
    DimVector nu23 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = edgeCoord1213-globalPos3;
    crossProductVector2 = edgeCoord1315-globalPos3;
    DimVector nu31 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1315-globalPos3;
    crossProductVector2 = globalPosFace13-globalPos3;
    DimVector nu32 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace13-globalPos3;
    crossProductVector2 = edgeCoord1213-globalPos3;
    DimVector nu33 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = edgeCoord1215-globalPos5;
    crossProductVector2 = globalPosFace15-globalPos5;
    DimVector nu51 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace15-globalPos5;
    crossProductVector2 = edgeCoord1315-globalPos5;
    DimVector nu52 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1315-globalPos5;
    crossProductVector2 = edgeCoord1215-globalPos5;
    DimVector nu53 = crossProduct(crossProductVector1, crossProductVector2);

    // compute T, i.e., the volume of the subtetrahedron
    crossProductVector1 = globalPosFace13-globalPos1;
    crossProductVector2 = globalPosFace15-globalPos1;
    Scalar T1 = (globalPosFace12-globalPos1) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1215-globalPos2;
    crossProductVector2 = edgeCoord1213-globalPos2;
    Scalar T2 = (globalPosFace12-globalPos2) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1213-globalPos3;
    crossProductVector2 = edgeCoord1315-globalPos3;
    Scalar T3 = (globalPosFace13-globalPos3) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1315-globalPos5;
    crossProductVector2 = edgeCoord1215-globalPos5;
    Scalar T5 = (globalPosFace15-globalPos5) * crossProduct(crossProductVector1, crossProductVector2);

    // compute components needed for flux calculation, denoted as 'omega' and 'r'
    DimVector K1nu11(0);
    K1.mv(nu11, K1nu11);
    DimVector K1nu12(0);
    K1.mv(nu12, K1nu12);
    DimVector K1nu13(0);
    K1.mv(nu13, K1nu13);

    DimVector K2nu21(0);
    K2.mv(nu21, K2nu21);
    DimVector K2nu22(0);
    K2.mv(nu22, K2nu22);
    DimVector K2nu23(0);
    K2.mv(nu23, K2nu23);

    DimVector K3nu31(0);
    K3.mv(nu31, K3nu31);
    DimVector K3nu32(0);
    K3.mv(nu32, K3nu32);
    DimVector K3nu33(0);
    K3.mv(nu33, K3nu33);

    DimVector K5nu51(0);
    K5.mv(nu51, K5nu51);
    DimVector K5nu52(0);
    K5.mv(nu52, K5nu52);
    DimVector K5nu53(0);
    K5.mv(nu53, K5nu53);

    Scalar omega111 = lambda12 * (outerNormaln1 * K1nu11) * faceArea12/T1;
    Scalar omega112 = lambda12 * (outerNormaln1 * K1nu12) * faceArea12/T1;
    Scalar omega113 = lambda12 * (outerNormaln1 * K1nu13) * faceArea12/T1;

    Scalar omega121 = lambda21 * (outerNormaln1 * K2nu21) * faceArea21/T2;
    Scalar omega122 = lambda21 * (outerNormaln1 * K2nu22) * faceArea21/T2;
    Scalar omega123 = lambda21 * (outerNormaln1 * K2nu23) * faceArea21/T2;

    Scalar omega211 = lambda13 * (outerNormaln2 * K1nu11) * faceArea13/T1;
    Scalar omega212 = lambda13 * (outerNormaln2 * K1nu12) * faceArea13/T1;
    Scalar omega213 = lambda13 * (outerNormaln2 * K1nu13) * faceArea13/T1;

    Scalar omega231 = lambda31 * (outerNormaln2 * K3nu31) * faceArea31/T3;
    Scalar omega232 = lambda31 * (outerNormaln2 * K3nu32) * faceArea31/T3;
    Scalar omega233 = lambda31 * (outerNormaln2 * K3nu33) * faceArea31/T3;

    Scalar omega311 = lambda15 * (outerNormaln3 * K1nu11) * faceArea15/T1;
    Scalar omega312 = lambda15 * (outerNormaln3 * K1nu12) * faceArea15/T1;
    Scalar omega313 = lambda15 * (outerNormaln3 * K1nu13) * faceArea15/T1;

    Scalar omega351 = lambda51 * (outerNormaln3 * K5nu51) * faceArea51/T5;
    Scalar omega352 = lambda51 * (outerNormaln3 * K5nu52) * faceArea51/T5;
    Scalar omega353 = lambda51 * (outerNormaln3 * K5nu53) * faceArea51/T5;

    Scalar r111 = (nu11 * (edgeCoord1213-globalPos1))/T1;
    Scalar r112 = (nu11 * (edgeCoord1215-globalPos1))/T1;
    Scalar r113 = (nu11 * (edgeCoord1315-globalPos1))/T1;

    Scalar r121 = (nu12 * (edgeCoord1213-globalPos1))/T1;
    Scalar r122 = (nu12 * (edgeCoord1215-globalPos1))/T1;
    Scalar r123 = (nu12 * (edgeCoord1315-globalPos1))/T1;

    Scalar r131 = (nu13 * (edgeCoord1213-globalPos1))/T1;
    Scalar r132 = (nu13 * (edgeCoord1215-globalPos1))/T1;
    Scalar r133 = (nu13 * (edgeCoord1315-globalPos1))/T1;

    // compute transmissibility matrix T = CA^{-1}B+D
    DimMatrix C(0), A(0);
    Dune::FieldMatrix<Scalar,dim,2*dim-dim+1> D(0), B(0);

    // evaluate matrix C, D, A, B
    C[0][0] = -omega111;
    C[0][1] = -omega112;
    C[0][2] = -omega113;
    C[1][0] = -omega211;
    C[1][1] = -omega212;
    C[1][2] = -omega213;
    C[2][0] = -omega311;
    C[2][1] = -omega312;
    C[2][2] = -omega313;

    D[0][0] = omega111 + omega112 + omega113;
    D[1][0] = omega211 + omega212 + omega213;
    D[2][0] = omega311 + omega312 + omega313;

    A[0][0] = omega111 - omega122 - omega121*r111 - omega123*r112;
    A[0][1] = omega112 - omega121*r121 - omega123*r122;
    A[0][2] = omega113 - omega121*r131 - omega123*r132;
    A[1][0] = omega211 - omega232*r111 - omega233*r113;
    A[1][1] = omega212 - omega231 - omega232*r121 - omega233*r123;
    A[1][2] = omega213 - omega232*r131 - omega233*r133;
    A[2][0] = omega311 - omega351*r113 - omega352*r112;
    A[2][1] = omega312 - omega351*r123 - omega352*r122;
    A[2][2] = omega313 - omega353 - omega351*r133 - omega352*r132;

    B[0][0] = omega111 + omega112 + omega113 + omega121*(1.0 - r111 - r121 -r131) + omega123*(1.0 - r112 - r122 - r132);
    B[0][1] = -omega121 - omega122 - omega123;
    B[1][0] = omega211 + omega212 + omega213 + omega232*(1.0 - r111 - r121 - r131) + omega233*(1.0 - r113 - r123 - r133);
    B[1][2] = -omega231 - omega232 - omega233;
    B[2][0] = omega311 + omega312 + omega313 + omega351*(1.0 - r113 - r123 - r133) + omega352*(1.0 - r112 -r122 -r132);
    B[2][3] = -omega351 - omega352 -omega353;

    // compute T
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));

    transmissibility = D;

    using std::isnan;
    if (isnan(transmissibility.frobenius_norm()))
    {
        std::cout<<"idx: "<<idx1<<idx2<<idx3<<idx5<<"\n";

        //        interactionVolume.printInteractionVolumeInfo();

        std::cout<<"case 1: transmissibility = "<<transmissibility<<"\n";
        std::cout<<"globalPos1 = "<<globalPos1<<"\n";
        std::cout<<"globalPos2 = "<<globalPos2<<"\n";
        std::cout<<"globalPos3 = "<<globalPos3<<"\n";
        std::cout<<"globalPos5 = "<<globalPos5<<"\n";
        std::cout<<"globalPosCenter = "<<globalPosCenter<<"\n";
        std::cout<<"outerNormaln1 = "<<outerNormaln1<<"\n";
        std::cout<<"outerNormaln2 = "<<outerNormaln2<<"\n";
        std::cout<<"outerNormaln3 = "<<outerNormaln3<<"\n";
        std::cout<<"xbar_1 = "<<globalPosFace12<<"\n";
        std::cout<<"xbar_2 = "<<globalPosFace13<<"\n";
        std::cout<<"xbar_3 = "<<globalPosFace15<<"\n";
        std::cout<<"xbar_4 = "<<edgeCoord1213<<"\n";
        std::cout<<"xbar_5 = "<<edgeCoord1215<<"\n";
        std::cout<<"xbar_6 = "<<edgeCoord1315<<"\n";
        std::cout<<"perm1 = "<<K1<<"\n";
        std::cout<<"perm2 = "<<K2<<"\n";
        std::cout<<"perm3 = "<<K3<<"\n";
        std::cout<<"perm5 = "<<K5<<"\n";
        std::cout<<"lambda = ";
        for (unsigned int i = 0; i < lambda.size(); i++)
        {
            std::cout<<lambda[i]<<" ";
        }
        std::cout<<"\n";
        std::cout<<"\n";
        DUNE_THROW(Dune::MathError,"T is nan");
    }

    return 1;
}

/*!
 * \brief Calculates the transmissibility matrix for the flux face between the cells with the local index idx1 and idx2 using the L-shape type 2
 *
 * For more details on L-shape type 2 see:
 *
 * Wolff, 2013 \cite wolff2013a <BR>
 * Wolff et al. 2013 \cite wolff2013b <BR>
 *
 *  \param transmissibility Reference to the local transmissibility matrix
 *  \param interactionVolume An interaction volume object
 *  \param lambda A vector containing the mobilities of the interaction volume cells in the order of the local interaction volume index
 *  \param idx1 Local index of cell 1 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx2 Local index of cell 2 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx4 Local index of cell 4 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa
 *  \param idx6 Local index of cell 6 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *
 *  \return 2 (L-shape 2)
 */
template<class TypeTag>
int FvMpfaL3dTransmissibilityCalculator<TypeTag>::transmissibilityCaseTwo(
                                                         Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                                         InteractionVolume& interactionVolume,
                                                         std::vector<DimVector >& lambda, int idx1, int idx2, int idx4, int idx6)
{
    auto element1 = interactionVolume.getSubVolumeElement(idx1);
    auto element2 = interactionVolume.getSubVolumeElement(idx2);
    auto element4 = interactionVolume.getSubVolumeElement(idx4);
    auto element6 = interactionVolume.getSubVolumeElement(idx6);

    // get global coordinate of cell centers
    const GlobalPosition& globalPos1 = element1.geometry().center();
    const GlobalPosition& globalPos2 = element2.geometry().center();
    const GlobalPosition& globalPos4 = element4.geometry().center();
    const GlobalPosition& globalPos6 = element6.geometry().center();

    const GlobalPosition& globalPosCenter = interactionVolume.getCenterPosition();

    const DimMatrix& K1 = problem_.spatialParams().intrinsicPermeability(element1);
    const DimMatrix& K2 = problem_.spatialParams().intrinsicPermeability(element2);
    const DimMatrix& K4 = problem_.spatialParams().intrinsicPermeability(element4);
    const DimMatrix& K6 = problem_.spatialParams().intrinsicPermeability(element6);

    // 1.use centered L-stencil 1 to compute the transmissibility of 1/4 face
    GlobalPosition globalPosFace12(0);
    GlobalPosition globalPosFace24(0);
    GlobalPosition globalPosFace26(0);

    DimVector outerNormaln1(0);
    DimVector outerNormaln2(0);
    DimVector outerNormaln3(0);

    GlobalPosition edgeCoord1224(0);
    GlobalPosition edgeCoord1226(0);
    GlobalPosition edgeCoord2426(0);

    // get lambda and face area
    Scalar lambda12 = 0;
    Scalar lambda21 = 0;
    Scalar lambda24 = 0;
    Scalar lambda42 = 0;
    Scalar lambda26 = 0;
    Scalar lambda62 = 0;
    Scalar faceArea12 = 0;
    Scalar faceArea21 = 0;
    Scalar faceArea24 = 0;
    Scalar faceArea42 = 0;
    Scalar faceArea26 = 0;
    Scalar faceArea62 = 0;

    if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 2) || (idx1 == 2 && idx2 == 0))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,0);

        outerNormaln1 = interactionVolume.getNormal(idx1, 0);

        globalPosFace24 = interactionVolume.getFacePosition(idx2,0);
        globalPosFace26 = interactionVolume.getFacePosition(idx2,2);

        edgeCoord1224 = interactionVolume.getEdgePosition(idx2, 1, 0);
        edgeCoord1226 = interactionVolume.getEdgePosition(idx2, 1, 1);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 0, 0);

        outerNormaln2 = interactionVolume.getNormal(idx2, 0);
        outerNormaln3 = interactionVolume.getNormal(idx2, 2);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][1];
        lambda24 = lambda[idx2][0];
        lambda42 = lambda[idx4][1];
        lambda26 = lambda[idx2][2];
        lambda62 = lambda[idx6][0];

        faceArea12 = interactionVolume.getFaceArea(idx1,0);
        faceArea21 = interactionVolume.getFaceArea(idx2,1);
        faceArea24 = interactionVolume.getFaceArea(idx2,0);
        faceArea42 = interactionVolume.getFaceArea(idx4,1);
        faceArea26 = interactionVolume.getFaceArea(idx2,2);
        faceArea62 = interactionVolume.getFaceArea(idx6,0);
    }
    else if ((idx1 == 5 && idx2 == 4) || (idx1 == 4 && idx2 == 6) || (idx1 == 6 && idx2 == 7) || (idx1 == 7 && idx2 == 5))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,2);

        outerNormaln1 = interactionVolume.getNormal(idx1, 2);

        globalPosFace24 = interactionVolume.getFacePosition(idx2,2);
        globalPosFace26 = interactionVolume.getFacePosition(idx2,0);

        edgeCoord1224 = interactionVolume.getEdgePosition(idx2, 1, 0);
        edgeCoord1226 = interactionVolume.getEdgePosition(idx2, 1, 1);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 2, 0);

        outerNormaln2 = interactionVolume.getNormal(idx2, 2);
        outerNormaln3 = interactionVolume.getNormal(idx2, 0);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][1];
        lambda24 = lambda[idx2][2];
        lambda42 = lambda[idx4][1];
        lambda26 = lambda[idx2][0];
        lambda62 = lambda[idx6][2];

        faceArea12 = interactionVolume.getFaceArea(idx1,2);
        faceArea21 = interactionVolume.getFaceArea(idx2,1);
        faceArea24 = interactionVolume.getFaceArea(idx2,2);
        faceArea42 = interactionVolume.getFaceArea(idx4,1);
        faceArea26 = interactionVolume.getFaceArea(idx2,0);
        faceArea62 = interactionVolume.getFaceArea(idx6,2);
    }
    else if ((idx1 == 4 && idx2 == 0) || (idx1 == 7 && idx2 == 3))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,0);

        outerNormaln1 = interactionVolume.getNormal(idx1, 0);
        outerNormaln2 = interactionVolume.getNormal(idx2, 1);
        outerNormaln3 = interactionVolume.getNormal(idx2, 0);

        globalPosFace24 = interactionVolume.getFacePosition(idx2,1);
        globalPosFace26 = interactionVolume.getFacePosition(idx2,0);

        edgeCoord1224 = interactionVolume.getEdgePosition(idx2, 2, 0);
        edgeCoord1226 = interactionVolume.getEdgePosition(idx2, 2, 1);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 1, 0);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][2];
        lambda24 = lambda[idx2][1];
        lambda42 = lambda[idx4][0];
        lambda26 = lambda[idx2][0];
        lambda62 = lambda[idx6][1];

        faceArea12 = interactionVolume.getFaceArea(idx1,0);
        faceArea21 = interactionVolume.getFaceArea(idx2,2);
        faceArea24 = interactionVolume.getFaceArea(idx2,1);
        faceArea42 = interactionVolume.getFaceArea(idx4,0);
        faceArea26 = interactionVolume.getFaceArea(idx2,0);
        faceArea62 = interactionVolume.getFaceArea(idx6,1);
    }
    else
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,2);

        outerNormaln1 = interactionVolume.getNormal(idx1, 2);

        globalPosFace24 = interactionVolume.getFacePosition(idx2,1);
        globalPosFace26 = interactionVolume.getFacePosition(idx2,2);

        edgeCoord1224 = interactionVolume.getEdgePosition(idx2, 0, 0);
        edgeCoord1226 = interactionVolume.getEdgePosition(idx2, 0, 1);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 1, 0);

        outerNormaln2 = interactionVolume.getNormal(idx2, 1);
        outerNormaln3 = interactionVolume.getNormal(idx2, 2);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][0];
        lambda24 = lambda[idx2][1];
        lambda42 = lambda[idx4][2];
        lambda26 = lambda[idx2][2];
        lambda62 = lambda[idx6][1];

        faceArea12 = interactionVolume.getFaceArea(idx1,2);
        faceArea21 = interactionVolume.getFaceArea(idx2,0);
        faceArea24 = interactionVolume.getFaceArea(idx2,1);
        faceArea42 = interactionVolume.getFaceArea(idx4,2);
        faceArea26 = interactionVolume.getFaceArea(idx2,2);
        faceArea62 = interactionVolume.getFaceArea(idx6,1);
    }


    // compute transmissibility matrix TC1 = CA^{-1}B+D
    DimMatrix C(0), A(0);
    Dune::FieldMatrix<Scalar,dim,2*dim-dim+1> D(0), B(0);

    // compute normal vectors for case 2 (idx1, idx2, idx4, idx6)
    DimVector crossProductVector1(edgeCoord1226-globalPos1);
    DimVector crossProductVector2(globalPosFace12-globalPos1);
    DimVector nu11 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1224-globalPos1;
    crossProductVector2 = edgeCoord1226-globalPos1;
    DimVector nu12 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace12-globalPos1;
    crossProductVector2 = edgeCoord1224-globalPos1;
    DimVector nu13 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = globalPosFace26-globalPos2;
    crossProductVector2 = globalPosFace24-globalPos2;
    DimVector nu21 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace12-globalPos2;
    crossProductVector2 = globalPosFace26-globalPos2;
    DimVector nu22 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace24-globalPos2;
    crossProductVector2 = globalPosFace12-globalPos2;
    DimVector nu23 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = edgeCoord2426-globalPos4;
    crossProductVector2 = edgeCoord1224-globalPos4;
    DimVector nu41 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace24-globalPos4;
    crossProductVector2 = edgeCoord2426-globalPos4;
    DimVector nu42 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1224-globalPos4;
    crossProductVector2 = globalPosFace24-globalPos4;
    DimVector nu43 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = globalPosFace26-globalPos6;
    crossProductVector2 = edgeCoord1226-globalPos6;
    DimVector nu61 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord2426-globalPos6;
    crossProductVector2 = globalPosFace26-globalPos6;
    DimVector nu62 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1226-globalPos6;
    crossProductVector2 = edgeCoord2426-globalPos6;
    DimVector nu63 = crossProduct(crossProductVector1, crossProductVector2);

    // compute T, i.e., the volume of the subtetrahedron
    crossProductVector1 = edgeCoord1224-globalPos1;
    crossProductVector2 = edgeCoord1226-globalPos1;
    Scalar T1 = (globalPosFace12-globalPos1) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace26-globalPos2;
    crossProductVector2 = globalPosFace24-globalPos2;
    Scalar T2 = (globalPosFace12-globalPos2) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord2426-globalPos4;
    crossProductVector2 = edgeCoord1224-globalPos4;
    Scalar T4 = (globalPosFace24-globalPos4) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1226-globalPos6;
    crossProductVector2 = edgeCoord2426-globalPos6;
    Scalar T6 = (globalPosFace26-globalPos6) * crossProduct(crossProductVector1, crossProductVector2);

    // compute components needed for flux calculation, denoted as 'omega' and 'r'
    DimVector K1nu11(0);
    K1.mv(nu11, K1nu11);
    DimVector K1nu12(0);
    K1.mv(nu12, K1nu12);
    DimVector K1nu13(0);
    K1.mv(nu13, K1nu13);

    DimVector K2nu21(0);
    K2.mv(nu21, K2nu21);
    DimVector K2nu22(0);
    K2.mv(nu22, K2nu22);
    DimVector K2nu23(0);
    K2.mv(nu23, K2nu23);

    DimVector K4nu41(0);
    K4.mv(nu41, K4nu41);
    DimVector K4nu42(0);
    K4.mv(nu42, K4nu42);
    DimVector K4nu43(0);
    K4.mv(nu43, K4nu43);

    DimVector K6nu61(0);
    K6.mv(nu61, K6nu61);
    DimVector K6nu62(0);
    K6.mv(nu62, K6nu62);
    DimVector K6nu63(0);
    K6.mv(nu63, K6nu63);

    Scalar omega111 = lambda12 * (outerNormaln1 * K1nu11) * faceArea12/T1;
    Scalar omega112 = lambda12 * (outerNormaln1 * K1nu12) * faceArea12/T1;
    Scalar omega113 = lambda12 * (outerNormaln1 * K1nu13) * faceArea12/T1;

    Scalar omega121 = lambda21 * (outerNormaln1 * K2nu21) * faceArea21/T2;
    Scalar omega122 = lambda21 * (outerNormaln1 * K2nu22) * faceArea21/T2;
    Scalar omega123 = lambda21 * (outerNormaln1 * K2nu23) * faceArea21/T2;

    Scalar omega221 = lambda24 * (outerNormaln2 * K2nu21) * faceArea24/T2;
    Scalar omega222 = lambda24 * (outerNormaln2 * K2nu22) * faceArea24/T2;
    Scalar omega223 = lambda24 * (outerNormaln2 * K2nu23) * faceArea24/T2;

    Scalar omega241 = lambda42 * (outerNormaln2 * K4nu41) * faceArea42/T4;
    Scalar omega242 = lambda42 * (outerNormaln2 * K4nu42) * faceArea42/T4;
    Scalar omega243 = lambda42 * (outerNormaln2 * K4nu43) * faceArea42/T4;

    Scalar omega321 = lambda26 * (outerNormaln3 * K2nu21) * faceArea26/T2;
    Scalar omega322 = lambda26 * (outerNormaln3 * K2nu22) * faceArea26/T2;
    Scalar omega323 = lambda26 * (outerNormaln3 * K2nu23) * faceArea26/T2;

    Scalar omega361 = lambda62 * (outerNormaln3 * K6nu61) * faceArea62/T6;
    Scalar omega362 = lambda62 * (outerNormaln3 * K6nu62) * faceArea62/T6;
    Scalar omega363 = lambda62 * (outerNormaln3 * K6nu63) * faceArea62/T6;

    Scalar r211 = (nu21 * (edgeCoord1224-globalPos2))/T2;
    Scalar r212 = (nu21 * (edgeCoord1226-globalPos2))/T2;
    Scalar r213 = (nu21 * (edgeCoord2426-globalPos2))/T2;

    Scalar r221 = (nu22 * (edgeCoord1224-globalPos2))/T2;
    Scalar r222 = (nu22 * (edgeCoord1226-globalPos2))/T2;
    Scalar r223 = (nu22 * (edgeCoord2426-globalPos2))/T2;

    Scalar r231 = (nu23 * (edgeCoord1224-globalPos2))/T2;
    Scalar r232 = (nu23 * (edgeCoord1226-globalPos2))/T2;
    Scalar r233 = (nu23 * (edgeCoord2426-globalPos2))/T2;

    // compute transmissibility matrix T = CA^{-1}B+D
    C = 0;
    A = 0;
    D = 0;
    B = 0;

    // evaluate matrix C, D, A, B
    C[0][0] = -omega121;
    C[0][1] = -omega122;
    C[0][2] = -omega123;
    C[1][0] = -omega221;
    C[1][1] = -omega222;
    C[1][2] = -omega223;
    C[2][0] = -omega321;
    C[2][1] = -omega322;
    C[2][2] = -omega323;

    D[0][1] = omega121 + omega122 + omega123;
    D[1][1] = omega221 + omega222 + omega223;
    D[2][1] = omega321 + omega322 + omega323;

    A[0][0] = omega121 - omega112 - omega111*r211 - omega113*r212;
    A[0][1] = omega122 - omega111*r221 - omega113*r222;
    A[0][2] = omega123 - omega111*r231 - omega113*r232;
    A[1][0] = omega221 - omega242*r211 - omega243*r213;
    A[1][1] = omega222 - omega241 - omega242*r221 - omega243*r223;
    A[1][2] = omega223 - omega242*r231 - omega243*r233;
    A[2][0] = omega321 - omega361*r213 - omega362*r212;
    A[2][1] = omega322 - omega361*r223 - omega362*r222;
    A[2][2] = omega323 - omega363 - omega361*r233 - omega362*r232;

    B[0][0] = -omega111 - omega112 - omega113;
    B[0][1] = omega121 + omega122 + omega123 + omega111*(1.0 - r211 - r221 -r231) + omega113*(1.0 - r212 - r222 - r232);
    B[1][1] = omega221 + omega222 + omega223 + omega242*(1.0 - r211 - r221 - r231) + omega243*(1.0 - r213 - r223 - r233);
    B[1][2] = -omega241 - omega242 - omega243;
    B[2][1] = omega321 + omega322 + omega323 + omega361*(1.0 - r213 - r223 - r233) + omega362*(1.0 - r212 -r222 -r232);
    B[2][3] = -omega361 - omega362 -omega363;

    // compute T
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));

    transmissibility = D;

    using std::isnan;
    if (isnan(transmissibility.frobenius_norm()))
    {
        std::cout<<"idx: "<<idx1<<idx2<<idx4<<idx6<<"\n";

        //        interactionVolume.printInteractionVolumeInfo();

        std::cout<<"case 2: transmissibility = "<<transmissibility<<"\n";
        std::cout<<"globalPos1 = "<<globalPos1<<"\n";
        std::cout<<"globalPos2 = "<<globalPos2<<"\n";
        std::cout<<"globalPos4 = "<<globalPos4<<"\n";
        std::cout<<"globalPos6 = "<<globalPos6<<"\n";
        std::cout<<"globalPosCenter = "<<globalPosCenter<<"\n";
        std::cout<<"outerNormaln1 = "<<outerNormaln1<<"\n";
        std::cout<<"outerNormaln2 = "<<outerNormaln2<<"\n";
        std::cout<<"outerNormaln3 = "<<outerNormaln3<<"\n";
        std::cout<<"xbar_1 = "<<globalPosFace12<<"\n";
        std::cout<<"xbar_2 = "<<globalPosFace24<<"\n";
        std::cout<<"xbar_3 = "<<globalPosFace26<<"\n";
        std::cout<<"xbar_6 = "<<edgeCoord2426<<"\n";
        std::cout<<"perm1 = "<<K1<<"\n";
        std::cout<<"perm2 = "<<K2<<"\n";
        std::cout<<"perm4 = "<<K4<<"\n";
        std::cout<<"perm6 = "<<K6<<"\n";
        std::cout<<"lambda = ";
        for (unsigned int i = 0; i < lambda.size(); i++)
        {
            std::cout<<lambda[i]<<" ";
        }
        std::cout<<"\n";
        std::cout<<"\n";
        DUNE_THROW(Dune::MathError,"T is nan");
    }

    return 2;
}

/*!
 * \brief Calculates the transmissibility matrix for the flux face between the cells with the local index idx1 and idx2 using the L-shape type 3
 *
 * For more details on L-shape type 3 see:
 *
 * Wolff, 2013 \cite wolff2013a <BR>
 * Wolff et al. 2013 \cite wolff2013b <BR>
 *
 *  \param transmissibility Reference to the local transmissibility matrix
 *  \param interactionVolume An interaction volume object
 *  \param lambda A vector containing the mobilities of the interaction volume cells in the order of the local interaction volume index
 *  \param idx1 Local index of cell 1 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx2 Local index of cell 2 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx4 Local index of cell 4 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa
 *  \param idx5 Local index of cell 5 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *
 *  \return 3 (L-shape 3)
 */
template<class TypeTag>
int FvMpfaL3dTransmissibilityCalculator<TypeTag>::transmissibilityCaseThree(
                                                         Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                                         InteractionVolume& interactionVolume,
                                                         std::vector<DimVector >& lambda,int idx1, int idx2, int idx4, int idx5)
{
    auto element1 = interactionVolume.getSubVolumeElement(idx1);
    auto element2 = interactionVolume.getSubVolumeElement(idx2);
    auto element4 = interactionVolume.getSubVolumeElement(idx4);
    auto element5 = interactionVolume.getSubVolumeElement(idx5);

    // get global coordinate of cell centers
    const GlobalPosition& globalPos1 = element1.geometry().center();
    const GlobalPosition& globalPos2 = element2.geometry().center();
    const GlobalPosition& globalPos4 = element4.geometry().center();
    const GlobalPosition& globalPos5 = element5.geometry().center();

    const GlobalPosition& globalPosCenter = interactionVolume.getCenterPosition();

    const DimMatrix& K1 = problem_.spatialParams().intrinsicPermeability(element1);
    const DimMatrix& K2 = problem_.spatialParams().intrinsicPermeability(element2);
    const DimMatrix& K4 = problem_.spatialParams().intrinsicPermeability(element4);
    const DimMatrix& K5 = problem_.spatialParams().intrinsicPermeability(element5);

    // 1.use centered L-stencil 1 to compute the transmissibility of 1/4 face
    GlobalPosition globalPosFace12(0);
    GlobalPosition globalPosFace15(0);
    GlobalPosition globalPosFace24(0);

    DimVector outerNormaln1(0);
    DimVector outerNormaln2(0);
    DimVector outerNormaln3(0);

    GlobalPosition edgeCoord1215(0);
    GlobalPosition edgeCoord1315(0);
    GlobalPosition edgeCoord1224(0);
    GlobalPosition edgeCoord2426(0);

    // get lambda and face area
    Scalar lambda12 = 0;
    Scalar lambda21 = 0;
    Scalar lambda24 = 0;
    Scalar lambda42 = 0;
    Scalar lambda15 = 0;
    Scalar lambda51 = 0;
    Scalar faceArea12 = 0;
    Scalar faceArea24 = 0;
    Scalar faceArea15 = 0;
    Scalar faceArea21 = 0;
    Scalar faceArea42 = 0;
    Scalar faceArea51 = 0;

    if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 2) || (idx1 == 2 && idx2 == 0))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace15 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace24 = interactionVolume.getFacePosition(idx2,0);

        outerNormaln1 = interactionVolume.getNormal(idx1, 0);
        outerNormaln2 = interactionVolume.getNormal(idx2, 0);
        outerNormaln3 = interactionVolume.getNormal(idx1, 2);

        edgeCoord1215 = interactionVolume.getEdgePosition(idx1, 0, 0);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 1, 1);
        edgeCoord1224 = interactionVolume.getEdgePosition(idx2, 1, 0);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 0, 0);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][1];
        lambda24 = lambda[idx2][0];
        lambda42 = lambda[idx4][1];
        lambda15 = lambda[idx1][2];
        lambda51 = lambda[idx5][0];

        faceArea12 = interactionVolume.getFaceArea(idx1,0);
        faceArea21 = interactionVolume.getFaceArea(idx2,1);
        faceArea24 = interactionVolume.getFaceArea(idx2,0);
        faceArea42 = interactionVolume.getFaceArea(idx4,1);
        faceArea15 = interactionVolume.getFaceArea(idx1,2);
        faceArea51 = interactionVolume.getFaceArea(idx5,0);
    }
    else if ((idx1 == 5 && idx2 == 4) || (idx1 == 4 && idx2 == 6) || (idx1 == 6 && idx2 == 7) || (idx1 == 7 && idx2 == 5))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace15 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace24 = interactionVolume.getFacePosition(idx2,2);

        outerNormaln1 = interactionVolume.getNormal(idx1, 2);
        outerNormaln2 = interactionVolume.getNormal(idx2, 2);
        outerNormaln3 = interactionVolume.getNormal(idx1, 0);

        edgeCoord1215 = interactionVolume.getEdgePosition(idx1, 2, 0);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 1, 1);
        edgeCoord1224 = interactionVolume.getEdgePosition(idx2, 1, 0);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 2, 0);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][1];
        lambda24 = lambda[idx2][2];
        lambda42 = lambda[idx4][1];
        lambda15 = lambda[idx1][0];
        lambda51 = lambda[idx5][2];

        faceArea12 = interactionVolume.getFaceArea(idx1,2);
        faceArea21 = interactionVolume.getFaceArea(idx2,1);
        faceArea24 = interactionVolume.getFaceArea(idx2,2);
        faceArea42 = interactionVolume.getFaceArea(idx4,1);
        faceArea15 = interactionVolume.getFaceArea(idx1,0);
        faceArea51 = interactionVolume.getFaceArea(idx5,2);
    }
    else if ((idx1 == 4 && idx2 == 0) || (idx1 == 7 && idx2 == 3))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace15 = interactionVolume.getFacePosition(idx1,1);
        globalPosFace24 = interactionVolume.getFacePosition(idx2,1);

        outerNormaln1 = interactionVolume.getNormal(idx1, 0);
        outerNormaln2 = interactionVolume.getNormal(idx2, 1);
        outerNormaln3 = interactionVolume.getNormal(idx1, 1);

        edgeCoord1215 = interactionVolume.getEdgePosition(idx1, 0, 0);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 2, 1);
        edgeCoord1224 = interactionVolume.getEdgePosition(idx2, 2, 0);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 1, 0);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][2];
        lambda24 = lambda[idx2][1];
        lambda42 = lambda[idx4][0];
        lambda15 = lambda[idx1][1];
        lambda51 = lambda[idx5][2];

        faceArea12 = interactionVolume.getFaceArea(idx1,0);
        faceArea21 = interactionVolume.getFaceArea(idx2,2);
        faceArea24 = interactionVolume.getFaceArea(idx2,1);
        faceArea42 = interactionVolume.getFaceArea(idx4,0);
        faceArea15 = interactionVolume.getFaceArea(idx1,1);
        faceArea51 = interactionVolume.getFaceArea(idx5,2);
    }
    else
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace15 = interactionVolume.getFacePosition(idx1,1);
        globalPosFace24 = interactionVolume.getFacePosition(idx2,1);

        outerNormaln1 = interactionVolume.getNormal(idx1, 2);
        outerNormaln2 = interactionVolume.getNormal(idx2, 1);
        outerNormaln3 = interactionVolume.getNormal(idx1, 1);

        edgeCoord1215 = interactionVolume.getEdgePosition(idx1, 2, 0);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 0, 1);
        edgeCoord1224 = interactionVolume.getEdgePosition(idx2, 0, 0);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 1, 0);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][0];
        lambda24 = lambda[idx2][1];
        lambda42 = lambda[idx4][2];
        lambda15 = lambda[idx1][1];
        lambda51 = lambda[idx5][0];

        faceArea12 = interactionVolume.getFaceArea(idx1,2);
        faceArea21 = interactionVolume.getFaceArea(idx2,0);
        faceArea24 = interactionVolume.getFaceArea(idx2,1);
        faceArea42 = interactionVolume.getFaceArea(idx4,2);
        faceArea15 = interactionVolume.getFaceArea(idx1,1);
        faceArea51 = interactionVolume.getFaceArea(idx5,0);
    }

    // 3.use non-centered L-stencil 1 to compute the transmissibility of 1/4 face
    // compute normal vectors for case 3 (idx1, idx2, idx4, idx5)
    DimVector crossProductVector1 = edgeCoord1224-globalPos1;
    DimVector crossProductVector2 = globalPosFace15-globalPos1;
    DimVector nu11 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace12-globalPos1;
    crossProductVector2 = edgeCoord1224-globalPos1;
    DimVector nu12 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace15-globalPos1;
    crossProductVector2 = globalPosFace12-globalPos1;
    DimVector nu13 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = edgeCoord1215-globalPos2;
    crossProductVector2 = globalPosFace24-globalPos2;
    DimVector nu21 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace12-globalPos2;
    crossProductVector2 = edgeCoord1215-globalPos2;
    DimVector nu22 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace24-globalPos2;
    crossProductVector2 = globalPosFace12-globalPos2;
    DimVector nu23 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = edgeCoord2426-globalPos4;
    crossProductVector2 = edgeCoord1224-globalPos4;
    DimVector nu41 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace24-globalPos4;
    crossProductVector2 = edgeCoord2426-globalPos4;
    DimVector nu42 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1224-globalPos4;
    crossProductVector2 = globalPosFace24-globalPos4;
    DimVector nu43 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = edgeCoord1315-globalPos5;
    crossProductVector2 = edgeCoord1215-globalPos5;
    DimVector nu51 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace15-globalPos5;
    crossProductVector2 = edgeCoord1315-globalPos5;
    DimVector nu52 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1215-globalPos5;
    crossProductVector2 = globalPosFace15-globalPos5;
    DimVector nu53 = crossProduct(crossProductVector1, crossProductVector2);

    // compute T, i.e., the volume of the subtetrahedron
    crossProductVector1 = edgeCoord1224-globalPos1;
    crossProductVector2 = globalPosFace15-globalPos1;
    Scalar T1 = (globalPosFace12-globalPos1) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1215-globalPos2;
    crossProductVector2 = globalPosFace24-globalPos2;
    Scalar T2 = (globalPosFace12-globalPos2) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord2426-globalPos4;
    crossProductVector2 = edgeCoord1224-globalPos4;
    Scalar T4 = (globalPosFace24-globalPos4) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1315-globalPos5;
    crossProductVector2 = edgeCoord1215-globalPos5;
    Scalar T5 = (globalPosFace15-globalPos5) * crossProduct(crossProductVector1, crossProductVector2);

    // compute components needed for flux calculation, denoted as 'omega' and 'r'
    DimVector K1nu11(0);
    K1.mv(nu11, K1nu11);
    DimVector K1nu12(0);
    K1.mv(nu12, K1nu12);
    DimVector K1nu13(0);
    K1.mv(nu13, K1nu13);

    DimVector K2nu21(0);
    K2.mv(nu21, K2nu21);
    DimVector K2nu22(0);
    K2.mv(nu22, K2nu22);
    DimVector K2nu23(0);
    K2.mv(nu23, K2nu23);

    DimVector K4nu41(0);
    K4.mv(nu41, K4nu41);
    DimVector K4nu42(0);
    K4.mv(nu42, K4nu42);
    DimVector K4nu43(0);
    K4.mv(nu43, K4nu43);

    DimVector K5nu51(0);
    K5.mv(nu51, K5nu51);
    DimVector K5nu52(0);
    K5.mv(nu52, K5nu52);
    DimVector K5nu53(0);
    K5.mv(nu53, K5nu53);


    Scalar omega111 = lambda12 * (outerNormaln1 * K1nu11) * faceArea12/T1;
    Scalar omega112 = lambda12 * (outerNormaln1 * K1nu12) * faceArea12/T1;
    Scalar omega113 = lambda12 * (outerNormaln1 * K1nu13) * faceArea12/T1;

    Scalar omega121 = lambda21 * (outerNormaln1 * K2nu21) * faceArea21/T2;
    Scalar omega122 = lambda21 * (outerNormaln1 * K2nu22) * faceArea21/T2;
    Scalar omega123 = lambda21 * (outerNormaln1 * K2nu23) * faceArea21/T2;

    Scalar omega221 = lambda24 * (outerNormaln2 * K2nu21) * faceArea24/T2;
    Scalar omega222 = lambda24 * (outerNormaln2 * K2nu22) * faceArea24/T2;
    Scalar omega223 = lambda24 * (outerNormaln2 * K2nu23) * faceArea24/T2;

    Scalar omega241 = lambda42 * (outerNormaln2 * K4nu41) * faceArea42/T4;
    Scalar omega242 = lambda42 * (outerNormaln2 * K4nu42) * faceArea42/T4;
    Scalar omega243 = lambda42 * (outerNormaln2 * K4nu43) * faceArea42/T4;

    Scalar omega311 = lambda15 * (outerNormaln3 * K1nu11) * faceArea15/T1;
    Scalar omega312 = lambda15 * (outerNormaln3 * K1nu12) * faceArea15/T1;
    Scalar omega313 = lambda15 * (outerNormaln3 * K1nu13) * faceArea15/T1;

    Scalar omega351 = lambda51 * (outerNormaln3 * K5nu51) * faceArea51/T5;
    Scalar omega352 = lambda51 * (outerNormaln3 * K5nu52) * faceArea51/T5;
    Scalar omega353 = lambda51 * (outerNormaln3 * K5nu53) * faceArea51/T5;

    Scalar r112 = (nu11 * (edgeCoord1215-globalPos1))/T1;
    Scalar r113 = (nu11 * (edgeCoord1315-globalPos1))/T1;
    Scalar r122 = (nu12 * (edgeCoord1215-globalPos1))/T1;
    Scalar r123 = (nu12 * (edgeCoord1315-globalPos1))/T1;
    Scalar r132 = (nu13 * (edgeCoord1215-globalPos1))/T1;
    Scalar r133 = (nu13 * (edgeCoord1315-globalPos1))/T1;
    Scalar r211 = (nu21 * (edgeCoord1224-globalPos2))/T2;
    Scalar r214 = (nu21 * (edgeCoord2426-globalPos2))/T2;
    Scalar r221 = (nu22 * (edgeCoord1224-globalPos2))/T2;
    Scalar r224 = (nu22 * (edgeCoord2426-globalPos2))/T2;
    Scalar r231 = (nu23 * (edgeCoord1224-globalPos2))/T2;
    Scalar r234 = (nu23 * (edgeCoord2426-globalPos2))/T2;

    Scalar coef = 1.0/(1.0 - r231 * r132);

    // compute transmissibility matrix T = CA^{-1}B+D
    DimMatrix C(0), A(0);
    Dune::FieldMatrix<Scalar,dim,2*dim-dim+1> D(0), B(0);

    // evaluate matrix C, D, A, B
    C[0][0] = -omega111 - coef * omega113 * r231 * r112 - coef * omega113 * r211;
    C[0][1] = -coef * omega113 * r221;
    C[0][2] = -omega112 - coef * omega113 * r231 * r122;
    C[1][0] = -omega221 - coef * omega223 * r112 - coef * omega223 * r132 * r211;
    C[1][1] = -omega222 - coef * omega223 * r132 * r221;
    C[1][2] = -coef * omega223 * r122;
    C[2][0] = -omega311 - coef * omega313 * r231 * r112 - coef * omega313 * r211;
    C[2][1] = -coef * omega313 * r221;
    C[2][2] = -omega312 - coef * omega313 * r231 * r122;

    D[0][0] = coef * omega113 * r231 * (r112 + r122 + r132 - 1.0) + omega111 + omega112 + omega113;
    D[0][1] = coef * omega113 * (r211 + r221 + r231 - 1.0);
    D[1][0] = coef * omega223 * (r112 + r122 + r132 - 1.0);
    D[1][1] = omega221 + omega222 + omega223 + coef * omega223 * r132 * (r211 + r221 + r231 - 1.0);
    D[2][0] = coef * omega313 * r231 * (r112 + r122 + r132 - 1.0) + omega311 + omega312 + omega313;
    D[2][1] = coef * omega313 * (r211 + r221 + r231 - 1.0);

    A[0][0] = omega111 - omega121 + coef * omega113 * (r231 * r112 + r211) - coef * omega123 * (r112 + r132 * r211);
    A[0][1] = coef * omega113 * r221 - omega122 - coef * omega123 * r132 * r221;
    A[0][2] = omega112 + coef * omega113 * r231 * r122 - coef * omega123 * r122;
    A[1][0] = omega221 - omega243 * r214 + coef * (omega223 - omega243 * r234) * (r112 + r132 * r211) - coef * omega242 * (r231 * r112 + r211);
    A[1][1] = omega222 - omega243 * r224 - omega241 + coef * (omega223 - omega243 * r234) * r132 * r221 - coef * omega242 * r221;
    A[1][2] = coef * r122 * (omega223 - omega243 * r234 - omega242 * r231);
    A[2][0] = omega311 - omega353 * r113 + coef * (omega313 - omega353 * r133) * (r231 * r112 + r211) - coef * omega352 * (r112 + r132 * r211);
    A[2][1] = coef * r221 * (omega313 - omega353 * r133 - omega352 * r132);
    A[2][2] = omega312 - omega351 - omega353 * r123 + coef * (omega313 - omega353 * r133) * r231 * r122 - coef * omega352 * r122;

    B[0][0] = coef * (omega113 * r231 - omega123) * (r112 + r122 + r132 - 1.0) + omega111 + omega112 + omega113;
    B[0][1] = coef * (omega113 - omega123 * r132) * (r211 + r221 + r231 - 1.0) - (omega121 + omega122 + omega123);
    B[1][0] = coef * (r112 + r122 + r132 - 1.0) * (omega223 - omega243 * r234 - omega242 * r231);
    B[1][1] = omega221 + omega222 + omega223 - omega243 * (r214 + r224 + r234 - 1.0) + coef * (r211 + r221 + r231 - 1.0)
        * (omega223 * r132 - omega243 * r234 * r132 - omega242);
    B[1][2] = -omega241 - omega242 - omega243;
    B[2][0] = coef * (r112 + r122 + r132 - 1.0) * (omega313 * r231 - omega353 * r133 * r231 - omega352) + omega311 + omega312 + omega313
        - omega353 * (r113 + r123 + r133 - 1.0);
    B[2][1] = coef * (r211 + r221 + r231 - 1.0) * (omega313 - omega353 * r133 - omega352 * r132);
    B[2][3] = -omega351 - omega352 - omega353;

    // compute T
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));

    transmissibility = D;

    if (std::isnan(transmissibility.frobenius_norm()))
    {
        std::cout<<"case 3: transmissibility = "<<transmissibility<<"\n";
        std::cout<<"globalPos1 = "<<globalPos1<<"\n";
        std::cout<<"globalPos2 = "<<globalPos2<<"\n";
        std::cout<<"globalPos4 = "<<globalPos4<<"\n";
        std::cout<<"globalPos5 = "<<globalPos5<<"\n";
        std::cout<<"globalPosCenter = "<<globalPosCenter<<"\n";
        std::cout<<"outerNormaln1 = "<<outerNormaln1<<"\n";
        std::cout<<"outerNormaln2 = "<<outerNormaln2<<"\n";
        std::cout<<"outerNormaln3 = "<<outerNormaln3<<"\n";
        std::cout<<"xbar_1 = "<<globalPosFace12<<"\n";
        std::cout<<"xbar_2 = "<<globalPosFace24<<"\n";
        std::cout<<"xbar_3 = "<<globalPosFace15<<"\n";
        std::cout<<"xbar_5 = "<<edgeCoord1215<<"\n";
        std::cout<<"xbar_6 = "<<edgeCoord1315<<"\n";
        std::cout<<"xbar_7 = "<<edgeCoord2426<<"\n";
        std::cout<<"perm1 = "<<K1<<"\n";
        std::cout<<"perm2 = "<<K2<<"\n";
        std::cout<<"perm4 = "<<K4<<"\n";
        std::cout<<"perm5 = "<<K5<<"\n";
        std::cout<<"lambda = ";
        for (unsigned int i = 0; i < lambda.size(); i++)
        {
            std::cout<<lambda[i]<<" ";
        }
        std::cout<<"\n";
        std::cout<<"\n";
        DUNE_THROW(Dune::MathError,"T is nan");
    }

    return 3;
}

/*!
 * \brief Calculates the transmissibility matrix for the flux face between the cells with the local index idx1 and idx2 using the L-shape type 4
 *
 * For more details on L-shape type 4 see:
 *
 * Wolff, 2013 \cite wolff2013a <BR>
 * Wolff et al. 2013 \cite wolff2013b <BR>
 *
 *  \param transmissibility Reference to the local transmissibility matrix
 *  \param interactionVolume An interaction volume object
 *  \param lambda A vector containing the mobilities of the interaction volume cells in the order of the local interaction volume index
 *  \param idx1 Local index of cell 1 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx2 Local index of cell 2 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *  \param idx3 Local index of cell 3 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa
 *  \param idx6 Local index of cell 6 of 6 cells for the calculation of the flux between cell idx1 and idx2  (see doc/docextra/3dmpfa)
 *
 *  \return 4 (L-shape 4)
 */
template<class TypeTag>
int FvMpfaL3dTransmissibilityCalculator<TypeTag>::transmissibilityCaseFour(
                                                         Dune::FieldMatrix<Scalar,dim,2*dim-dim+1>& transmissibility,
                                                         InteractionVolume& interactionVolume,
                                                         std::vector<DimVector >& lambda,int idx1, int idx2, int idx3, int idx6)
{
    auto element1 = interactionVolume.getSubVolumeElement(idx1);
    auto element2 = interactionVolume.getSubVolumeElement(idx2);
    auto element3 = interactionVolume.getSubVolumeElement(idx3);
    auto element6 = interactionVolume.getSubVolumeElement(idx6);

    // get global coordinate of cell centers
    const GlobalPosition& globalPos1 = element1.geometry().center();
    const GlobalPosition& globalPos2 = element2.geometry().center();
    const GlobalPosition& globalPos3 = element3.geometry().center();
    const GlobalPosition& globalPos6 = element6.geometry().center();

    const GlobalPosition& globalPosCenter = interactionVolume.getCenterPosition();

    const DimMatrix& K1 = problem_.spatialParams().intrinsicPermeability(element1);
    const DimMatrix& K2 = problem_.spatialParams().intrinsicPermeability(element2);
    const DimMatrix& K3 = problem_.spatialParams().intrinsicPermeability(element3);
    const DimMatrix& K6 = problem_.spatialParams().intrinsicPermeability(element6);

    // 1.use centered L-stencil 1 to compute the transmissibility of 1/4 face
    GlobalPosition globalPosFace12(0);
    GlobalPosition globalPosFace13(0);
    GlobalPosition globalPosFace26(0);

    DimVector outerNormaln1(0);
    DimVector outerNormaln2(0);
    DimVector outerNormaln3(0);

    GlobalPosition edgeCoord1213(0);
    GlobalPosition edgeCoord1315(0);
    GlobalPosition edgeCoord1226(0);
    GlobalPosition edgeCoord2426(0);

    // get lambda and face area
    Scalar lambda12 = 0;
    Scalar lambda21 = 0;
    Scalar lambda13 = 0;
    Scalar lambda31 = 0;
    Scalar lambda26 = 0;
    Scalar lambda62 = 0;
    Scalar faceArea12 = 0;
    Scalar faceArea21 = 0;
    Scalar faceArea13 = 0;
    Scalar faceArea31 = 0;
    Scalar faceArea26 = 0;
    Scalar faceArea62 = 0;

    if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 3) || (idx1 == 3 && idx2 == 2) || (idx1 == 2 && idx2 == 0))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace13 = interactionVolume.getFacePosition(idx1,1);
        globalPosFace26 = interactionVolume.getFacePosition(idx2,2);

        outerNormaln1 = interactionVolume.getNormal(idx1, 0);
        outerNormaln2 = interactionVolume.getNormal(idx1, 1);
        outerNormaln3 = interactionVolume.getNormal(idx2, 2);

        edgeCoord1213 = interactionVolume.getEdgePosition(idx1, 0, 1);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 1, 1);
        edgeCoord1226 = interactionVolume.getEdgePosition(idx2, 1, 1);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 0, 0);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][1];
        lambda13 = lambda[idx1][1];
        lambda31 = lambda[idx3][0];
        lambda26 = lambda[idx2][2];
        lambda62 = lambda[idx6][0];

        faceArea12 = interactionVolume.getFaceArea(idx1,0);
        faceArea21 = interactionVolume.getFaceArea(idx2,1);
        faceArea13 = interactionVolume.getFaceArea(idx1,1);
        faceArea31 = interactionVolume.getFaceArea(idx3,0);
        faceArea26 = interactionVolume.getFaceArea(idx2,2);
        faceArea62 = interactionVolume.getFaceArea(idx6,0);
    }
    else if ((idx1 == 5 && idx2 == 4) || (idx1 == 4 && idx2 == 6) || (idx1 == 6 && idx2 == 7) || (idx1 == 7 && idx2 == 5))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace13 = interactionVolume.getFacePosition(idx1,1);
        globalPosFace26 = interactionVolume.getFacePosition(idx2,0);

        outerNormaln1 = interactionVolume.getNormal(idx1, 2);
        outerNormaln2 = interactionVolume.getNormal(idx1, 1);
        outerNormaln3 = interactionVolume.getNormal(idx2, 0);

        edgeCoord1213 = interactionVolume.getEdgePosition(idx1, 2, 1);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 1, 1);
        edgeCoord1226 = interactionVolume.getEdgePosition(idx2, 1, 1);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 2, 0);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][1];
        lambda13 = lambda[idx1][1];
        lambda31 = lambda[idx3][2];
        lambda26 = lambda[idx2][0];
        lambda62 = lambda[idx6][2];

        faceArea12 = interactionVolume.getFaceArea(idx1,2);
        faceArea21 = interactionVolume.getFaceArea(idx2,1);
        faceArea13 = interactionVolume.getFaceArea(idx1,1);
        faceArea31 = interactionVolume.getFaceArea(idx3,2);
        faceArea26 = interactionVolume.getFaceArea(idx2,0);
        faceArea62 = interactionVolume.getFaceArea(idx6,2);
    }
    else if ((idx1 == 4 && idx2 == 0) || (idx1 == 7 && idx2 == 3))
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace13 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace26 = interactionVolume.getFacePosition(idx2,0);

        outerNormaln1 = interactionVolume.getNormal(idx1, 0);
        outerNormaln2 = interactionVolume.getNormal(idx1, 2);
        outerNormaln3 = interactionVolume.getNormal(idx2, 0);

        edgeCoord1213 = interactionVolume.getEdgePosition(idx1, 0, 1);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 2, 1);
        edgeCoord1226 = interactionVolume.getEdgePosition(idx2, 2, 1);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 1, 0);

        lambda12 = lambda[idx1][0];
        lambda21 = lambda[idx2][2];
        lambda13 = lambda[idx1][2];
        lambda31 = lambda[idx3][1];
        lambda26 = lambda[idx2][0];
        lambda62 = lambda[idx6][1];

        faceArea12 = interactionVolume.getFaceArea(idx1,0);
        faceArea21 = interactionVolume.getFaceArea(idx2,2);
        faceArea13 = interactionVolume.getFaceArea(idx1,2);
        faceArea31 = interactionVolume.getFaceArea(idx3,1);
        faceArea26 = interactionVolume.getFaceArea(idx2,0);
        faceArea62 = interactionVolume.getFaceArea(idx6,1);
    }
    else
    {
        globalPosFace12 = interactionVolume.getFacePosition(idx1,2);
        globalPosFace13 = interactionVolume.getFacePosition(idx1,0);
        globalPosFace26 = interactionVolume.getFacePosition(idx2,2);

        outerNormaln1 = interactionVolume.getNormal(idx1, 2);
        outerNormaln2 = interactionVolume.getNormal(idx1, 0);
        outerNormaln3 = interactionVolume.getNormal(idx2, 2);

        edgeCoord1213 = interactionVolume.getEdgePosition(idx1, 2, 1);
        edgeCoord1315 = interactionVolume.getEdgePosition(idx1, 0, 1);
        edgeCoord1226 = interactionVolume.getEdgePosition(idx2, 0, 1);
        edgeCoord2426 = interactionVolume.getEdgePosition(idx2, 1, 0);

        lambda12 = lambda[idx1][2];
        lambda21 = lambda[idx2][0];
        lambda13 = lambda[idx1][0];
        lambda31 = lambda[idx3][1];
        lambda26 = lambda[idx2][2];
        lambda62 = lambda[idx6][1];

        faceArea12 = interactionVolume.getFaceArea(idx1,2);
        faceArea21 = interactionVolume.getFaceArea(idx2,0);
        faceArea13 = interactionVolume.getFaceArea(idx1,0);
        faceArea31 = interactionVolume.getFaceArea(idx3,1);
        faceArea26 = interactionVolume.getFaceArea(idx2,2);
        faceArea62 = interactionVolume.getFaceArea(idx6,1);
    }

    // 4.use non-centered L-stencil 2 to compute the transmissibility of 1/4 face
    // compute normal vectors for case 4 (idx1, idx2, idx3, idx6)
    DimVector crossProductVector1 = globalPosFace13-globalPos1;
    DimVector crossProductVector2 = edgeCoord1226-globalPos1;
    DimVector nu11 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1226-globalPos1;
    crossProductVector2 = globalPosFace12-globalPos1;
    DimVector nu12 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace12-globalPos1;
    crossProductVector2 = globalPosFace13-globalPos1;
    DimVector nu13 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = globalPosFace26-globalPos2;
    crossProductVector2 = edgeCoord1213-globalPos2;
    DimVector nu21 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1213-globalPos2;
    crossProductVector2 = globalPosFace12-globalPos2;
    DimVector nu22 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace12-globalPos2;
    crossProductVector2 = globalPosFace26-globalPos2;
    DimVector nu23 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = edgeCoord1213-globalPos3;
    crossProductVector2 = edgeCoord1315-globalPos3;
    DimVector nu31 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1315-globalPos3;
    crossProductVector2 = globalPosFace13-globalPos3;
    DimVector nu32 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace13-globalPos3;
    crossProductVector2 = edgeCoord1213-globalPos3;
    DimVector nu33 = crossProduct(crossProductVector1, crossProductVector2);

    crossProductVector1 = edgeCoord1226-globalPos6;
    crossProductVector2 = edgeCoord2426-globalPos6;
    DimVector nu61 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord2426-globalPos6;
    crossProductVector2 = globalPosFace26-globalPos6;
    DimVector nu62 = crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace26-globalPos6;
    crossProductVector2 = edgeCoord1226-globalPos6;
    DimVector nu63 = crossProduct(crossProductVector1, crossProductVector2);

    // compute T, i.e., the volume of the subtetrahedron
    crossProductVector1 = globalPosFace13-globalPos1;
    crossProductVector2 = edgeCoord1226-globalPos1;
    Scalar T1 = (globalPosFace12-globalPos1) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = globalPosFace26-globalPos2;
    crossProductVector2 = edgeCoord1213-globalPos2;
    Scalar T2 = (globalPosFace12-globalPos2) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1213-globalPos3;
    crossProductVector2 = edgeCoord1315-globalPos3;
    Scalar T3 = (globalPosFace13-globalPos3) * crossProduct(crossProductVector1, crossProductVector2);
    crossProductVector1 = edgeCoord1226-globalPos6;
    crossProductVector2 = edgeCoord2426-globalPos6;
    Scalar T6 = (globalPosFace26-globalPos6) * crossProduct(crossProductVector1, crossProductVector2);

    // compute components needed for flux calculation, denoted as 'omega' and 'r'
    DimVector K1nu11(0);
    K1.mv(nu11, K1nu11);
    DimVector K1nu12(0);
    K1.mv(nu12, K1nu12);
    DimVector K1nu13(0);
    K1.mv(nu13, K1nu13);

    DimVector K2nu21(0);
    K2.mv(nu21, K2nu21);
    DimVector K2nu22(0);
    K2.mv(nu22, K2nu22);
    DimVector K2nu23(0);
    K2.mv(nu23, K2nu23);

    DimVector K3nu31(0);
    K3.mv(nu31, K3nu31);
    DimVector K3nu32(0);
    K3.mv(nu32, K3nu32);
    DimVector K3nu33(0);
    K3.mv(nu33, K3nu33);

    DimVector K6nu61(0);
    K6.mv(nu61, K6nu61);
    DimVector K6nu62(0);
    K6.mv(nu62, K6nu62);
    DimVector K6nu63(0);
    K6.mv(nu63, K6nu63);

    //    std::cout<<"outerNormaln2 = "<<outerNormaln2<<"\n";
    //    std::cout<<"outerNormaln3 = "<<outerNormaln3<<"\n";

    Scalar omega111 = lambda12 * (outerNormaln1 * K1nu11) * faceArea12/T1;
    Scalar omega112 = lambda12 * (outerNormaln1 * K1nu12) * faceArea12/T1;
    Scalar omega113 = lambda12 * (outerNormaln1 * K1nu13) * faceArea12/T1;

    Scalar omega121 = lambda21 * (outerNormaln1 * K2nu21) * faceArea21/T2;
    Scalar omega122 = lambda21 * (outerNormaln1 * K2nu22) * faceArea21/T2;
    Scalar omega123 = lambda21 * (outerNormaln1 * K2nu23) * faceArea21/T2;

    Scalar omega211 = lambda13 * (outerNormaln2 * K1nu11) * faceArea13/T1;
    Scalar omega212 = lambda13 * (outerNormaln2 * K1nu12) * faceArea13/T1;
    Scalar omega213 = lambda13 * (outerNormaln2 * K1nu13) * faceArea13/T1;

    Scalar omega231 = lambda31 * (outerNormaln2 * K3nu31) * faceArea31/T3;
    Scalar omega232 = lambda31 * (outerNormaln2 * K3nu32) * faceArea31/T3;
    Scalar omega233 = lambda31 * (outerNormaln2 * K3nu33) * faceArea31/T3;

    Scalar omega321 = lambda26 * (outerNormaln3 * K2nu21) * faceArea26/T2;
    Scalar omega322 = lambda26 * (outerNormaln3 * K2nu22) * faceArea26/T2;
    Scalar omega323 = lambda26 * (outerNormaln3 * K2nu23) * faceArea26/T2;

    Scalar omega361 = lambda62 * (outerNormaln3 * K6nu61) * faceArea62/T6;
    Scalar omega362 = lambda62 * (outerNormaln3 * K6nu62) * faceArea62/T6;
    Scalar omega363 = lambda62 * (outerNormaln3 * K6nu63) * faceArea62/T6;

    Scalar r111 = (nu11 * (edgeCoord1213-globalPos1))/T1;
    Scalar r114 = (nu11 * (edgeCoord1315-globalPos1))/T1;
    Scalar r121 = (nu12 * (edgeCoord1213-globalPos1))/T1;
    Scalar r124 = (nu12 * (edgeCoord1315-globalPos1))/T1;
    Scalar r131 = (nu13 * (edgeCoord1213-globalPos1))/T1;
    Scalar r134 = (nu13 * (edgeCoord1315-globalPos1))/T1;

    Scalar r212 = (nu21 * (edgeCoord1226-globalPos2))/T2;
    Scalar r213 = (nu21 * (edgeCoord2426-globalPos2))/T2;
    Scalar r222 = (nu22 * (edgeCoord1226-globalPos2))/T2;
    Scalar r223 = (nu22 * (edgeCoord2426-globalPos2))/T2;
    Scalar r232 = (nu23 * (edgeCoord1226-globalPos2))/T2;
    Scalar r233 = (nu23 * (edgeCoord2426-globalPos2))/T2;

    Scalar coef = 1.0/(1.0 - r232 * r131);

    // compute transmissibility matrix TC1 = CA^{-1}B+D
    DimMatrix C(0), A(0);
    Dune::FieldMatrix<Scalar,dim,2*dim-dim+1> D(0), B(0);

    // evaluate matrix C, D, A, B
    C[0][0] = -omega111 - coef * omega113 * r212 - coef * omega113 * r232 * r111;
    C[0][1] = -omega112 - coef * omega113 * r232 * r121;
    C[0][2] = -coef * omega113 * r222;
    C[1][0] = -omega211 - coef * omega213 * r212 - coef * omega213 * r232 * r111;
    C[1][1] = -omega212 - coef * omega213 * r232 * r121;
    C[1][2] = -coef * omega213 * r222;
    C[2][0] = -omega321 - coef * omega323 * r111 - coef * omega323 * r131 * r212;
    C[2][1] = -coef * omega323 * r121;
    C[2][2] = -omega322 - coef * omega323 * r131 * r222;

    D[0][0] = coef * omega113 * r232 * (r111 + r121 + r131 - 1.0) + omega111 + omega112 + omega113;
    D[0][1] = coef * omega113 * (r212 + r222 + r232 - 1.0);
    D[1][0] = coef * omega213 * r232 * (r111 + r121 + r131 - 1.0) + omega211 + omega212 + omega213;
    D[1][1] = coef * omega213 * (r212 + r222 + r232 - 1.0);
    D[2][0] = coef * omega323 * (r111 + r121 + r131 - 1.0);
    D[2][1] = coef * omega323 * r131 * (r212 + r222 + r232 - 1.0) + omega321 + omega322 + omega323;

    A[0][0] = omega111 - omega121 + coef * omega113 * (r212 + r232 * r111) - coef * omega123 * (r111 + r131 * r212);
    A[0][1] = omega112 + coef * omega113 * r232 * r121 - coef * omega123 * r121;
    A[0][2] = coef * omega113 * r222 - omega122 - coef * omega123 * r131 * r222;
    A[1][0] = omega211 - omega233 * r114 + coef * (omega213 - omega233 * r134) * (r212 + r232 * r111) - coef * omega232 * (r111 + r131 * r212);
    A[1][1] = omega212 - omega231 - omega233 * r124 + coef * (omega213 - omega233 * r134) * r232 * r121 - coef * omega232 * r121;
    A[1][2] = coef * r222 * (omega213 - omega233 * r134 - omega232 * r131);
    A[2][0] = omega321 - omega363 * r213 + coef * (omega323 - omega363 * r233) * (r111 + r131 * r212) - coef * omega362 * (r212 + r232 * r111);
    A[2][1] = coef * r121 * (omega323 - omega363 * r233 - omega362 * r232);
    A[2][2] = omega322 - omega361 - omega363 * r223 + coef * (omega323 - omega363 * r233) * r131 * r222 - coef * omega362 * r222;

    B[0][0] = coef * (omega113 * r232 - omega123) * (r111 + r121 + r131 - 1.0) + omega111 + omega112 + omega113;
    B[0][1] = coef * (omega113 - omega123 * r131) * (r212 + r222 + r232 - 1.0) - (omega121 + omega122 + omega123);
    B[1][0] = omega211 + omega212 + omega213 - omega233 * (r114 + r124 + r134 - 1.0) + coef * (r111 + r121 + r131 - 1.0)
        * (omega213 * r232 - omega233 * r134 * r232 - omega232);
    B[1][1] = coef * (r212 + r222 + r232 - 1.0) * (omega213 - omega233 * r134 - omega232 * r131);
    B[1][2] = -omega231 - omega232 - omega233;
    B[2][0] = coef * (r111 + r121 + r131 - 1.0) * (omega323 - omega363 * r233 - omega362 * r232);
    B[2][1] = omega321 + omega322 + omega323 + coef * (r212 + r222 + r232 - 1.0) * (omega323 * r131 - omega363 * r233 * r131 - omega362)
        - omega363 * (r213 + r223 + r233 - 1.0);
    B[2][3] = -omega361 - omega362 - omega363;

    // compute T
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));

    transmissibility = D;

    if (std::isnan(transmissibility.frobenius_norm()))
    {
        std::cout<<"case 4: transmissibility = "<<transmissibility<<"\n";

        std::cout<<"globalPos1 = "<<globalPos1<<"\n";
        std::cout<<"globalPos2 = "<<globalPos2<<"\n";
        std::cout<<"globalPos3 = "<<globalPos3<<"\n";
        std::cout<<"globalPos6 = "<<globalPos6<<"\n";

        std::cout<<"globalPosCenter = "<<globalPosCenter<<"\n";
        std::cout<<"outerNormaln1 = "<<outerNormaln1<<"\n";
        std::cout<<"outerNormaln2 = "<<outerNormaln2<<"\n";
        std::cout<<"outerNormaln3 = "<<outerNormaln3<<"\n";
        std::cout<<"xbar_1 = "<<globalPosFace12<<"\n";
        std::cout<<"xbar_2 = "<<globalPosFace13<<"\n";
        std::cout<<"xbar_3 = "<<globalPosFace26<<"\n";
        std::cout<<"xbar_4 = "<<edgeCoord1213<<"\n";
        std::cout<<"xbar_5 = "<<edgeCoord1226<<"\n";
        std::cout<<"xbar_6 = "<<edgeCoord2426<<"\n";
        std::cout<<"xbar_7 = "<<edgeCoord1315<<"\n";
        std::cout<<"perm1 = "<<K1<<"\n";
        std::cout<<"perm2 = "<<K2<<"\n";
        std::cout<<"perm3 = "<<K3<<"\n";
        std::cout<<"perm6 = "<<K6<<"\n";
        std::cout<<"lambda = ";
        for (unsigned int i = 0; i < lambda.size(); i++)
        {
            std::cout<<lambda[i]<<" ";
        }
        std::cout<<"\n";
        std::cout<<"\n";
        DUNE_THROW(Dune::MathError,"T is nan");
    }

    return 4;
}
} // end namespace Dumux
#endif
