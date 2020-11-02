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
 * \ingroup SequentialTwoPModel
 * \brief  Velocity calculation on adaptive grids using a 2-d MPFA L-method.
 */
#ifndef DUMUX_FVMPFAL2DVELOCITY2P_ADAPTIVE_HH
#define DUMUX_FVMPFAL2DVELOCITY2P_ADAPTIVE_HH

#include "2dvelocity.hh"

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Class for calculating 2-d velocities from cell-wise constant pressure values.
 *
 * Calculates phase velocities or total velocity from a known pressure field applying a
 * finite volume discretization and a MPFA L-method.
 * At Dirichlet boundaries a two-point flux approximation is used.
 * The pressure has to be given as piecewise constant cell values.
 * The velocities are calculated as
 *
 *\f[ \boldsymbol v_\alpha = - \lambda_\alpha \boldsymbol K \textbf{grad}\, \Phi_\alpha, \f]
 * and,
 * \f[ \boldsymbol v_t = \boldsymbol v_w + \boldsymbol v_n,\f]
 *
 * where \f$ \Phi_\alpha \f$ denotes the potential of phase \f$ \alpha \f$,
 * \f$ \boldsymbol K \f$ the intrinsic permeability,
 * and \f$ \lambda_\alpha \f$ a phase mobility.
 *
 * Remark1: only for quadrilaterals!
 *
 * \tparam TypeTag The problem Type Tag
 */
template<class TypeTag> class FvMpfaL2dVelocity2pAdaptive : public FvMpfaL2dVelocity2p<TypeTag>
{
    using ParentType = FvMpfaL2dVelocity2p<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Grid = typename GridView::Grid;
    using IndexSet = typename GridView::IndexSet;

    using Geometry = typename Element::Geometry;
    using JacobianTransposed = typename Geometry::JacobianTransposed ;

    using GridTypeIndices = GetPropType<TypeTag, Properties::GridTypeIndices>;

    using InteractionVolume = FVMPFALInteractionVolume<TypeTag>;
    using InnerBoundaryVolumeFaces = std::vector<Dune::FieldVector<bool, 2*dim> >;
    using TransmissibilityCalculator = FvMpfaL2dTransmissibilityCalculator<TypeTag>;

    enum
        {
            pw = Indices::pressureW,
            pn = Indices::pressureNw,
            pGlobal = Indices::pressureGlobal,
            sw = Indices::saturationW,
            sn = Indices::saturationNw,
            vw = Indices::velocityW,
            vn = Indices::velocityNw,
            vt = Indices::velocityTotal
        };
    enum
        {
            wPhaseIdx = Indices::wPhaseIdx,
            nPhaseIdx = Indices::nPhaseIdx,
            pressureIdx = Indices::pressureIdx,
            saturationIdx = Indices::saturationIdx,
            pressureEqIdx = Indices::pressureEqIdx,
            satEqIdx = Indices::satEqIdx,
            numPhases = getPropValue<TypeTag, Properties::NumPhases>()
        };

    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using DimVector = Dune::FieldVector<Scalar, dim>;

public:
    /*!
     * \brief Constructs a FvMpfaL2dVelocity2pAdaptive object
     *
     * \param problem A problem class object
     */
    FvMpfaL2dVelocity2pAdaptive(Problem& problem) :
        ParentType(problem), problem_(problem)
    {}

    //! Calculatex velocities for flux faces of a hanging node interaction volume
    void calculateHangingNodeInteractionVolumeVelocity(InteractionVolume& interactionVolume, CellData& cellData1,
                                                       CellData& cellData2, CellData& cellData4,
                                                       InnerBoundaryVolumeFaces& innerBoundaryVolumeFaces);

private:
    Problem& problem_;
};
// end of template

/*!
 * \brief Calculates the velocities at the flux faces of an interation volume around a hanging node vertex.
 *
 *  Calculates the velocities at the flux faces of an interation volume around a hanging
 *  node vertex and adds them to the face velocity vectors in the <tt>CellData</tt> objects.
 *
 * \param interactionVolume An <tt>InteractionVolume</tt> object including the information for calculating the MPFA transmissibilities
 * \param cellData1  <tt>CellData</tt> object of an IMPES model for sub-volume 1
 * \param cellData2  <tt>CellData</tt> object of an IMPES model for sub-volume 2
 * \param cellData4  <tt>CellData</tt> object of an IMPES model for sub-volume 4
 * \param innerBoundaryVolumeFaces container including information about faces intersecting a boundary
 */
template<class TypeTag>
void FvMpfaL2dVelocity2pAdaptive<TypeTag>::calculateHangingNodeInteractionVolumeVelocity(InteractionVolume& interactionVolume,
                                                                                         CellData& cellData1, CellData& cellData2,
                                                                                         CellData& cellData4,
                                                                                         InnerBoundaryVolumeFaces& innerBoundaryVolumeFaces)
{
    // cell index
    int globalIdx1 = problem_.variables().index(interactionVolume.getSubVolumeElement(0));
    int globalIdx2 = problem_.variables().index(interactionVolume.getSubVolumeElement(1));

    // get pressure values
    Dune::FieldVector < Scalar, 2 * dim > potW(0);
    Dune::FieldVector < Scalar, 2 * dim > potNw(0);

    potW[0] = cellData1.potential(wPhaseIdx);
    potW[1] = cellData2.potential(wPhaseIdx);
    potW[2] = cellData4.potential(wPhaseIdx);

    potNw[0] = cellData1.potential(nPhaseIdx);
    potNw[1] = cellData2.potential(nPhaseIdx);
    potNw[2] = cellData4.potential(nPhaseIdx);

    //get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda1(cellData1.mobility(wPhaseIdx));
    lambda1[nPhaseIdx] = cellData1.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal1 = lambda1[wPhaseIdx] + lambda1[nPhaseIdx];

    //get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda2(cellData2.mobility(wPhaseIdx));
    lambda2[nPhaseIdx] = cellData2.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal2 = lambda2[wPhaseIdx] + lambda2[nPhaseIdx];

    //get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda4(cellData4.mobility(wPhaseIdx));
    lambda4[nPhaseIdx] = cellData4.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal4 = lambda4[wPhaseIdx] + lambda4[nPhaseIdx];

    std::vector < DimVector > lambda(4);

    lambda[0][0] = lambdaTotal1;
    lambda[0][1] = lambdaTotal1;
    lambda[1][0] = lambdaTotal2;
    lambda[1][1] = lambdaTotal2;
    lambda[3][0] = lambdaTotal4;
    lambda[3][1] = lambdaTotal4;

    Scalar potentialDiffW12 = 0;
    Scalar potentialDiffW14 = 0;
    Scalar potentialDiffW24 = 0;

    Scalar potentialDiffNw12 = 0;
    Scalar potentialDiffNw14 = 0;
    Scalar potentialDiffNw24 = 0;

    //flux vector
    Dune::FieldVector<Scalar, 3> fluxW(0);
    Dune::FieldVector<Scalar, 3> fluxNw(0);

    Dune::FieldMatrix<Scalar, dim, 2 * dim - dim + 1> T(0);
    DimVector Tu(0);
    Dune::FieldVector<Scalar, 2 * dim - dim + 1> u(0);

    int transmissibilityType = this->transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 0, 1, 3, 3);

    if (transmissibilityType == TransmissibilityCalculator::rightTriangle)
    {
        u[0] = potW[1];
        u[1] = potW[2];
        u[2] = potW[0];

        T.mv(u, Tu);

        fluxW[0] = Tu[1];
        potentialDiffW12 = Tu[1];

        u[0] = potNw[1];
        u[1] = potNw[2];
        u[2] = potNw[0];

        T.mv(u, Tu);

        fluxNw[0] = Tu[1];
        potentialDiffNw12 = Tu[1];
    }
    else if (transmissibilityType == TransmissibilityCalculator::leftTriangle)
    {
        u[0] = potW[0];
        u[1] = potW[2];
        u[2] = potW[1];

        T.mv(u, Tu);

        fluxW[0] = Tu[1];
        potentialDiffW12 = Tu[1];

        u[0] = potNw[0];
        u[1] = potNw[2];
        u[2] = potNw[1];

        T.mv(u, Tu);

        fluxNw[0] = Tu[1];
        potentialDiffNw12 = Tu[1];
    }

    transmissibilityType = this->transmissibilityCalculator_.calculateLeftHNTransmissibility(T, interactionVolume, lambda, 1, 3, 0);

    if (transmissibilityType == TransmissibilityCalculator::leftTriangle)
    {
        u[0] = potW[1];
        u[1] = potW[0];
        u[2] = potW[2];

        T.mv(u, Tu);

        fluxW[1] = Tu[1];
        potentialDiffW24 = Tu[1];

        u[0] = potNw[1];
        u[1] = potNw[0];
        u[2] = potNw[2];

        T.mv(u, Tu);

        fluxNw[1] = Tu[1];
        potentialDiffNw24 = Tu[1];
    }

    transmissibilityType = this->transmissibilityCalculator_.calculateRightHNTransmissibility(T, interactionVolume, lambda, 3, 0, 1);

    if (transmissibilityType == TransmissibilityCalculator::rightTriangle)
    {
        u[0] = potW[0];
        u[1] = potW[1];
        u[2] = potW[2];

        T.mv(u, Tu);

        fluxW[2] = Tu[1];
        potentialDiffW14 = -Tu[1];

        u[0] = potNw[0];
        u[1] = potNw[1];
        u[2] = potNw[2];

        T.mv(u, Tu);

        fluxNw[2] = Tu[1];
        potentialDiffNw14 = -Tu[1];
    }

    //store potentials for further calculations (saturation, ...) -> maybe add new potential to old one!!
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 0), potentialDiffW12);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 0), potentialDiffNw12);
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 1), potentialDiffW14);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 1), potentialDiffNw14);
    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 0), potentialDiffW24);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 0), potentialDiffNw24);
    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -potentialDiffW12);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -potentialDiffNw12);
    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 0), -potentialDiffW14);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 0), -potentialDiffNw14);
    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 1), -potentialDiffW24);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 1), -potentialDiffNw24);

    //compute mobilities of face 1
    Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
    lambda12Upw[wPhaseIdx] = (potentialDiffW12 >= 0) ? lambda1[wPhaseIdx] : lambda2[wPhaseIdx];
    lambda12Upw[nPhaseIdx] = (potentialDiffNw12 >= 0) ? lambda1[nPhaseIdx] : lambda2[nPhaseIdx];

    //compute mobilities of face 4
    Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
    lambda14Upw[wPhaseIdx] = (potentialDiffW14 >= 0) ? lambda1[wPhaseIdx] : lambda4[wPhaseIdx];
    lambda14Upw[nPhaseIdx] = (potentialDiffNw14 >= 0) ? lambda1[nPhaseIdx] : lambda4[nPhaseIdx];

    //compute mobilities of face 2
    Dune::FieldVector<Scalar, numPhases> lambda24Upw(0.0);
    lambda24Upw[wPhaseIdx] = (potentialDiffW24 >= 0) ? lambda2[wPhaseIdx] : lambda4[wPhaseIdx];
    lambda24Upw[nPhaseIdx] = (potentialDiffNw24 >= 0) ? lambda2[nPhaseIdx] : lambda4[nPhaseIdx];

    for (int i = 0; i < numPhases; i++)
    {
        // evaluate parts of velocity --> always take the normal for which the flux is calculated!
        DimVector vel12 = interactionVolume.getNormal(0, 0);
        DimVector vel14 = interactionVolume.getNormal(3, 0);
        DimVector vel24 = interactionVolume.getNormal(1, 0);
        DimVector vel21 = interactionVolume.getNormal(0, 0);
        DimVector vel41 = interactionVolume.getNormal(3, 0);
        DimVector vel42 = interactionVolume.getNormal(1, 0);

        Dune::FieldVector<Scalar, 3> flux(0);
        switch (i)
        {
        case wPhaseIdx:
        {
            flux = fluxW;
            break;
        }
        case nPhaseIdx:
        {
            flux = fluxNw;
            break;
        }
        }

        vel12 *= flux[0] / (2 * interactionVolume.getFaceArea(0, 0)); //divide by 2 because the flux is related to the half face!
        vel14 *= flux[2] / (2 * interactionVolume.getFaceArea(3, 0));
        vel24 *= flux[1] / (2 * interactionVolume.getFaceArea(1, 0));
        vel21 *= flux[0] / (2 * interactionVolume.getFaceArea(0, 0));
        vel41 *= flux[2] / (4 * interactionVolume.getFaceArea(3, 0));
        vel42 *= flux[1] / (4 * interactionVolume.getFaceArea(1, 0));

        Scalar lambdaT12 = lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx];
        Scalar lambdaT14 = lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx];
        Scalar lambdaT24 = lambda24Upw[wPhaseIdx] + lambda24Upw[nPhaseIdx];
        Scalar fracFlow12 = (lambdaT12 > ParentType::threshold_) ? lambda12Upw[i] / (lambdaT12) : 0.0;
        Scalar fracFlow14 = (lambdaT14 > ParentType::threshold_) ? lambda14Upw[i] / (lambdaT14) : 0.0;
        Scalar fracFlow24 = (lambdaT24 > ParentType::threshold_) ? lambda24Upw[i] / (lambdaT24) : 0.0;

        vel12 *= fracFlow12;
        vel14 *= fracFlow14;
        vel24 *= fracFlow24;
        vel21 *= fracFlow12;
        vel41 *= fracFlow14;
        vel42 *= fracFlow24;

        if (innerBoundaryVolumeFaces[globalIdx1][interactionVolume.getIndexOnElement(0, 0)])
        {
            vel12 *= 2;
            vel21 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx1][interactionVolume.getIndexOnElement(0, 1)])
        {
            vel14 *= 2;
            vel41 *= 2;
        }
        if (innerBoundaryVolumeFaces[globalIdx2][interactionVolume.getIndexOnElement(1, 0)])
        {
            vel24 *= 2;
            vel42 *= 2;
        }

        //store velocities
        //set velocity
        cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 0), vel12);
        cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 1), vel14);
        cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 0), vel24);
        cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 1), vel21);
        cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 0), vel41);
        cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 1), vel42);
    }
    //set velocity marker
    cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 0));
    cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 1));
    cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 0));
    cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 1));
    cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 0));
    cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 1));
}

} // end namespace Dumux
#endif
