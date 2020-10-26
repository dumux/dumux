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
 * \brief  Velocity field from a finite volume solution of a pressure equation using a grid adaptive MPFA L-method.
 */
#ifndef DUMUX_MPFAL2DPRESSUREVELOCITIES2P_ADAPTIVE_HH
#define DUMUX_MPFAL2DPRESSUREVELOCITIES2P_ADAPTIVE_HH

#include <dune/common/float_cmp.hh>

#include "2dpressureadaptive.hh"
#include "2dvelocityadaptive.hh"

#include <dumux/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Class for the calculation of velocities from the  pressure solution of an IMPES scheme using a grid adaptive MPFA L-method.
 *
 * Can be used for calculating the complete velocity field before the solution of the transport equation (more efficient),
 * or for face-wise velocity calculation directly in the transport solution (less efficient).
 *
 * Remark1: only for 2-D quadrilateral grids!
 *
 * Remark2: can use UGGrid, ALUGrid
 *
 * Remark3: Allowed difference in grid levels of two neighboring cells: 1
 *
 * \tparam TypeTag The problem Type Tag
 */
template<class TypeTag> class FvMpfaL2dPressureVelocity2pAdaptive: public FvMpfaL2dPressure2pAdaptive<TypeTag>
{
    using ParentType = FvMpfaL2dPressure2pAdaptive<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using Grid = typename GridView::Grid;
    using IndexSet = typename GridView::IndexSet;
    using Intersection = typename GridView::Intersection;

    using Geometry = typename Grid::template Codim<0>::Entity::Geometry;
    using JacobianTransposed = typename Geometry::JacobianTransposed ;

    using InteractionVolume = typename ParentType::InteractionVolume;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        sw = Indices::saturationW,
        sn = Indices::saturationNw
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressEqIdx = Indices::pressureEqIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };

    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using DimVector = Dune::FieldVector<Scalar, dim>;

public:
    /*!
     * \brief Constructs a FvMpfaL2dPressureVelocity2pAdaptive object
     *
     * \param problem A problem class object
     */
    FvMpfaL2dPressureVelocity2pAdaptive(Problem& problem) :
            ParentType(problem), problem_(problem), velocity_(problem)
    {
        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;

        calcVelocityInTransport_ = getParam<bool>("MPFA.CalcVelocityInTransport");
    }

    //! Calculates the velocities at all cell-cell interfaces.
    void calculateVelocity();

    //! Calculates the velocity at a cell-cell interface.
    void calculateVelocity(const Intersection&, CellData&);
    void calculateVelocityOnBoundary(const Intersection& intersection, CellData& cellData);

    //! Function for updating the velocity field if iterations are necessary in the transport solution
    void updateVelocity()
    {
        this->updateMaterialLaws();

        this->storePressureSolution();

        if (!calculateVelocityInTransport())
            calculateVelocity();
    }

    /*!
     * \brief Initializes pressure and velocity
     *
     * \copydetails FVPressure::initialize()
     */
    void initialize()
    {
        const auto element = *problem_.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(element));
        fluidState.setTemperature(problem_.temperature(element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);

        ParentType::initialize();
        velocity_.initialize();
        if (!calculateVelocityInTransport())
            calculateVelocity();

        return;
    }

    /*!
     * \brief Pressure and velocity update
     *
     * \copydetails FVPressure::update()
     */
    void update()
    {
        ParentType::update();
        if (!calculateVelocityInTransport())
            calculateVelocity();
    }

    /*!
     * \brief Indicates if velocity is reconstructed in the pressure step or in the transport step
     *
     * Returns true (In the standard finite volume discretization the velocity is calculated during the saturation transport.)
     */
    bool calculateVelocityInTransport()
    {
        return calcVelocityInTransport_;
    }

    /*!
     * \brief Adds velocity output to the output file
     *
     * Adds the phase velocities or a total velocity (depending on the formulation) to the output.
     * If the VtkOutputLevel is equal to zero (default) only primary variables are written,
     * if it is larger than zero also secondary variables are written.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);
        velocity_.addOutputVtkFields(writer);
    }

private:
     Problem& problem_;
     FvMpfaL2dVelocity2pAdaptive<TypeTag> velocity_;

     Scalar density_[numPhases];
     Scalar viscosity_[numPhases];
     bool calcVelocityInTransport_;

     //! gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
     static const int pressureType_ = getPropValue<TypeTag, Properties::PressureFormulation>();
     //! gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
     static const int saturationType_ = getPropValue<TypeTag, Properties::SaturationFormulation>();
};
// end of template

/*!
 * \brief Calculates the velocities at a cell-cell interfaces for the entire simulation grid.
 *
 * Calculates the velocities at a cell-cell interfaces from a given pressure field for the entire simulation grid.
 */
template<class TypeTag>
void FvMpfaL2dPressureVelocity2pAdaptive<TypeTag>::calculateVelocity()
{
    // run through all elements
    for (const auto& vertex : vertices(problem_.gridView()))
    {
        int vIdxGlobal = problem_.variables().index(vertex);

        InteractionVolume& interactionVolume = this->interactionVolumes_[vIdxGlobal];

        if (interactionVolume.isInnerVolume())
        {
            if (interactionVolume.getElementNumber() == 4)
            {
                auto element1 = interactionVolume.getSubVolumeElement(0);
                auto element2 = interactionVolume.getSubVolumeElement(1);
                auto element3 = interactionVolume.getSubVolumeElement(2);
                auto element4 = interactionVolume.getSubVolumeElement(3);

                // cell index
                int eIdxGlobal1 = problem_.variables().index(element1);
                int eIdxGlobal2 = problem_.variables().index(element2);
                int eIdxGlobal3 = problem_.variables().index(element3);
                int eIdxGlobal4 = problem_.variables().index(element4);

                //get the cell Data
                CellData& cellData1 = problem_.variables().cellData(eIdxGlobal1);
                CellData& cellData2 = problem_.variables().cellData(eIdxGlobal2);
                CellData& cellData3 = problem_.variables().cellData(eIdxGlobal3);
                CellData& cellData4 = problem_.variables().cellData(eIdxGlobal4);

                velocity_.calculateInnerInteractionVolumeVelocity(interactionVolume, cellData1, cellData2, cellData3,
                                                                  cellData4, this->innerBoundaryVolumeFaces_);

            }
            else if (interactionVolume.getElementNumber() == 3)
            {
                auto element1 = interactionVolume.getSubVolumeElement(0);
                  auto element2 = interactionVolume.getSubVolumeElement(1);
                  auto element4 = interactionVolume.getSubVolumeElement(3);

                  // cell index
                  int eIdxGlobal1 = problem_.variables().index(element1);
                  int eIdxGlobal2 = problem_.variables().index(element2);
                  int eIdxGlobal4 = problem_.variables().index(element4);

                  //get the cell Data
                  CellData& cellData1 = problem_.variables().cellData(eIdxGlobal1);
                  CellData& cellData2 = problem_.variables().cellData(eIdxGlobal2);
                  CellData& cellData4 = problem_.variables().cellData(eIdxGlobal4);

                  velocity_.calculateHangingNodeInteractionVolumeVelocity(interactionVolume, cellData1, cellData2,
                                                                          cellData4, this->innerBoundaryVolumeFaces_);

            }
            else
            {
                DUNE_THROW(Dune::NotImplemented, "Unknown interactionvolume type!");
            }
        }

        // at least one face on boundary!
        else
        {
            for (int elemIdx = 0; elemIdx < 2 * dim; elemIdx++)
            {
                bool isOutside = false;
                for (int fIdx = 0; fIdx < dim; fIdx++)
                {
                    int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, fIdx);
                    if (interactionVolume.isOutsideFace(intVolFaceIdx))
                    {
                        isOutside = true;
                        break;
                    }
                }
                if (isOutside)
                {
                    continue;
                }
                // cell index
                int eIdxGlobal = problem_.variables().index(interactionVolume.getSubVolumeElement(elemIdx));
                //get the cell Data
                CellData& cellData = problem_.variables().cellData(eIdxGlobal);

                velocity_.calculateBoundaryInteractionVolumeVelocity(interactionVolume, cellData, elemIdx);
            }

        } // end boundaries

    } // end vertex iterator
    return;
} // end method calcTotalVelocity

/*!
 * \brief Calculates the velocity at a cell-cell interface.
 *
 * Calculates the velocity at a cell-cell interface from a given pressure field.
 *
 * \param intersection Intersection of two grid cells
 * \param cellData Object containing all model relevant cell data
 */
template<class TypeTag>
void FvMpfaL2dPressureVelocity2pAdaptive<TypeTag>::calculateVelocity(const Intersection& intersection, CellData& cellData)
{
    int numVertices = intersection.geometry().corners();

    auto elementI = intersection.inside();
    auto elementJ = intersection.outside();

    int levelI = elementI.level();
    int levelJ = elementJ.level();

    int eIdxGlobalI = problem_.variables().index(elementI);
    int eIdxGlobalJ = problem_.variables().index(elementJ);

    CellData& cellDataJ = problem_.variables().cellData(eIdxGlobalJ);

    const auto refElement = referenceElement(elementI);

    int indexInInside = intersection.indexInInside();
    int indexInOutside = intersection.indexInOutside();

    int fIdx = indexInInside;

    if (levelI < levelJ)
        fIdx = indexInOutside;

    std::vector<CellData> cellDataTemp(0);

    if (levelI != levelJ)
    {
        DimVector vel(0);
        cellData.fluxData().setVelocity(wPhaseIdx, indexInInside, vel);
        cellData.fluxData().setVelocity(nPhaseIdx, indexInInside, vel);
        cellData.fluxData().setUpwindPotential(wPhaseIdx, indexInInside, 0);
        cellData.fluxData().setUpwindPotential(nPhaseIdx, indexInInside, 0);

        cellDataJ.fluxData().setVelocity(wPhaseIdx, indexInOutside, vel);
        cellDataJ.fluxData().setVelocity(nPhaseIdx, indexInOutside, vel);
        cellDataJ.fluxData().setUpwindPotential(wPhaseIdx, indexInOutside, 0);
        cellDataJ.fluxData().setUpwindPotential(nPhaseIdx, indexInOutside, 0);
    }

    for (int vIdx = 0; vIdx < numVertices; vIdx++)
    {
        int localVertIdx = refElement.subEntity(fIdx, dim - 1, vIdx, dim);

        int vIdxGlobal = 0;
                if (levelI >= levelJ)
                {
                    vIdxGlobal = problem_.variables().index(elementI.template subEntity<dim>(localVertIdx));
                }
                else
                {
                    vIdxGlobal = problem_.variables().index(elementJ.template subEntity<dim>(localVertIdx));
                }

        InteractionVolume& interactionVolume = this->interactionVolumes_[vIdxGlobal];

        if (interactionVolume.isInnerVolume())
        {
            // cell index vector
            std::vector<int> eIdxGlobal(0);

            if (interactionVolume.getElementNumber() == 4)
            {
                auto element1 = interactionVolume.getSubVolumeElement(0);
                auto element2 = interactionVolume.getSubVolumeElement(1);
                auto element3 = interactionVolume.getSubVolumeElement(2);
                auto element4 = interactionVolume.getSubVolumeElement(3);

                eIdxGlobal.resize(4);

                eIdxGlobal[0] = problem_.variables().index(element1);
                eIdxGlobal[1] = problem_.variables().index(element2);
                eIdxGlobal[2] = problem_.variables().index(element3);
                eIdxGlobal[3] = problem_.variables().index(element4);

                //cell Data vector
                cellDataTemp.resize(4);

                cellDataTemp[0] = problem_.variables().cellData(eIdxGlobal[0]);
                cellDataTemp[1] = problem_.variables().cellData(eIdxGlobal[1]);
                cellDataTemp[2] = problem_.variables().cellData(eIdxGlobal[2]);
                cellDataTemp[3] = problem_.variables().cellData(eIdxGlobal[3]);

                velocity_.calculateInnerInteractionVolumeVelocity(interactionVolume, cellDataTemp[0], cellDataTemp[1],
                                                                  cellDataTemp[2], cellDataTemp[3], this->innerBoundaryVolumeFaces_);
            }
            else if (interactionVolume.getElementNumber() == 3)
            {
                auto element1 = interactionVolume.getSubVolumeElement(0);
                auto element2 = interactionVolume.getSubVolumeElement(1);
                auto element4 = interactionVolume.getSubVolumeElement(3);

                eIdxGlobal.resize(3);

                eIdxGlobal[0] = problem_.variables().index(element1);
                eIdxGlobal[1] = problem_.variables().index(element2);
                eIdxGlobal[2] = problem_.variables().index(element4);

                //cell Data vector
                cellDataTemp.resize(3);

                cellDataTemp[0] = problem_.variables().cellData(eIdxGlobal[0]);
                cellDataTemp[1] = problem_.variables().cellData(eIdxGlobal[1]);
                cellDataTemp[2] = problem_.variables().cellData(eIdxGlobal[2]);


                velocity_.calculateHangingNodeInteractionVolumeVelocity(interactionVolume, cellDataTemp[0], cellDataTemp[1],
                                                                        cellDataTemp[2], this->innerBoundaryVolumeFaces_);
            }
            else
            {
                DUNE_THROW(Dune::NotImplemented, "Unknown interactionvolume type!");
            }

            int size = cellDataTemp.size();
        for (int i = 0; i < size; i++)
        {
            if (eIdxGlobal[i] == eIdxGlobalI)
            {
                if (levelI >= levelJ)
                {
                 cellData.fluxData().setVelocity(wPhaseIdx, indexInInside,
                                                 cellDataTemp[i].fluxData().velocity(wPhaseIdx, indexInInside));
                 cellData.fluxData().setVelocity(nPhaseIdx, indexInInside,
                                                 cellDataTemp[i].fluxData().velocity(nPhaseIdx, indexInInside));
                 cellData.fluxData().setUpwindPotential(wPhaseIdx, indexInInside,
                                                        cellDataTemp[i].fluxData().upwindPotential(wPhaseIdx, indexInInside));
                 cellData.fluxData().setUpwindPotential(nPhaseIdx, indexInInside,
                                                        cellDataTemp[i].fluxData().upwindPotential(nPhaseIdx, indexInInside));

                 if (levelI > levelJ)
                 {
                     cellDataJ.fluxData().setVelocity(wPhaseIdx, indexInOutside,
                                                      cellDataTemp[i].fluxData().velocity(wPhaseIdx, indexInInside));
                     cellDataJ.fluxData().setVelocity(nPhaseIdx, indexInOutside,
                                                      cellDataTemp[i].fluxData().velocity(nPhaseIdx, indexInInside));
                     cellDataJ.fluxData().setUpwindPotential(wPhaseIdx, indexInOutside,
                                                             cellDataTemp[i].fluxData().upwindPotential(wPhaseIdx, indexInInside));
                     cellDataJ.fluxData().setUpwindPotential(nPhaseIdx, indexInOutside,
                                                             cellDataTemp[i].fluxData().upwindPotential(nPhaseIdx, indexInInside));

                 }
                }

            }
            else if (eIdxGlobal[i] == eIdxGlobalJ)
            {
                if (levelJ >= levelI)
                {
                cellDataJ.fluxData().setVelocity(wPhaseIdx, indexInOutside,
                                                 cellDataTemp[i].fluxData().velocity(wPhaseIdx, indexInOutside));
                cellDataJ.fluxData().setVelocity(nPhaseIdx, indexInOutside,
                                                 cellDataTemp[i].fluxData().velocity(nPhaseIdx, indexInOutside));
                cellDataJ.fluxData().setUpwindPotential(wPhaseIdx, indexInOutside,
                                                        cellDataTemp[i].fluxData().upwindPotential(wPhaseIdx, indexInOutside));
                cellDataJ.fluxData().setUpwindPotential(nPhaseIdx, indexInOutside,
                                                        cellDataTemp[i].fluxData().upwindPotential(nPhaseIdx, indexInOutside));

                if (levelJ > levelI)
                {
                    cellData.fluxData().setVelocity(wPhaseIdx, indexInInside,
                                                    cellDataTemp[i].fluxData().velocity(wPhaseIdx, indexInOutside));
                    cellData.fluxData().setVelocity(nPhaseIdx, indexInInside,
                                                    cellDataTemp[i].fluxData().velocity(nPhaseIdx, indexInOutside));
                    cellData.fluxData().setUpwindPotential(wPhaseIdx, indexInInside,
                                                           cellDataTemp[i].fluxData().upwindPotential(wPhaseIdx, indexInOutside));
                    cellData.fluxData().setUpwindPotential(nPhaseIdx, indexInInside,
                                                           cellDataTemp[i].fluxData().upwindPotential(nPhaseIdx, indexInOutside));

                }
                }
            }
        }
    }
    }
    if (levelI == levelJ)
    {
        cellData.fluxData().setVelocityMarker(indexInInside);
        cellDataJ.fluxData().setVelocityMarker(indexInOutside);
    }
}

/*!
 * \brief Calculates the velocity at a boundary.
 *
 * Calculates the velocity at a boundary from a given pressure field.
 *
 * \param intersection Boundary intersection
 * \param cellData Object containing all model relevant cell data
 */
template<class TypeTag>
void FvMpfaL2dPressureVelocity2pAdaptive<TypeTag>::calculateVelocityOnBoundary(const Intersection& intersection, CellData& cellData)
{
    auto element = intersection.inside();

    //get face index
    int isIndex = intersection.indexInInside();

    //get face normal
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

    BoundaryTypes bcType;
    //get boundary type
    problem_.boundaryTypes(bcType, intersection);
    PrimaryVariables boundValues(0.0);

    if (bcType.isDirichlet(pressEqIdx))
    {
        problem_.dirichlet(boundValues, intersection);

        // get global coordinates of cell centers
        const GlobalPosition& globalPosI = element.geometry().center();

        // center of face in global coordinates
        const GlobalPosition& globalPosJ = intersection.geometry().center();

        // get mobilities and fractional flow factors
        Scalar lambdaWI = cellData.mobility(wPhaseIdx);
        Scalar lambdaNwI = cellData.mobility(nPhaseIdx);

        // get capillary pressure
        Scalar pcI = cellData.capillaryPressure();

        // distance vector between barycenters
        GlobalPosition distVec = globalPosJ - globalPosI;

        // compute distance between cell centers
        Scalar dist = distVec.two_norm();

        //permeability vector at boundary
        // compute vectorized permeabilities
        DimMatrix meanPermeability(0);

        problem_.spatialParams().meanK(meanPermeability, problem_.spatialParams().intrinsicPermeability(element));

        Dune::FieldVector<Scalar, dim> permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);

        //determine saturation at the boundary -> if no saturation is known directly at the boundary use the cell saturation
        Scalar satW = 0;
        if (bcType.isDirichlet(satEqIdx))
        {
            switch (saturationType_)
            {
            case sw:
            {
                satW = boundValues[saturationIdx];
                break;
            }
            case sn:
            {
                satW = 1 - boundValues[saturationIdx];
                break;
            }
            }
        }
        else
        {
            satW = cellData.saturation(wPhaseIdx);
        }

        const Scalar pressBound = boundValues[pressureIdx];

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

        const Scalar pcBound = fluidMatrixInteraction.pc(satW);

        //determine phase pressures from primary pressure variable
        Scalar pressWBound = 0;
        Scalar pressNwBound = 0;
        if (pressureType_ == pw)
        {
            pressWBound = pressBound;
            pressNwBound = pressBound + pcBound;
        }
        else if (pressureType_ == pn)
        {
            pressWBound = pressBound - pcBound;
            pressNwBound = pressBound;
        }

        const Scalar lambdaWBound = fluidMatrixInteraction.krw(satW) / viscosity_[wPhaseIdx];
        const Scalar lambdaNwBound = fluidMatrixInteraction.krn(satW) / viscosity_[nPhaseIdx];

        Scalar potentialDiffW = cellData.fluxData().upwindPotential(wPhaseIdx, isIndex);
        Scalar potentialDiffNw = cellData.fluxData().upwindPotential(nPhaseIdx, isIndex);

        //calculate potential gradient
        potentialDiffW = (cellData.pressure(wPhaseIdx) - pressWBound);
        potentialDiffNw = (cellData.pressure(nPhaseIdx) - pressNwBound);

        potentialDiffW += density_[wPhaseIdx] * (distVec * problem_.gravity());
        potentialDiffNw += density_[nPhaseIdx] * (distVec * problem_.gravity());

        //store potential gradients for further calculations
        cellData.fluxData().setUpwindPotential(wPhaseIdx, isIndex, potentialDiffW);
        cellData.fluxData().setUpwindPotential(nPhaseIdx, isIndex, potentialDiffNw);

        //do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaW = (potentialDiffW > 0.) ? lambdaWI : lambdaWBound;
        lambdaW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaW;
        Scalar lambdaNw = (potentialDiffNw > 0.) ? lambdaNwI : lambdaNwBound;
        lambdaNw = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (lambdaNwI + lambdaNwBound) : lambdaNw;


        Scalar scalarPerm = permeability.two_norm();

        //calculate the gravity term
        Dune::FieldVector<Scalar, dimWorld> velocityW(unitOuterNormal);
        Dune::FieldVector<Scalar, dimWorld> velocityNw(unitOuterNormal);

        //calculate unit distVec
        distVec /= dist;
        Scalar areaScaling = (unitOuterNormal * distVec);
        //this treatment of g allows to account for gravity flux through faces where the face normal
        //has no z component (e.g. parallelepiped grids)
        Scalar gravityTermW = (problem_.gravity() * distVec) * density_[wPhaseIdx] * areaScaling;
        Scalar gravityTermNw = (problem_.gravity() * distVec) * density_[nPhaseIdx] * areaScaling;

        //calculate velocity depending on the pressure used -> use pc = pn - pw
        switch (pressureType_)
        {
        case pw:
        {
            velocityW *= lambdaW * scalarPerm * ((cellData.pressure(wPhaseIdx) - pressBound) / dist + gravityTermW);
            velocityNw *= lambdaNw * scalarPerm * ((cellData.pressure(wPhaseIdx) - pressBound) / dist + gravityTermNw)
                    + 0.5 * (lambdaNwI + lambdaNwBound) * scalarPerm * (pcI - pcBound) / dist;
            break;
        }
        case pn:
        {
            velocityW *= lambdaW * scalarPerm * ((cellData.pressure(nPhaseIdx) - pressBound) / dist + gravityTermW)
                    - 0.5 * (lambdaWI + lambdaWBound) * scalarPerm * (pcI - pcBound) / dist;
            velocityNw *= lambdaNw * scalarPerm * ((cellData.pressure(nPhaseIdx) - pressBound) / dist + gravityTermNw);
            break;
        }
        }

        //store velocities
        cellData.fluxData().setVelocity(wPhaseIdx, isIndex, velocityW);
        cellData.fluxData().setVelocity(nPhaseIdx, isIndex, velocityNw);
        cellData.fluxData().setVelocityMarker(isIndex);

    } //end dirichlet boundary

    else if (bcType.isNeumann(pressEqIdx))
    {
        problem_.neumann(boundValues, intersection);

        Dune::FieldVector<Scalar, dimWorld> velocityW(unitOuterNormal);
        Dune::FieldVector<Scalar, dimWorld> velocityNw(unitOuterNormal);

        velocityW *= boundValues[wPhaseIdx];
        velocityNw *= boundValues[nPhaseIdx];

            velocityW /= density_[wPhaseIdx];
            velocityNw /= density_[nPhaseIdx];

        //store potential gradients for further calculations
        cellData.fluxData().setUpwindPotential(wPhaseIdx, isIndex, boundValues[wPhaseIdx]);
        cellData.fluxData().setUpwindPotential(nPhaseIdx, isIndex, boundValues[nPhaseIdx]);

        cellData.fluxData().setVelocity(wPhaseIdx, isIndex, velocityW);
        cellData.fluxData().setVelocity(nPhaseIdx, isIndex, velocityNw);
        cellData.fluxData().setVelocityMarker(isIndex);
    } //end neumann boundary
    else
    {
        DUNE_THROW(Dune::NotImplemented, "No valid boundary condition type defined for pressure equation!");
    }
}

}
// end of Dune namespace
#endif
