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
 * \brief Base class for all models which use the ZeroEq box model.
 */
#ifndef DUMUX_ZEROEQ_MODEL_HH
#define DUMUX_ZEROEQ_MODEL_HH

#include "zeroeqindices.hh"
#include "zeroeqfluxvariables.hh"
#include "zeroeqproblem.hh"
#include "zeroeqproperties.hh"
#include <dumux/freeflow/stokes/stokesmodel.hh>

namespace Dumux
{
/*!
 * \ingroup BoxZeroEqModel
 * \brief Adaptation of the box scheme to the ZeroEq model.
 *
 * This model implements an single-phase isothermal free flow
 * solving the mass and the momentum balance. For the momentum balance
 * the Reynolds-averaged Navier-Stokes (RANS) equation with zero equation
 * (algebraic) turbulence model is used.
 *
 * Mass balance:
 * \f[
 *  \frac{\partial \varrho_\textrm{g}}{\partial t}
 *  + \boldsymbol{\nabla}\boldsymbol{\cdot}\left(\varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g}\right)
 *  - q_\textrm{g} = 0
 * \f]
 *
 * Momentum Balance:
 * \f[
 *   \frac{\partial \left(\varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g}\right)}{\partial t}
 *   + \boldsymbol{\nabla} \boldsymbol{\cdot} \left(
 *     \varrho_\textrm{g} {\boldsymbol{v}_\textrm{g} {\boldsymbol{v}}_\textrm{g}}
 *     - \left[ \mu_\textrm{g} + \mu_\textrm{g,t} \right]
 *       \left(\boldsymbol{\nabla} \boldsymbol{v}_\textrm{g}
 *             + \boldsymbol{\nabla} \boldsymbol{v}_\textrm{g}^T \right)
 *   \right)
 *   + \left(p_\textrm{g} {\bf {I}} \right)
 *   - \varrho_\textrm{g} {\bf g} = 0,
 * \f]
 *
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
 */
template<class TypeTag>
class ZeroEqModel : public GET_PROP_TYPE(TypeTag, BaseStokesModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        intervals = GET_PROP_VALUE(TypeTag, NumberOfIntervals),
        bboxMinIsWall = GET_PROP_VALUE(TypeTag, BBoxMinIsWall),
        bboxMaxIsWall = GET_PROP_VALUE(TypeTag, BBoxMaxIsWall),
        walls = (bboxMinIsWall ? 1 : 0) + (bboxMaxIsWall ? 1 : 0),
        prec = Indices::scvDataPrecision, // precision of scv data
        width = Indices::scvDataWidth // width of column
    };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;


public:
    ZeroEqModel()
        : flowNormal_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, FlowNormal))
        , wallNormal_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, WallNormal))
        , writeAllSCVData_(GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, WriteAllSCVData))
    {
        eps_ = 1e-6;

        // check whether sand grain roughness may be used
        if ((GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMinSandGrainRoughness) > 0
              || GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMaxSandGrainRoughness) > 0)
            && surfaceRoughnessNotImplemented())
        {
            Dune::dwarn << "warning: surface roughness is not implemented for eddy viscosity model "
                        << GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyViscosityModel)
                        << "." << std::endl;
        }
    }

    /*!
     * \brief Calculate the fluxes across a certain layer in the domain.
     *        The layer is situated perpendicular to the coordinate axis "coord" and cuts
     *        the axis at the value "coordVal".
     *
     * \param globalSol The global solution vector.
     * \param flux A vector to store the flux.
     * \param axis The dimension, perpendicular to which the layer is situated.
     * \param coordVal The (Scalar) coordinate on the axis, at which the layer is situated.
     */
    void calculateFluxAcrossLayer(const SolutionVector &globalSol, Dune::FieldVector<Scalar, numEq> &flux, int axis, Scalar coordVal)
    {
        GlobalPosition globalI, globalJ;
        PrimaryVariables tmpFlux(0.0);

        FVElementGeometry fvGeometry;
        ElementVolumeVariables elemVolVars;

        // Loop over elements
        ElementIterator eIt = this->problem_.gridView().template begin<0>();
        ElementIterator eEndIt = this->problem_.gridView().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            if (eIt->partitionType() != Dune::InteriorEntity)
                continue;

            fvGeometry.update(this->gridView_(), *eIt);
            elemVolVars.update(this->problem_(), *eIt, fvGeometry);
            this->localResidual().evalFluxes(*eIt, elemVolVars);

            bool hasLeft = false;
            bool hasRight = false;
            for (int i = 0; i < fvGeometry.numVertices; i++) {
                const GlobalPosition &globalPos = fvGeometry.subContVol[i].global;
                if (globalI[axis] < coordVal)
                    hasLeft = true;
                else if (globalI[axis] >= coordVal)
                    hasRight = true;
            }
            if (!hasLeft || !hasRight)
                continue;

            for (int i = 0; i < fvGeometry.numVertices; i++) {
                const GlobalPosition &globalPos = fvGeometry.subContVol[i].global;
                if (globalI[axis] < coordVal)
                    flux += this->localResidual().residual(i);
            }
        }

        flux = this->problem_.gridView().comm().sum(flux);
    }

    /*!
     * \brief Write vtk and additional plain text scv-data.
     *
     * The routine for the vtk data should be the same as in dumux-stable.
     *
     * The scv-data files contain mainly information for the turbulence model.
     *
     * \param sol The solution vector.
     * \param writer The writer for multi-file VTK datasets.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VelocityField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        ScalarField &pN = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;

        ElementBoundaryTypes elemBcTypes;

        VertexIterator vIt = this->gridView_().template begin<dim>();
        VertexIterator vEndIt = this->gridView_().template end<dim>();
        for (; vIt != vEndIt; ++vIt)
        {}


        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();

        for (; eIt != eEndIt; ++eIt)
        {
            int idx = this->elementMapper().map(*eIt);
            rank[idx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt, fvGeometry);


            int numLocalVerts = eIt->template count<dim>();
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int vIdxGlobal = this->vertexMapper().map(*eIt, i, dim);
                volVars.update(sol[vIdxGlobal],
                               this->problem_(),
                               *eIt,
                               fvGeometry,
                               i,
                               false);

                pN[vIdxGlobal] = volVars.pressure();
                delP[vIdxGlobal] = volVars.pressure() - 1e5;
                rho[vIdxGlobal] = volVars.density();
                mu[vIdxGlobal] = volVars.dynamicViscosity();
                velocity[vIdxGlobal] = volVars.velocity();
            };
        }

        writer.attachVertexData(pN, "P");
        writer.attachVertexData(delP, "delP");
        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachVertexData(velocity, "v", dim);


        // ensure that the actual values are given out
        asImp_().updateWallProperties();

        eIt = this->gridView_().template begin<0>();
        eEndIt = this->gridView_().template end<0>();

        for (; eIt != eEndIt; ++eIt)
        {
            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt, fvGeometry);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false);

            IntersectionIterator isIt = this->gridView_().ibegin(*eIt);
            const IntersectionIterator &endIt = this->gridView_().iend(*eIt);

            for (; isIt != endIt; ++isIt)
            {
                int fIdx = isIt->indexInInside();

                FluxVariables fluxVars(this->problem_(),
                                                    *eIt,
                                                    fvGeometry,
                                                    fIdx,
                                                    elemVolVars,
                                                    false);

                GlobalPosition globalPos = fvGeometry.subContVolFace[fIdx].ipGlobal;

                if (asImp_().shouldWriteSCVData(globalPos))
                    asImp_().writeSCVData(volVars, fluxVars, globalPos);
            }
        }
    }


    /**********************************
     * SCVDATA
     * -------
     * This section contains functions to zeroeq output related functions.
     */

    /*!
     * \brief Returns true if a scv-data file for gnuplot evaluation
     *        should be given out.
     *
     * \param globalPos Global Position.
     */
    const bool shouldWriteSCVData(const GlobalPosition &globalPos) const
    {
        if (writeAllSCVData_ > 1.0)
            return true;
        else if (writeAllSCVData_ < 0.0)
            return false;
        else
        {
            Scalar relGlobalPosMin = (globalPos[flowNormal_] - 1e-3 - this->problem_().bBoxMin()[flowNormal_])
                                     / (this->problem_().bBoxMax()[flowNormal_] - this->problem_().bBoxMin()[flowNormal_]);
            Scalar relGlobalPosMax = (globalPos[flowNormal_] + 1e-3 - this->problem_().bBoxMin()[flowNormal_])
                                     / (this->problem_().bBoxMax()[flowNormal_] - this->problem_().bBoxMin()[flowNormal_]);
            return (relGlobalPosMax > writeAllSCVData_ && relGlobalPosMin < writeAllSCVData_);
        }
    }

    /*!
     * \brief Writes scv-data file for gnuplot evaluation.
     *
     * \param volVars Volume Variables of current scv.
     * \param fluxVars Flux Variables of current element.
     * \param globalPos Global Position.
     */
    void writeSCVData(const VolumeVariables &volVars, const FluxVariables &fluxVars, const GlobalPosition &globalPos)
    {
        std::ofstream fileAct ("");
        std::stringstream stream ("");

        char flowNormalChar = 'x';
        if (flowNormal_ == 1) flowNormalChar = 'y';
        stream << "scvData" << "-" << flowNormalChar << std::setprecision(4)
               << (globalPos[flowNormal_] - this->problem_().bBoxMin()[flowNormal_])
                  / (this->problem_().bBoxMax()[flowNormal_] - this->problem_().bBoxMin()[flowNormal_])
               << ".dat";
        std::string fileString2(stream.str());
        stream.str(""); stream.clear();
        const char* fileNameAct = fileString2.c_str();

        int posIdx = getPosIdx(globalPos);
        int wallIdx = getWallIdx(globalPos, posIdx);
        if (!wall[wallIdx].headerWritten[posIdx])
        {
            fileAct.open(fileNameAct, std::ios::out);

            asImp_().writeSCVHeader(stream, fluxVars, globalPos);
            stream << std::endl;
            if (bboxMinIsWall)
            {
                asImp_().writeWallHeader(stream, posIdx, 0);
                stream << std::endl;
            }
            if (bboxMaxIsWall)
            {
                asImp_().writeWallHeader(stream, posIdx, 1);
                stream << std::endl;
            }
            asImp_().writeDataHeader(stream, posIdx);
            stream << std::endl;

            fileAct << stream.str();
            stream.str(""); stream.clear();
            fileAct.close();
            for (int wIdx = 0; wIdx < walls; ++wIdx)
                wall[wIdx].headerWritten[posIdx] = true;
        }
        fileAct.open(fileNameAct, std::ios::app);

        asImp_().writeSCVDataValues(stream, volVars, fluxVars, globalPos);
        stream << std::endl;
        fileAct << stream.str();
        fileAct.close();
    }

    /*!
     * \brief Writes scv-data header, containing basic properties, which are constant
     *        at this point.
     *
     * \param stream Output Filestream.
     * \param fluxVars Flux Variables of current element.
     * \param globalPos Global Position.
     */
    void writeSCVHeader(std::stringstream &stream, const FluxVariables &fluxVars, const GlobalPosition &globalPos)
    {
        int posIdx = getPosIdx(globalPos);
        stream << std::setprecision(prec) << "# "
               << std::setw(width-2) << "over x[m]:"
               << std::setw(width) << globalPos[flowNormal_]
               << std::setw(width) << "at t[s]:"
               << std::setw(width) << this->problem_().timeManager().time()
               << std::setw(width) << "Re_d[-]:"
               << std::setw(width) << wall[0].maxVelocityAbs[posIdx][flowNormal_]
                                      * (this->problem_().bBoxMax()[wallNormal_] - this->problem_().bBoxMin()[wallNormal_])
                                      / fluxVars.kinematicViscosity()
               << std::setw(width) << "Re_x[-]:"
               << std::setw(width) << wall[0].maxVelocityAbs[posIdx][flowNormal_]
                                      * (this->problem_().bBoxMax()[flowNormal_] - this->problem_().bBoxMin()[flowNormal_])
                                      / fluxVars.kinematicViscosity()
               << std::endl;
        stream << std::setprecision(prec) << "# "
               << std::setw(width-2)<< "eddyVMod"
               << std::setw(width) << GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyViscosityModel)
               << " - " << std::setw(width*2-3) << std::left  << eddyViscosityModelName() << std::right;
    }

    /*!
     * \brief Writes property names to wall-data header.
     *
     * \param stream Output Filestream.
     * \param posIdx Position Index of current Global Position.
     * \param wallIdx Wall Index of current Global Position.
     */
    void writeWallHeader(std::stringstream &stream, int posIdx, int wallIdx)
    {
        stream << std::setprecision(prec) << "# "
               << std::setw(width-4) << "wall[" << wallIdx << "]"
               << std::setw(width) << "u_*[m/s]:"
               << std::setw(width) << std::sqrt(wall[wallIdx].wallShearStress[posIdx] / wall[wallIdx].wallDensity[posIdx])
               << std::setw(width) << "vs_5[m]:"
               << std::setw(width) << 5.0 * wall[wallIdx].wallKinematicViscosity[posIdx]
                                      / std::sqrt(wall[wallIdx].wallShearStress[posIdx] / wall[wallIdx].wallDensity[posIdx])
               << std::setw(width) << "bl_tot[m]:"
               << std::setw(width) << std::abs(wall[wallIdx].boundaryLayerThickness[posIdx])
               << std::setw(width) << "ks_+[-]:"
               << std::setw(width) << wall[wallIdx].sandGrainRoughness[posIdx]
                                      * std::sqrt(std::abs(wall[wallIdx].wallVelGrad[posIdx]) * wall[wallIdx].wallKinematicViscosity[posIdx])
                                      / wall[wallIdx].wallKinematicViscosity[posIdx]
               << std::setw(width) << "velGrad[1/s]:"
               << std::setw(width) << wall[wallIdx].wallVelGrad[posIdx]
               << std::setw(width) << "nu[m^2/s]:"
               << std::setw(width) << wall[wallIdx].wallKinematicViscosity[posIdx]
               << std::setw(width) << "rho[kg/m^3]:"
               << std::setw(width) << wall[wallIdx].wallDensity[posIdx]
               << std::setw(width) << "tau[Pa]:"
               << std::setw(width) << wall[wallIdx].wallShearStress[posIdx];
    }

    /*!
     * \brief Writes property names to scv-data header.
     *
     * \param stream Output Filestream.
     * \param posIdx Position Index of current Global Position.
     */
    void writeDataHeader(std::stringstream &stream, int posIdx)
    {
        stream << std::setw(1) << "#" << std::setw(width-1);
        for (int numIdx = 1; numIdx < 30; ++numIdx)
            stream << numIdx << std::setw(width);
        stream << std::endl;
        stream << std::setw(1) << "#" << std::setw(width-1) << "x [m]"
               << std::setw(width) << "y [m]"
               << std::setw(width) << "u [m/s]"
               << std::setw(width) << "v [m/s]"
               << std::setw(width) << "u/u_max [-]"
               << std::setw(width) << "v/v_max [-]"
               << std::setw(width) << "y^+ [-]"
               << std::setw(width) << "u^+ [-]"
               << std::setw(width) << "velGrad [1/s]"
               << std::setw(width) << "rho [kg/m^3]"
               << std::setw(width) << "nu [m^2/s]"
               << std::setw(width) << "nu_t [m^2/s]";
    }

    /*!
     * \brief Writes values to scv-data.
     *
     * \param stream Output Filestream.
     * \param volVars Volume Variables of current scv.
     * \param fluxVars Flux Variables of current element.
     * \param globalPos Global Position.
     */
    void writeSCVDataValues(std::stringstream &stream, const VolumeVariables &volVars, const FluxVariables &fluxVars, const GlobalPosition &globalPos)
    {
        int posIdx = getPosIdx(globalPos);
        int wallIdx = getWallIdx(globalPos, posIdx);
        stream << std::setw(width) << std::setprecision(prec+1) << globalPos[0]
               << std::setw(width) << globalPos[1]
               << std::setw(width) << std::setprecision(prec) << fluxVars.velocity()[0]
               << std::setw(width) << fluxVars.velocity()[1]
               << std::setw(width) << fluxVars.velocity()[0] / wall[wallIdx].maxVelocityAbs[posIdx][0]
               << std::setw(width) << fluxVars.velocity()[1] / wall[wallIdx].maxVelocityAbs[posIdx][1]
               << std::setw(width) << fluxVars.frictionVelocityWall() * std::abs(wall[wallIdx].wallPos[posIdx] - globalPos[wallNormal_]) / fluxVars.kinematicViscosity()
               << std::setw(width) << fluxVars.velocity()[flowNormal_] / fluxVars.frictionVelocityWall()
               << std::setw(width) << fluxVars.velocityGrad()[flowNormal_][wallNormal_]
               << std::setw(width) << fluxVars.density()
               << std::setw(width) << fluxVars.dynamicViscosity() / fluxVars.density()
               << std::setw(width) << fluxVars.dynamicEddyViscosity() / fluxVars.density();
    }



    /**********************************
     * WALL PROPERTIES
     * ---------------
     * This section contains functions which treat the wall handling.
     */

     /*!
     * \brief Container for all necessary information, needed to calculate the
     *        eddy viscosity at a certain point in relation to the wall distance
     *        and fluid/flow properties at the wall.
     */
    struct WallProperties
    {
    public:
        bool isBBoxMinWall;                    //!< Actual wall properties are located on bboxmin or bboxmax.
        Scalar wallPos[intervals];             //!< Position of the wall interval in global coordinates.
        Scalar sandGrainRoughness[intervals];  //!< Sand grain roughness.
        Scalar boundaryLayerThickness[intervals];             //!< Domain influenced by this wall.
        Scalar boundaryLayerThicknessCalculated[intervals];   //!< Boundary layer thickness based on v = v_99.
        Scalar viscousSublayerThicknessCalculated[intervals]; //!< Viscous sublayer thickness based on y^+ <= 5.
        Scalar crossLength[intervals];         //!< Switching point for the Baldwin-Lomax model.
        Scalar maxVelocity[intervals][dim];    //!< Max velocity related to this wall.
        Scalar maxVelocityAbs[intervals][dim]; //!< Max velocity at this interval position.
        Scalar minVelocity[intervals][dim];    //!< Min velocity related to this wall.
        Scalar wallDensity[intervals];         //!< Fluid density at the wall.
        Scalar wallKinematicViscosity[intervals];             //!< Kinematic viscosity at the wall.
        Scalar wallVelGrad[intervals];         //!< Velocity gradient at the wall, in wall normal direction.
        Scalar wallShearStress[intervals];     //!< Shear stress at the wall.
        Scalar fMax[intervals];                //!< Max value of the f function for Baldwin-Lomax model.
        Scalar yMax[intervals];                //!< Distance of position where fMax occurs.
        Scalar maxMassFraction[intervals];     //!< Max mass fraction related to this wall.
        Scalar maxMoleFraction[intervals];     //!< Max mole fraction related to this wall.
        Scalar maxTemperature[intervals];      //!< Max temperature related to this wall.
        int isInterpolated[intervals];         //!< Value is only interpolated between neighboring wall intervals.
        int fluxValuesCount[intervals];        //!< Number of flux values contributing to the interval.
        int wallValuesCount[intervals];        //!< Number of values contributing properties directly at the wall.
        bool headerWritten[intervals];         //!< Header of scv-data file was already written.
        WallProperties() {}                    //!< Constructor for wall properties.
    };
    WallProperties wall[walls];

    /*!
     * \brief Initializes the wall structure with values.
     *
     * This should only be done before updating (preTimeStepFunction and
     * only once per time Step), because the crossLength
     * for the Baldwin Lomax model is reset.
     */
    void resetWallProperties()
    {
        for (int wallIdx = 0; wallIdx < walls; ++wallIdx)
        {
            for (int posIdx = 0; posIdx < intervals; ++posIdx)
            {
                if (walls == 1)
                {
                    if (bboxMinIsWall)
                    {
                        wall[wallIdx].wallPos[posIdx] = this->problem_().bBoxMin()[wallNormal_];
                        wall[wallIdx].boundaryLayerThickness[posIdx] = this->problem_().bBoxMax()[wallNormal_] - this->problem_().bBoxMin()[wallNormal_] + eps_;
                        wall[wallIdx].isBBoxMinWall = true;
                        wall[wallIdx].sandGrainRoughness[posIdx] = GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMinSandGrainRoughness);
                    }
                    if (bboxMaxIsWall)
                    {
                        wall[wallIdx].wallPos[posIdx] = this->problem_().bBoxMax()[wallNormal_];
                        wall[wallIdx].boundaryLayerThickness[posIdx] = this->problem_().bBoxMin()[wallNormal_] - this->problem_().bBoxMax()[wallNormal_] - eps_;
                        wall[wallIdx].isBBoxMinWall = false;
                        wall[wallIdx].sandGrainRoughness[posIdx] = GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMaxSandGrainRoughness);
                    }
                }
                if (walls == 2 && wallIdx == 0)
                {
                    wall[0].wallPos[posIdx] = this->problem_().bBoxMin()[wallNormal_];
                    wall[1].wallPos[posIdx] = this->problem_().bBoxMax()[wallNormal_];
                    wall[0].boundaryLayerThickness[posIdx] = (wall[1].wallPos[posIdx] - wall[0].wallPos[posIdx]) / 2  + eps_;
                    wall[1].boundaryLayerThickness[posIdx] = - wall[0].boundaryLayerThickness[posIdx];
                    wall[0].isBBoxMinWall = true;
                    wall[1].isBBoxMinWall = false;
                    wall[0].sandGrainRoughness[posIdx] = GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMinSandGrainRoughness);
                    wall[1].sandGrainRoughness[posIdx] = GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMaxSandGrainRoughness);
                }
            wall[wallIdx].crossLength[posIdx] = wall[wallIdx].boundaryLayerThickness[posIdx];
            wall[wallIdx].boundaryLayerThicknessCalculated[posIdx] = wall[wallIdx].boundaryLayerThickness[posIdx]; // != 0, da kleinster wert gwünscht
            wall[wallIdx].viscousSublayerThicknessCalculated[posIdx] = 0; // == 0, da größter wert gewünscht
            for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                wall[wallIdx].maxVelocity[posIdx][dimIdx] = 0.0;
                wall[wallIdx].maxVelocityAbs[posIdx][dimIdx] = 0.0;
                wall[wallIdx].minVelocity[posIdx][dimIdx] = 0.0;
            }
            wall[wallIdx].fMax[posIdx] = 0.0;
            wall[wallIdx].yMax[posIdx] = 0.0;
            wall[wallIdx].isInterpolated[posIdx] = 0.0;
            wall[wallIdx].fluxValuesCount[posIdx] = 0.0;
            wall[wallIdx].wallValuesCount[posIdx] = 0.0;
            wall[wallIdx].headerWritten[posIdx] = false;
            }
        }
    }

    /*!
     * \brief Initializes the wall fluid properties with 0.
     */
    void resetWallFluidProperties()
    {
        for (int wallIdx = 0; wallIdx < walls; ++wallIdx)
            for (int posIdx = 0; posIdx < intervals; ++posIdx)
            {
                wall[wallIdx].wallDensity[posIdx] = 0.0;
                wall[wallIdx].wallKinematicViscosity[posIdx] = 0.0;
                wall[wallIdx].wallVelGrad[posIdx] = 0.0;
                wall[wallIdx].wallShearStress[posIdx] = 0.0;
            }
    }

    /*!
     * \brief Complete calculation for all relevant wall values.
     *
     * If the order is changed, errors may occur because of unknown or zero values
     */
    void updateWallProperties()
    {
        asImp_().resetWallProperties();
        asImp_().updateMaxFluxVars();
        asImp_().updateCrossLength();
        asImp_().resetWallFluidProperties();
        asImp_().updateWallFluidProperties();
        for (int wallIdx = 0; wallIdx < walls; ++wallIdx)
        {
            for (int posIdx = 0; posIdx < intervals; ++posIdx)
                if (wall[wallIdx].viscousSublayerThicknessCalculated[posIdx] > wall[wallIdx].boundaryLayerThicknessCalculated[posIdx])
                    wall[wallIdx].viscousSublayerThicknessCalculated[posIdx] = wall[wallIdx].boundaryLayerThicknessCalculated[posIdx];
            asImp_().interpolateWallProperties(wallIdx);
        }
    }

    /*!
     * \brief Get position index (interval section) a point belongs to.
     *
     * \param globalPos Global Position.
     */
    const int getPosIdx(const GlobalPosition &globalPos) const
    {
        int posIdx = int (intervals * (globalPos[flowNormal_] - this->problem_().bBoxMin()[flowNormal_]))
                         / (this->problem_().bBoxMax()[flowNormal_] - this->problem_().bBoxMin()[flowNormal_]);
        if (posIdx >= intervals)
            posIdx = intervals -1;
        return posIdx;
    }

    /*!
     * \brief Return wall a point belongs to.
     *
     * \param posIdx Position Index of current Global Position.
     * \param globalPos Global Position.
     */
    const int getWallIdx(const GlobalPosition &globalPos, const int posIdx) const
    {
        if (walls == 0)
            DUNE_THROW(Dune::NotImplemented, "Eddy viscosity models are not implemented for use without walls.");

        for (int wallIdx = 0; wallIdx < walls; ++wallIdx)
            if ((wall[wallIdx].isBBoxMinWall && globalPos[wallNormal_] < wall[wallIdx].wallPos[posIdx] + wall[wallIdx].boundaryLayerThickness[posIdx])
                 || (!wall[wallIdx].isBBoxMinWall && globalPos[wallNormal_] > wall[wallIdx].wallPos[posIdx] + wall[wallIdx].boundaryLayerThickness[posIdx]))
            {
                return wallIdx;
            }

        static bool alreadyPrintedError = false;
        static Scalar timePrintedError = 0.0;
        if (this->problem_().timeManager().time() > timePrintedError)
        {
            alreadyPrintedError = false;
        }

        if (!alreadyPrintedError)
        {
            Dune::dinfo << "info: point " << globalPos << " in interval " << posIdx
                        << " does not belong to wall -> now belongs to wall 0" << std::endl;
            alreadyPrintedError = true;
            timePrintedError = this->problem_().timeManager().time();
        }
        return 0;
    }


    /*!
     * \brief Return the distance to corresponding wall including roughness effects
     *
     * For ksPlus = 0, this is the normal wall distance.
     * For ksPlus > 0, an additional length is added to the real wall distance.
     * \param globalPos Global Position.
     * \param wallIdx Wall Index of current Global Position.
     * \param posIdx Position Index of current Global Position.
     */
    const Scalar distanceToWallRough(const GlobalPosition &globalPos, const int wallIdx, const int posIdx) const
    {
        if (surfaceRoughnessNotImplemented())
            { return distanceToWallReal (globalPos, wallIdx, posIdx); }
        if (GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMinSandGrainRoughness) < eps_&& wall[wallIdx].isBBoxMinWall)
            { return distanceToWallReal (globalPos, wallIdx, posIdx); }
        if (GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMaxSandGrainRoughness) < eps_&& !wall[wallIdx].isBBoxMinWall)
            { return distanceToWallReal (globalPos, wallIdx, posIdx); }

        Scalar ksPlus = wall[wallIdx].sandGrainRoughness[posIdx] * sqrt(wall[wallIdx].wallShearStress[posIdx] / wall[wallIdx].wallDensity[posIdx])
                       / wall[wallIdx].wallKinematicViscosity[posIdx];

        if (ksPlus > 2000)
        {
            std::cout << "info: equivalent sand grain roughness ks+=" << ksPlus << " at " << globalPos
                      << " is not in the valid range (ksPlus < 2000),"
                      << " for high ksPlus values the roughness function reaches a turning point."<< std::endl;
            DUNE_THROW(Dune::NotImplemented, "Unphysical roughness behavior.");
            ksPlus = 2000;
        }
        else if (ksPlus < 4.535)
        {
            Dune::dinfo << "info: equivalent sand grain roughness ks+=" << ksPlus << " at " << globalPos
                        << " is not in the valid range (ksPlus > 4.535) and now set to 0.0"<< std::endl;
            ksPlus = 0.0;
        }

        Scalar delta = 0.0; // delta is only changed for ksPlus, which are in valid range, otherwise smooth wall is assumed
        if (ksPlus >= 4.535)
        {
            delta = 0.9 * wall[wallIdx].wallKinematicViscosity[posIdx] / sqrt(wall[wallIdx].wallShearStress[posIdx] / wall[wallIdx].wallDensity[posIdx]) * (sqrt(ksPlus) - ksPlus * exp(- 1.0 * ksPlus / 6.0));
        }

        int sign = std::abs(globalPos[wallNormal_] - wall[wallIdx].wallPos[posIdx]) / (globalPos[wallNormal_] - wall[wallIdx].wallPos[posIdx]);
        return globalPos[wallNormal_] - wall[wallIdx].wallPos[posIdx] + sign * delta;
    }

    /*!
     * \brief Return the distance to corresponding wall.
     *
     * \param globalPos Global Position.
     * \param wallIdx Wall Index of current Global Position.
     * \param posIdx Position Index of current Global Position.
     */
    const Scalar distanceToWallReal(const GlobalPosition &globalPos, const int wallIdx, const int posIdx) const
    {   return globalPos[wallNormal_] - wall[wallIdx].wallPos[posIdx];  }

    /*!
     * \brief Return true if function for muInner should be used.
     *
     * \param globalPos Global Position.
     * \param posIdx Position Index of current Global Position.
     */
    const bool useViscosityInner(const GlobalPosition &globalPos, const int posIdx) const
    {
        for (int wallIdx = 0; wallIdx < walls; ++wallIdx)
            if ((wall[wallIdx].isBBoxMinWall && globalPos[wallNormal_] < wall[wallIdx].wallPos[posIdx] + wall[wallIdx].crossLength[posIdx])
                || (!wall[wallIdx].isBBoxMinWall && globalPos[wallNormal_] > wall[wallIdx].wallPos[posIdx] + wall[wallIdx].crossLength[posIdx]))
                return true;
        return false;
    }

    /*!
     * \brief Loop over all elements, used to calculate maximum flux values
     *        in the whole domain.
     */
    void updateMaxFluxVars()
    {
        FVElementGeometry fvGeometry;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();

        for (; eIt != eEndIt; ++eIt)
        {
            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt, fvGeometry);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false);

            for (int fIdx = 0; fIdx < fvGeometry.numScvf; ++fIdx)
            {

                FluxVariables fluxVars(this->problem_(),
                                    *eIt,
                                    fvGeometry,
                                    fIdx,
                                    elemVolVars,
                                    false);

                GlobalPosition globalPos = fvGeometry.subContVolFace[fIdx].ipGlobal;

                asImp_().calculateMaxFluxVars(fluxVars, globalPos);
            }
        }
    }

    /*!
     * \brief Update maximum values in the domain and
     *        set them to the corresponding wall.
     *
     * \param fluxVars Flux Variables of current element.
     * \param globalPos Global Position.
     */
    void calculateMaxFluxVars(const FluxVariables &fluxVars, const GlobalPosition globalPos)
    {
        int posIdx = getPosIdx(globalPos);
        int wallIdx = getWallIdx(globalPos, posIdx);
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            if (std::abs(wall[wallIdx].maxVelocity[posIdx][dimIdx]) < std::abs(fluxVars.velocity()[dimIdx]))
            {
                wall[wallIdx].maxVelocity[posIdx][dimIdx] = fluxVars.velocity()[dimIdx];
//                 // if the values in the middle should be set on both wall
//                 for (int wIdx = 0; wIdx < walls; ++wIdx)
//                     if (std::abs(distanceToWallReal(globalPos, wallIdx, posIdx)) < std::abs(wall[wIdx].boundaryLayerThickness[posIdx] + 1e-5))
//                         wall[wIdx].maxVelocity[posIdx][dimIdx] = fluxVars.velocity()[dimIdx];
                // set it as maxVelocityAbs
                if (std::abs(wall[wallIdx].maxVelocityAbs[posIdx][dimIdx]) < std::abs(fluxVars.velocity()[dimIdx]))
                    for (int wIdx = 0; wIdx < walls; ++wIdx)
                        wall[wIdx].maxVelocityAbs[posIdx][dimIdx] = fluxVars.velocity()[dimIdx];
                wall[wallIdx].fluxValuesCount[posIdx]++;
            }
            if (std::abs(wall[wallIdx].minVelocity[posIdx][dimIdx]) > std::abs(fluxVars.velocity()[dimIdx]))
                wall[wallIdx].minVelocity[posIdx][dimIdx] = fluxVars.velocity()[dimIdx];
        }

        // fMax and yMax
        if (wall[wallIdx].fMax[posIdx] < fluxVars.fz())
        {
            wall[wallIdx].fMax[posIdx] = fluxVars.fz();
            wall[wallIdx].yMax[posIdx] = distanceToWallReal(globalPos, wallIdx, posIdx);
//            // if the values in the middle should be set on both wall
//            for (int wIdx = 0; wIdx < walls; ++wIdx)
//                if (std::abs(distanceToWall(globalPos, wIdx, posIdx)) < std::abs(wall[wIdx].boundaryLayerThickness[posIdx] + 1e-4))
//                    {
//                        wall[wIdx].fMax[posIdx] = fluxVars.fz();
//                        wall[wIdx].yMax[posIdx] = distanceToWall(globalPos, wallIdx, posIdx);
//                    }
        }
    }

    /*!
     * \brief Get maximum values of the f(z) function in the Baldwin-Lomax model.
     */
    void updateCrossLength()
    {
        FVElementGeometry fvGeometry;
        ElementBoundaryTypes elemBcTypes;
        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();

        for (; eIt != eEndIt; ++eIt)
        {
            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt, fvGeometry);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false);

            IntersectionIterator isIt = this->gridView_().ibegin(*eIt);
            const IntersectionIterator &endIt = this->gridView_().iend(*eIt);

            for (; isIt != endIt; ++isIt)
            {
                int fIdx = isIt->indexInInside();
                FluxVariables fluxVars(this->problem_(),
                                                    *eIt,
                                                    fvGeometry,
                                                    fIdx,
                                                    elemVolVars,
                                                    false);

                GlobalPosition globalPos = fvGeometry.subContVolFace[fIdx].ipGlobal;
                int posIdx = getPosIdx(globalPos);
                int wallIdx = getWallIdx(globalPos, posIdx);

                // muCross
                if (fluxVars.dynamicEddyViscosityOuter() < fluxVars.dynamicEddyViscosityInner() && useViscosityInner(globalPos, posIdx))
                    wall[wallIdx].crossLength[posIdx] = distanceToWallReal(globalPos, wallIdx, posIdx) - eps_ * std::abs(distanceToWallReal(globalPos, wallIdx, posIdx) / distanceToWallReal(globalPos, wallIdx, posIdx));

                if (std::abs(fluxVars.velocity()[flowNormal_]) >= 0.99 * std::abs(wall[wallIdx].maxVelocity[posIdx][flowNormal_])
                    && std::abs(wall[wallIdx].boundaryLayerThicknessCalculated[posIdx]) > std::abs(distanceToWallReal(globalPos, wallIdx, posIdx)))
                    wall[wallIdx].boundaryLayerThicknessCalculated[posIdx] = distanceToWallReal(globalPos, wallIdx, posIdx);

                if (wall[wallIdx].maxVelocity[posIdx][flowNormal_] * globalPos[flowNormal_] / fluxVars.kinematicViscosity() < 2300)
                    wall[wallIdx].viscousSublayerThicknessCalculated[posIdx] = wall[wallIdx].boundaryLayerThicknessCalculated[posIdx];
                else if ((fluxVars.velocity()[flowNormal_] / fluxVars.frictionVelocityWall() <= 5.0)
                         && (std::abs(wall[wallIdx].viscousSublayerThicknessCalculated[posIdx])
                             < std::abs(distanceToWallReal(globalPos, wallIdx, posIdx))))
                    wall[wallIdx].viscousSublayerThicknessCalculated[posIdx] = distanceToWallReal(globalPos, wallIdx, posIdx);
            }
        }
    }

    /*!
     * \brief Loop over all elements to update the values at the wall.
     */
    void updateWallFluidProperties()
    {
        FVElementGeometry fvGeometry;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();

        for (; eIt != eEndIt; ++eIt)
        {
            fvGeometry.update(this->gridView_(), *eIt);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false);

            const ReferenceElement &refElement = ReferenceElements::general(eIt->geometry().type());
            IntersectionIterator isIt = this->gridView_().ibegin(*eIt);
            const IntersectionIterator &endIt = this->gridView_().iend(*eIt);
            for (; isIt != endIt; ++isIt)
            {
                // handle only faces on the boundary
                if (!isIt->boundary())
                    continue;

                // Assemble the boundary for all vertices of the current
                // face
                int fIdx = isIt->indexInInside();
                int numFaceVerts = refElement.size(fIdx, 1, dim);
                for (int faceVertIdx = 0;
                     faceVertIdx < numFaceVerts;
                     ++faceVertIdx)
                {
                    int boundaryFaceIdx = fvGeometry.boundaryFaceIndex(fIdx, faceVertIdx);




                    const FluxVariables boundaryVars(this->problem_(),
                                                     *eIt,
                                                     fvGeometry,
                                                     boundaryFaceIdx,
                                                     elemVolVars,
                                                     true);

                    GlobalPosition globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;
                    if (
                        globalPos[wallNormal_] > this->problem_().bBoxMin()[wallNormal_] + 1e-15
                        && globalPos[wallNormal_] < this->problem_().bBoxMax()[wallNormal_] - 1e-15
                        )
                    continue;

                    asImp_().calculateWallFluidProperties(boundaryVars, globalPos);
                } // end loop over intersections
            } // end loop over element vertices
        }
    }

    /*!
     * \brief Calculate / average the values at the wall.
     *
     * \param boundaryVars Flux Variables at boundary segment.
     * \param globalPos Global Position.
     */
    void calculateWallFluidProperties(const FluxVariables &boundaryVars, const GlobalPosition &globalPos)
    {
        int posIdx = getPosIdx(globalPos);
        int wallIdx = getWallIdx(globalPos, posIdx);
        if (globalPos[wallNormal_] > wall[wallIdx].wallPos[posIdx] - 1e-8
            && globalPos[wallNormal_] < wall[wallIdx].wallPos[posIdx] + 1e-8)
        {
            wall[wallIdx].wallValuesCount[posIdx] += 1;
            wall[wallIdx].wallDensity[posIdx] =
                (wall[wallIdx].wallDensity[posIdx] * (wall[wallIdx].wallValuesCount[posIdx] - 1) + boundaryVars.density())
                / wall[wallIdx].wallValuesCount[posIdx];
            wall[wallIdx].wallKinematicViscosity[posIdx] =
                (wall[wallIdx].wallKinematicViscosity[posIdx] * (wall[wallIdx].wallValuesCount[posIdx] - 1) + boundaryVars.kinematicViscosity())
                / wall[wallIdx].wallValuesCount[posIdx];
            wall[wallIdx].wallVelGrad[posIdx] =
                (boundaryVars.velocityGrad()[flowNormal_][wallNormal_] * (wall[wallIdx].wallValuesCount[posIdx] - 1) + boundaryVars.velocityGrad()[flowNormal_][wallNormal_])
                / wall[wallIdx].wallValuesCount[posIdx];
            wall[wallIdx].wallShearStress[posIdx] =
                (wall[wallIdx].wallShearStress[posIdx] * (wall[wallIdx].wallValuesCount[posIdx] - 1)
                    + std::abs(boundaryVars.velocityGrad()[flowNormal_][wallNormal_]) * boundaryVars.dynamicViscosity())
                / wall[wallIdx].wallValuesCount[posIdx];
        }
    }


    /*!
     * \brief Find points with given values and start interpolation.
     *
     * \param wallIdx Wall Index of current Global Position.
     */
    const void interpolateWallProperties(const int wallIdx)
    {
        const int startInterpolation = 0;
        const int endInterpolation = intervals;

        for (int posIdx = startInterpolation; posIdx < endInterpolation; ++posIdx)
        {
            if (wall[wallIdx].fluxValuesCount[posIdx] == 0)
            {
                int prevIdx = posIdx;
                int nextIdx = posIdx;
                // Getting previous value, if 0 is reached (and tau is still < eps_), get next value
                for (prevIdx = posIdx-1; prevIdx >= startInterpolation; --prevIdx)
                    if (wall[wallIdx].fluxValuesCount[prevIdx] != 0)
                        break;
                if (prevIdx < startInterpolation) // interpolation at x=0, prevIdx is -1
                    for (prevIdx = posIdx+1; prevIdx < endInterpolation; ++prevIdx)
                        if (wall[wallIdx].fluxValuesCount[prevIdx] != 0)
                            break;
                if (prevIdx == endInterpolation && this->problem_().timeManager().time() > this->problem_().timeManager().timeStepSize())
                {
                    Dune::dinfo << "info: for posIdx " << posIdx << "on wall " << wallIdx
                                << "there are no fluxValues for interpolation." << std::endl;
                    break;
                }

                // Getting next value, if intervals is reached get prev value
                for (nextIdx = posIdx+1; nextIdx < endInterpolation; ++nextIdx)
                    if (wall[wallIdx].fluxValuesCount[nextIdx] != 0 && nextIdx != prevIdx)
                        break;
                if (nextIdx == endInterpolation)
                    for (nextIdx = posIdx-1; nextIdx >= startInterpolation; --nextIdx)
                        if (wall[wallIdx].fluxValuesCount[nextIdx] != 0 && nextIdx != prevIdx)
                            break;

                asImp_().doInterpolationFluxValues(wallIdx, posIdx, prevIdx, nextIdx);
            }


            if (wall[wallIdx].wallValuesCount[posIdx] == 0) // || wall[wallIdx].wallValuesCount[posIdx] > 50)
            {
                int prevIdx = posIdx;
                int nextIdx = posIdx;
                // Getting previous value, if 0 is reached (and tau is still < eps_), get next value
                for (prevIdx = posIdx-1; prevIdx >= startInterpolation; --prevIdx)
                    if (wall[wallIdx].wallValuesCount[prevIdx] != 0)
                        break;
                if (prevIdx < startInterpolation) // interpolation at x=0, prevIdx is -1
                    for (prevIdx = posIdx+1; prevIdx < endInterpolation; ++prevIdx)
                        if (wall[wallIdx].wallValuesCount[prevIdx] != 0)
                            break;
                if (prevIdx == endInterpolation && this->problem_().timeManager().time() > this->problem_().timeManager().timeStepSize())
                {
                    Dune::dinfo << "info: for posIdx " << posIdx << "on wall " << wallIdx
                                << "there are no wallValues for interpolation." << std::endl;
                    break;
                }

                // Getting next value, if intervals is reached get prev value
                for (nextIdx = posIdx+1; nextIdx < endInterpolation; ++nextIdx)
                    if (wall[wallIdx].wallValuesCount[nextIdx] != 0 && nextIdx != prevIdx)
                        break;
                if (nextIdx == endInterpolation)
                    for (nextIdx = posIdx-1; nextIdx >= startInterpolation; --nextIdx)
                        if (wall[wallIdx].wallValuesCount[nextIdx] != 0 && nextIdx != prevIdx)
                            break;

                asImp_().doInterpolationWallValues(wallIdx, posIdx, prevIdx, nextIdx);
            }
        }
    }

    /*!
     * \brief Interpolate flux Values, so that flux related Properties
     *        are given in every Interval.
     *
     * \param wallIdx Wall Index for interpolation.
     * \param posIdx Position Index for interpolation (no given value).
     * \param prevIdx Position Index with value.
     * \param nextIdx Position Index with value.
     */
    const void doInterpolationFluxValues(const int wallIdx, const int posIdx, const int prevIdx, const int nextIdx)
    {
        wall[wallIdx].boundaryLayerThickness[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].boundaryLayerThickness[prevIdx], nextIdx, wall[wallIdx].boundaryLayerThickness[nextIdx]);
        wall[wallIdx].boundaryLayerThicknessCalculated[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].boundaryLayerThicknessCalculated[prevIdx], nextIdx, wall[wallIdx].boundaryLayerThicknessCalculated[nextIdx]);
        wall[wallIdx].viscousSublayerThicknessCalculated[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].viscousSublayerThicknessCalculated[prevIdx], nextIdx, wall[wallIdx].viscousSublayerThicknessCalculated[nextIdx]);
        wall[wallIdx].crossLength[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].crossLength[prevIdx], nextIdx, wall[wallIdx].crossLength[nextIdx]);
        wall[wallIdx].maxVelocity[posIdx][0] = interpolation(posIdx, prevIdx, wall[wallIdx].maxVelocity[prevIdx][0], nextIdx, wall[wallIdx].maxVelocity[nextIdx][0]);
        wall[wallIdx].maxVelocity[posIdx][1] = interpolation(posIdx, prevIdx, wall[wallIdx].maxVelocity[prevIdx][1], nextIdx, wall[wallIdx].maxVelocity[nextIdx][1]);
        wall[wallIdx].minVelocity[posIdx][0] = interpolation(posIdx, prevIdx, wall[wallIdx].minVelocity[prevIdx][0], nextIdx, wall[wallIdx].minVelocity[nextIdx][0]);
        wall[wallIdx].minVelocity[posIdx][1] = interpolation(posIdx, prevIdx, wall[wallIdx].minVelocity[prevIdx][1], nextIdx, wall[wallIdx].minVelocity[nextIdx][1]);
        wall[wallIdx].maxVelocityAbs[posIdx][0] = interpolation(posIdx, prevIdx, wall[wallIdx].maxVelocityAbs[prevIdx][0], nextIdx, wall[wallIdx].maxVelocityAbs[nextIdx][0]);
        wall[wallIdx].maxVelocityAbs[posIdx][1] = interpolation(posIdx, prevIdx, wall[wallIdx].maxVelocityAbs[prevIdx][1], nextIdx, wall[wallIdx].maxVelocityAbs[nextIdx][1]);
        wall[wallIdx].fMax[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].fMax[prevIdx], nextIdx, wall[wallIdx].fMax[nextIdx]);
        wall[wallIdx].yMax[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].yMax[prevIdx], nextIdx, wall[wallIdx].yMax[nextIdx]);
        wall[wallIdx].isInterpolated[posIdx] += 1;
    }

    /*!
     * \brief Interpolate wall Values, so that wall related Properties
     *        are given in every Interval.
     *
     * \param wallIdx Wall Index for interpolation.
     * \param posIdx Position Index for interpolation (no given value).
     * \param prevIdx Position Index with value.
     * \param nextIdx Position Index with value.
     */
    const void doInterpolationWallValues(const int wallIdx, const int posIdx, const int prevIdx, const int nextIdx)
    {
        wall[wallIdx].wallDensity[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].wallDensity[prevIdx], nextIdx, wall[wallIdx].wallDensity[nextIdx]);
        wall[wallIdx].wallKinematicViscosity[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].wallKinematicViscosity[prevIdx], nextIdx, wall[wallIdx].wallKinematicViscosity[nextIdx]);
        wall[wallIdx].wallVelGrad[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].wallVelGrad[prevIdx], nextIdx, wall[wallIdx].wallVelGrad[nextIdx]);
        wall[wallIdx].wallShearStress[posIdx] = interpolation(posIdx, prevIdx, wall[wallIdx].wallShearStress[prevIdx], nextIdx, wall[wallIdx].wallShearStress[nextIdx]);
        wall[wallIdx].isInterpolated[posIdx] += 2;
    }

    /*!
     * \brief Mathematical linear interpolation routine.
     *
     * Return the linear interpolated value at any point between to values.
     *
     * \param position Position for which interpolation is made.
     * \param prev First Position for which the value is known.
     * \param prevValue Known value at prev position.
     * \param next Second Position for which the value is known.
     * \param nextValue Known value at next position.
     */
    const Scalar interpolation(const Scalar position, const Scalar prev, const Scalar prevValue, const Scalar next, const Scalar nextValue)
    {
        return (prevValue + (nextValue - prevValue) / (next - prev) * (position - prev));
    }

    /*!
     * \brief Returns whether the actual eddy viscosity model includes surface roughness approach.
     *
     * Surface roughness is not included in the Baldwin Lomax model
     */
    const bool surfaceRoughnessNotImplemented() const
    {
        switch (GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyViscosityModel))
        {
            case EddyViscosityIndices::noEddyViscosityModel: // 0
                return true;
                break;
            default:
                return false;
        }
    }

    /*!
     * \brief Returns the name of the used eddy viscosity model.
     */
    const char *eddyViscosityModelName() const
    {
        switch (GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyViscosityModel))
        {
            case EddyViscosityIndices::noEddyViscosityModel: // 0
                return "noEddyViscosityModel";
                break;
            case EddyViscosityIndices::prandtl: // 1
                return "prandtl";
                break;
            case EddyViscosityIndices::modifiedVanDriest: // 2
                return "modifiedVanDriest";
                break;
            case EddyViscosityIndices::baldwinLomax: // 3
                return "baldwinLomax";
                break;
            default:
                DUNE_THROW(Dune::NotImplemented, "This eddy viscosity model is not implemented.");
        }
    }


protected:
    //! Current implementation.
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    //! Current implementation.
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

private:
    const int flowNormal_;
    const int wallNormal_;
    Scalar writeAllSCVData_;
    Scalar scvDataCoordinateDelta_;
    Scalar eps_;
};

}

#include "zeroeqpropertydefaults.hh"

#endif // DUMUX_ZEROEQ_MODEL_HH
