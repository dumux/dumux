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
 * \brief  Finite volume discretization of a saturation transport equation.
 */
#ifndef DUMUX_FVSATURATION2P_HH
#define DUMUX_FVSATURATION2P_HH

#include <dumux/porousmediumflow/2p/sequential/transport/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/transport.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief The finite volume discretization of a saturation transport equation.
 *
 * This model solves equations of the form
 *
 *  \f[
 *  \phi \frac{\partial (\varrho_\alpha S_\alpha)}{\partial t} + \text{div}\, (\varrho_\alpha \boldsymbol{v_\alpha}) = q_\alpha,
 *  \f]
 *
 *  where \f$ S_\alpha \f$ is the saturation of phase \f$\alpha \in \{ w, n \}\f$
 *  and \f$ \boldsymbol v_\alpha \f$ is the phase velocity defined by
 *  the multi-phase Darcy equation.
 *  If a phase velocity is reconstructed from the pressure solution it can be directly inserted into
 *  the previous equation. In the incompressible case the equation is further divided by the phase density
 *  \f$ \varrho_\alpha \f$. If a total velocity is reconstructed the saturation equation is reformulated into:
 *
 * \f[
 *  \phi \frac{\partial S_w}{\partial t} + f_w \text{div}\, \boldsymbol{v}_{t} + f_w \lambda_n \boldsymbol{K}\left(\textbf{grad}\,
 *  p_c - (\varrho_n-\varrho_w) {\textbf g} \right)= q_\alpha,
 * \f]
 * to get a wetting phase saturation or
 * \f[
 * \phi \frac{\partial S_n}{\partial t} + f_n \text{div}\, \boldsymbol{v}_{t} - f_n \lambda_w \boldsymbol{K}\left(\textbf{grad}\,
 * p_c - (\varrho_n-\varrho_w) {\textbf g} \right)= q_\alpha,
 * \f]
 * if the nonwetting phase saturation is the primary transport variable.
 *
 *  The total velocity formulation is only implemented for incompressible fluids and \f$ f_\alpha \f$
 *  is the fractional flow function, \f$ \lambda_\alpha \f$ is the mobility, \f$ \boldsymbol K \f$
 *  the absolute permeability tensor,\f$ p_c \f$ the capillary pressure, \f$ \varrho_\alpha \f$ the phase density,
 *  \f$ {\textbf g} \f$ the gravitational acceleration vector, and \f$ q_\alpha \f$ the source term.
 *
 *
 *  In the IMPES models the default setting is:
 *
 * formulation: \f$ p_w \f$ - \f$ S_w \f$ (Property: \a Formulation defined as \a SequentialTwoPCommonIndices::pwsw)
 *
 * compressibility: disabled (Property: \a EnableCompressibility set to \a false)
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVSaturation2P: public FVTransport<TypeTag>
{
    using ParentType = FVTransport<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::TransportModel>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Velocity = GetPropType<TypeTag, Properties::Velocity>;
    using CapillaryFlux = GetPropType<TypeTag, Properties::CapillaryFlux>;
    using GravityFlux = GetPropType<TypeTag, Properties::GravityFlux>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using CellData = GetPropType<TypeTag, Properties::CellData>;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pGlobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNw,
        vt = Indices::velocityTotal,
        sw = Indices::saturationW,
        sn = Indices::saturationNw
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        saturationIdx = Indices::saturationIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };

    using TransportSolutionType = GetPropType<TypeTag, Properties::TransportSolutionType>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dim>;

protected:
    CapillaryFlux& capillaryFlux()
    {
        return *capillaryFlux_;
    }

    const CapillaryFlux& capillaryFlux() const
    {
        return *capillaryFlux_;
    }

    GravityFlux& gravityFlux()
    {
        return *gravityFlux_;
    }

    const GravityFlux& gravityFlux() const
    {
        return *gravityFlux_;
    }

public:
    Velocity& velocity()
    {
        return *velocity_;
    }

    Velocity& velocity() const
    {
        return *velocity_;
    }

    //! Function which calculates the flux update
    void getFlux(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    //! Function which calculates the boundary flux update
    void getFluxOnBoundary(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    //! Function which calculates the source update
    void getSource(Scalar& update, const Element& element, CellData& cellDataI);

    //! Sets the initial solution
    void initialize();

    //! Update the values of the material laws and constitutive relations.
    void updateMaterialLaws();


    /*!
     * \brief Writes the current values of the primary transport variable
     *  into the <tt>transportedQuantity</tt>-vector (comes as function argument)
     *
     * \copydetails FVTransport::getTransportedQuantity(TransportSolutionType&)
     */
    void getTransportedQuantity(TransportSolutionType& transportedQuantity)
    {
        int size = problem_.gridView().size(0);
        transportedQuantity.resize(size);
        for (int i = 0; i < size; i++)
        {
            switch (saturationType_)
            {
            case sw:
            {
                transportedQuantity[i] = problem_.variables().cellData(i).saturation(wPhaseIdx);
                break;
            }
            case sn:
            {
                transportedQuantity[i] = problem_.variables().cellData(i).saturation(nPhaseIdx);
                break;
            }
            }
        }
    }

    /*!
     * \brief Writes the current values of the primary transport variable into the variable container
     *
     * \copydetails FVTransport::setTransportedQuantity(TransportSolutionType&)
     *
     */
    void setTransportedQuantity(TransportSolutionType& transportedQuantity)
    {
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            switch (saturationType_)
            {
            case sw:
            {
                Scalar sat = transportedQuantity[i];

                    cellData.setSaturation(wPhaseIdx, sat);
                    cellData.setSaturation(nPhaseIdx, 1 - sat);
                break;
            }
            case sn:
            {
                Scalar sat = transportedQuantity[i];

                    cellData.setSaturation(wPhaseIdx,1 -sat);
                    cellData.setSaturation(nPhaseIdx, sat);
                break;
            }
            }
        }
    }

    /*!
     * \brief Check if saturation is in physical range.
     *
     * \tparam DataEntry Data class
     * \param entry Entry which is checked
     */
    template<class DataEntry>
    bool inPhysicalRange(DataEntry& entry)
    {
        return (entry > -1e-6 && entry < 1.0 + 1e-6);
    }

    /*!
     * \brief Updates the primary transport variable.
     *
     * \copydetails FVTransport::updateTransportedQuantity(TransportSolutionType&)
     */
    void updateTransportedQuantity(TransportSolutionType& updateVec)
    {
        if (this->enableLocalTimeStepping())
            this->innerUpdate(updateVec);
        else
            asImp_().updateSaturationSolution(updateVec);
    }

    /*!
     * \brief Updates the primary transport variable.
     *
     * \param updateVec Vector containing the global update
     * \param dt time step for update
     */
    void updateTransportedQuantity(TransportSolutionType& updateVec, Scalar dt)
    {
        asImp_().updateSaturationSolution(updateVec, dt);
    }

    /*!
     * \brief Globally updates the saturation solution
     *
     * \param updateVec Vector containing the global update.
     */
    void updateSaturationSolution(TransportSolutionType& updateVec)
    {
        Scalar dt = problem_.timeManager().timeStepSize();
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            asImp_().updateSaturationSolution(i, updateVec[i][0], dt);
        }
    }

    /*!
     * \brief Globally updates the saturation solution
     *
     * \param updateVec Vector containing the global update.
     * \param dt time step for update
     */
    void updateSaturationSolution(TransportSolutionType& updateVec, Scalar dt)
    {
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            asImp_().updateSaturationSolution(i, updateVec[i][0], dt);
        }
    }

    /*!
     * \brief Updates the saturation solution of a cell
     *
     * Calculates secondary saturation variables and stores saturations.
     *
     * \param eIdxGlobal Global cell index
     * \param update Cell saturation update
     * \param dt Current time step
     */
    void updateSaturationSolution(int eIdxGlobal, Scalar update, Scalar dt)
    {
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        switch (saturationType_)
        {
        case sw:
        {
            Scalar sat = cellData.saturation(wPhaseIdx) + dt*update;

                cellData.setSaturation(wPhaseIdx, sat);
                cellData.setSaturation(nPhaseIdx, 1 - sat);
            break;
        }
        case sn:
        {
            Scalar sat = cellData.saturation(nPhaseIdx) + dt*update;

                cellData.setSaturation(wPhaseIdx,1 -sat);
                cellData.setSaturation(nPhaseIdx, sat);
            break;
        }
        }
    }

    /*!
     * \brief Adds saturation output to the output file
     *
     * Adds the phase saturation to the output.
     * If the velocity is calculated in the transport model it is also added to the output.
     * If the VtkOutputLevel is equal to zero (default) only primary variables are written,
     * if it is larger than zero also secondary variables are written.
     *
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        TransportSolutionType *saturation = writer.allocateManagedBuffer(size);
        TransportSolutionType *saturationSecond = 0;

        if (vtkOutputLevel_ > 0)
        {
            saturationSecond = writer.allocateManagedBuffer(size);
        }

        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);

            if (saturationType_ == sw)
            {
                (*saturation)[i] = cellData.saturation(wPhaseIdx);
                if (vtkOutputLevel_ > 0)
                {
                    (*saturationSecond)[i] = cellData.saturation(nPhaseIdx);
                }
            }
            else if (saturationType_ == sn)
            {
                (*saturation)[i] = cellData.saturation(nPhaseIdx);
                if (vtkOutputLevel_ > 0)
                {
                    (*saturationSecond)[i] = cellData.saturation(wPhaseIdx);
                }
            }
        }

        if (saturationType_ == sw)
        {
            writer.attachCellData(*saturation, "wetting saturation");
            if (vtkOutputLevel_ > 0)
            {
                writer.attachCellData(*saturationSecond, "nonwetting saturation");
            }
        }
        else if (saturationType_ == sn)
        {
            writer.attachCellData(*saturation, "nonwetting saturation");
            if (vtkOutputLevel_ > 0)
            {
                writer.attachCellData(*saturationSecond, "wetting saturation");
            }
        }

        if (velocity().calculateVelocityInTransport())
        velocity().addOutputVtkFields(writer);

        return;
    }

    /*!
     * \brief  Function for serialization of the primary transport variable.
     *
     *\copydetails FVTransport::serializeEntity(std::ostream&,const Element&)
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int eIdxGlobal = problem_.variables().index(element);

        Scalar sat = 0.0;
        switch (saturationType_)
        {
        case sw:
            sat = problem_.variables().cellData(eIdxGlobal).saturation(wPhaseIdx);
        break;
        case sn:
            sat = problem_.variables().cellData(eIdxGlobal).saturation(nPhaseIdx);
            break;
        }

        outstream << sat;
    }

    /*!
     * \brief  Function for deserialization of the primary transport variable.
     *
     *\copydetails FVTransport::deserializeEntity(std::istream&,const Element&)
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int eIdxGlobal = problem_.variables().index(element);

        Scalar sat = 0.;
        instream >> sat;

        switch (saturationType_)
        {
        case sw:
            problem_.variables().cellData(eIdxGlobal).setSaturation(wPhaseIdx, sat);
            problem_.variables().cellData(eIdxGlobal).setSaturation(nPhaseIdx, 1-sat);
        break;
        case sn:
            problem_.variables().cellData(eIdxGlobal).setSaturation(nPhaseIdx, sat);
            problem_.variables().cellData(eIdxGlobal).setSaturation(wPhaseIdx, 1-sat);
            break;
        }
    }


    /*!
     * \brief Constructs a FVSaturation2P object
     *
     * \param problem A problem class object
     */
    FVSaturation2P(Problem& problem) :
            ParentType(problem), problem_(problem), threshold_(1e-6),
            switchNormals_(getParam<bool>("Impet.SwitchNormals"))
    {
        if (compressibility_ && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented,
                    "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (saturationType_ != sw && saturationType_ != sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (pressureType_ != pw && pressureType_ != pn && pressureType_ != pGlobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }

        capillaryFlux_ = std::make_shared<CapillaryFlux>(problem);
        gravityFlux_ = std::make_shared<GravityFlux>(problem);
        velocity_ = std::make_shared<Velocity>(problem);

        vtkOutputLevel_ = getParam<int>("Vtk.OutputLevel");
        porosityThreshold_ = getParam<Scalar>("Impet.PorosityThreshold");
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    Problem& problem_;
    std::shared_ptr<Velocity> velocity_;
    std::shared_ptr<CapillaryFlux> capillaryFlux_;
    std::shared_ptr<GravityFlux> gravityFlux_;

    int vtkOutputLevel_;
    Scalar porosityThreshold_;

    static const bool compressibility_ = getPropValue<TypeTag, Properties::EnableCompressibility>();
    static const int saturationType_ = getPropValue<TypeTag, Properties::SaturationFormulation>();
    static const int velocityType_ = getPropValue<TypeTag, Properties::VelocityFormulation>();
    static const int pressureType_ = getPropValue<TypeTag, Properties::PressureFormulation>();

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    const Scalar threshold_;
    const bool switchNormals_;
};

/*!
 * \brief Function which calculates the flux update
 *
 * \copydetails FVTransport::getFlux(Scalar&,const Intersection&,CellData&)
 *
 * If a total velocity formulation is used this functions calculates not only the advective flux
 * but also fluxes due to gravity and capillary diffusion.
 * These have to be defined separately as implementation of a DiffusivePart or ConvectivePart
 * (e.g. GravityPart / CapillaryDiffusion ) and added to the property system via properties
 * <tt>CapillaryFlux</tt> and <tt>GravityFlux</tt>.
 */
template<class TypeTag>
void FVSaturation2P<TypeTag>::getFlux(Scalar& update, const Intersection& intersection, CellData& cellDataI)
{
    auto elementI = intersection.inside();
    auto elementJ = intersection.outside();

    const CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(elementJ));

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = elementI.geometry().center();
    const GlobalPosition& globalPosJ = elementJ.geometry().center();

    // cell volume, assume linear map here
    Scalar volume = elementI.geometry().volume();
    using std::max;
    Scalar porosity = max(problem_.spatialParams().porosity(elementI), porosityThreshold_);

    if (compressibility_)
    {
        viscosity_[wPhaseIdx] = cellDataI.viscosity(wPhaseIdx);
        viscosity_[nPhaseIdx] = cellDataI.viscosity(nPhaseIdx);
    }

    // local number of faces
    int isIndex = intersection.indexInInside();

    GlobalPosition unitOuterNormal = intersection.centerUnitOuterNormal();
    if (switchNormals_)
        unitOuterNormal *= -1.0;

    Scalar faceArea = intersection.geometry().volume();

    if (velocity().calculateVelocityInTransport() && !cellDataI.fluxData().haveVelocity(isIndex))
        velocity().calculateVelocity(intersection, cellDataI);

    //get velocity*normalvector*facearea/(volume*porosity)
    Scalar factorW = (cellDataI.fluxData().velocity(wPhaseIdx, isIndex) * unitOuterNormal) * faceArea;
    Scalar factorNw = (cellDataI.fluxData().velocity(nPhaseIdx, isIndex) * unitOuterNormal) * faceArea;
    Scalar factorTotal = factorW + factorNw;

    // distance vector between barycenters
    GlobalPosition distVec = globalPosJ - globalPosI;
    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    bool takeNeighbor = (elementI.level() < elementJ.level());
    //get phase potentials
    bool upwindWI =
            (takeNeighbor) ? !cellDataJ.fluxData().isUpwindCell(wPhaseIdx, intersection.indexInOutside()) :
                    cellDataI.fluxData().isUpwindCell(wPhaseIdx, isIndex);
    bool upwindNwI =
            (takeNeighbor) ? !cellDataJ.fluxData().isUpwindCell(nPhaseIdx, intersection.indexInOutside()) :
                    cellDataI.fluxData().isUpwindCell(nPhaseIdx, isIndex);

    Scalar lambdaW = 0;
    Scalar lambdaNw = 0;

    //upwinding of lambda dependend on the phase potential gradients
    if (upwindWI)
    {
        lambdaW = cellDataI.mobility(wPhaseIdx);
        if (compressibility_)
        {
            lambdaW /= cellDataI.density(wPhaseIdx);
        } //divide by density because lambda is saved as lambda*density
    }
    else
    {
        lambdaW = cellDataJ.mobility(wPhaseIdx);
        if (compressibility_)
        {
            lambdaW /= cellDataJ.density(wPhaseIdx);
        } //divide by density because lambda is saved as lambda*density
    }

    if (upwindNwI)
    {
        lambdaNw = cellDataI.mobility(nPhaseIdx);
        if (compressibility_)
        {
            lambdaNw /= cellDataI.density(nPhaseIdx);
        } //divide by density because lambda is saved as lambda*density
    }
    else
    {
        lambdaNw = cellDataJ.mobility(nPhaseIdx);
        if (compressibility_)
        {
            lambdaNw /= cellDataJ.density(nPhaseIdx);
        } //divide by density because lambda is saved as lambda*density
    }

    switch (velocityType_)
    {
    case vt:
    {
        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                            viscosity_[nPhaseIdx], factorTotal, intersection);

        //determine phase saturations from primary saturation variable
        Scalar satWI = cellDataI.saturation(wPhaseIdx);
        Scalar satWJ = cellDataJ.saturation(wPhaseIdx);

        Scalar pcI = cellDataI.capillaryPressure();
        Scalar pcJ = cellDataJ.capillaryPressure();

        // calculate the saturation gradient
        GlobalPosition pcGradient = unitOuterNormal;
        pcGradient *= (pcI - pcJ) / dist;

        // get the diffusive part
        DimVector flux(0.);
        capillaryFlux().getFlux(flux, intersection, satWI, satWJ, pcGradient);
        Scalar capillaryFlux = (flux * unitOuterNormal * faceArea);

        flux = 0.0;
        gravityFlux().getFlux(flux, intersection, satWI, satWJ);
        Scalar gravityFlux =  (flux * unitOuterNormal * faceArea);

        switch (saturationType_)
        {
        case sw:
        {
            //vt*fw
            factorTotal *= lambdaW / (lambdaW + lambdaNw);
            break;
        }
        case sn:
        {
            //vt*fn
            factorTotal *= lambdaNw / (lambdaW + lambdaNw);
            capillaryFlux *= -1;
            gravityFlux *= -1;
            break;
        }
        }
        factorTotal -= capillaryFlux;
        factorTotal += gravityFlux;

        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                            viscosity_[nPhaseIdx], 10 * capillaryFlux, intersection);
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                            viscosity_[nPhaseIdx], 10 * gravityFlux, intersection);

        break;
    }
    default:
    {
        if (compressibility_)
        {
            factorW /= cellDataI.density(wPhaseIdx);
            factorNw /= cellDataI.density(nPhaseIdx);
        }

        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                            viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                            viscosity_[nPhaseIdx], factorNw, intersection, nPhaseIdx);

        break;
    }
    }

    switch (velocityType_)
    {
    case vt:
        update -= factorTotal / (volume * porosity); //-:v>0, if flow leaves the cell
        break;
    default:
        switch (saturationType_)
        {
        case sw:
            update -= factorW / (volume * porosity);//-:v>0, if flow leaves the cell
            break;
        case sn:
            update -= factorNw / (volume * porosity); //-:v>0, if flow leaves the cell
            break;
        }
        break;
    }
}

/*!
 * \brief Function which calculates the boundary flux update
 *
 * \copydetails FVTransport::getFluxOnBoundary(Scalar&,const Intersection&,CellData&)
 *
 * Dirichlet boundary condition is a phase saturation depending on the formulation (\f$ S_w \f$ (default) or \f$ S_n \f$),
 * Neumann boundary condition are phase mass fluxes (\f$ q_w \f$ (default) or \f$ q_n \f$  [\f$\text{kg}/(\text{m}^2 \text{s}\f$])
 */
template<class TypeTag>
void FVSaturation2P<TypeTag>::getFluxOnBoundary(Scalar& update, const Intersection& intersection, CellData& cellDataI)
{
    auto elementI = intersection.inside();

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = elementI.geometry().center();

    // center of face in global coordinates
    const GlobalPosition& globalPosJ = intersection.geometry().center();

    // cell volume, assume linear map here
    Scalar volume = elementI.geometry().volume();
    using std::max;
    Scalar porosity = max(problem_.spatialParams().porosity(elementI), porosityThreshold_);

    if (compressibility_)
    {
        viscosity_[wPhaseIdx] = cellDataI.viscosity(wPhaseIdx);
        viscosity_[nPhaseIdx] = cellDataI.viscosity(nPhaseIdx);
    }

    // local number of faces
    int isIndex = intersection.indexInInside();

    GlobalPosition unitOuterNormal = intersection.centerUnitOuterNormal();
    if (switchNormals_)
        unitOuterNormal *= -1.0;

    Scalar faceArea = intersection.geometry().volume();

    if (velocity().calculateVelocityInTransport())
        velocity().calculateVelocityOnBoundary(intersection, cellDataI);

    //get velocity*normalvector*facearea/(volume*porosity)
    Scalar factorW = (cellDataI.fluxData().velocity(wPhaseIdx, isIndex) * unitOuterNormal) * faceArea;
    Scalar factorNw = (cellDataI.fluxData().velocity(nPhaseIdx, isIndex) * unitOuterNormal) * faceArea;
    Scalar factorTotal = factorW + factorNw;

    // distance vector between barycenters
    GlobalPosition distVec = globalPosJ - globalPosI;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    //get boundary type
    BoundaryTypes bcType;
    problem_.boundaryTypes(bcType, intersection);
    PrimaryVariables boundValues(0.0);

    if (bcType.isDirichlet(satEqIdx))
    {
        problem_.dirichlet(boundValues, intersection);

        Scalar satBound = boundValues[saturationIdx];

        //determine phase saturations from primary saturation variable
        Scalar satWI = cellDataI.saturation(wPhaseIdx);
        Scalar satWBound = 0;
        switch (saturationType_)
        {
        case sw:
        {
            satWBound = satBound;
            break;
        }
        case sn:
        {
            satWBound = 1 - satBound;
            break;
        }
        }

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(elementI.geometry().center());
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), elementI);

        const Scalar pcBound = fluidMatrixInteraction.pc(satWBound);

        Scalar lambdaW = 0;
        Scalar lambdaNw = 0;

        //upwinding of lambda dependend on the phase potential gradients
        if (cellDataI.fluxData().isUpwindCell(wPhaseIdx, isIndex))
        {
            lambdaW = cellDataI.mobility(wPhaseIdx);
            if (compressibility_)
            {
                lambdaW /= cellDataI.density(wPhaseIdx);
            } //divide by density because lambda is saved as lambda*density
        }
        else
        {
            if (compressibility_)
            {
                lambdaW = fluidMatrixInteraction.krw(satWBound)
                        / FluidSystem::viscosity(cellDataI.fluidState(), wPhaseIdx);
            }
            else
            {
                lambdaW = fluidMatrixInteraction.krw(satWBound)
                        / viscosity_[wPhaseIdx];
            }
        }

        if (cellDataI.fluxData().isUpwindCell(nPhaseIdx, isIndex))
        {
            lambdaNw = cellDataI.mobility(nPhaseIdx);
            if (compressibility_)
            {
                lambdaNw /= cellDataI.density(nPhaseIdx);
            } //divide by density because lambda is saved as lambda*density
        }
        else
        {
            if (compressibility_)
            {
                lambdaNw = fluidMatrixInteraction.krn(satWBound)
                        / FluidSystem::viscosity(cellDataI.fluidState(), nPhaseIdx);
            }
            else
            {
                lambdaNw = fluidMatrixInteraction.krn(satWBound)
                        / viscosity_[nPhaseIdx];
            }
        }

        switch (velocityType_)
        {
        case vt:
        {
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], factorTotal, intersection);

            Scalar pcI = cellDataI.capillaryPressure();

            // calculate the saturation gradient
            GlobalPosition pcGradient = unitOuterNormal;
            pcGradient *= (pcI - pcBound) / dist;

            // get the diffusive part -> give 1-sat because sat = S_n and lambda = lambda(S_w) and pc = pc(S_w)
            DimVector flux(0.);
            capillaryFlux().getFlux(flux, intersection, satWI, satWBound, pcGradient);
            Scalar capillaryFlux = flux * unitOuterNormal * faceArea;

            flux = 0.0;
            gravityFlux().getFlux(flux, intersection, satWI, satWBound);
            Scalar gravityFlux = flux * unitOuterNormal * faceArea;

            switch (saturationType_)
            {
            case sw:
            {
                //vt*fw
                factorTotal *= lambdaW / (lambdaW + lambdaNw);
                break;
            }
            case sn:
            {
                //vt*fn
                factorTotal *= lambdaNw / (lambdaW + lambdaNw);
                capillaryFlux *= -1; //add cflFlux for time-stepping
                this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                    viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
                this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                    viscosity_[nPhaseIdx], factorNw, intersection, nPhaseIdx);
                gravityFlux *= -1;
                break;
            }
            }
            //vt*fw
            factorTotal -= capillaryFlux;
            factorTotal += gravityFlux;

            //add cflFlux for time-stepping
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], 10 * capillaryFlux, intersection);
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], 10 * gravityFlux, intersection);

            break;
        }
        default:
        {
            if (compressibility_)
            {
                factorW /= cellDataI.density(wPhaseIdx);
                factorNw /= cellDataI.density(nPhaseIdx);
            }

            //add cflFlux for time-stepping
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], factorNw, intersection, nPhaseIdx);

            break;
        }
        }
    } //end dirichlet boundary

    if (bcType.isNeumann(satEqIdx))
    {
        problem_.neumann(boundValues, intersection);
        factorW = boundValues[wPhaseIdx];
        factorNw = boundValues[nPhaseIdx];
        factorW *= faceArea;
        factorNw *= faceArea;

        //get mobilities
        Scalar lambdaW, lambdaNw;

        lambdaW = cellDataI.mobility(wPhaseIdx);
        lambdaNw = cellDataI.mobility(nPhaseIdx);
        if (compressibility_)
        {
            lambdaW /= cellDataI.density(wPhaseIdx);
            lambdaNw /= cellDataI.density(nPhaseIdx);
            factorW /= cellDataI.density(wPhaseIdx);
            factorNw /= cellDataI.density(nPhaseIdx);
        }
        else
        {
            factorW /= density_[wPhaseIdx];
            factorNw /= density_[nPhaseIdx];
        }

        switch (velocityType_)
        {
        case vt:
        {
            //add cflFlux for time-stepping
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], factorW + factorNw, intersection);
            break;
        }
        default:
        {
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], factorNw, intersection, nPhaseIdx);

            break;
        }
        }

    } //end neumann boundary
    if (bcType.isOutflow(satEqIdx))
    {
        //get mobilities
        Scalar lambdaW = cellDataI.mobility(wPhaseIdx);
        Scalar lambdaNw = cellDataI.mobility(nPhaseIdx);
        if (compressibility_)
        {
            lambdaW /= cellDataI.density(wPhaseIdx);
            lambdaNw /= cellDataI.density(nPhaseIdx);
        }

        if (velocityType_ == vt)
        {
            switch (saturationType_)
            {
            case sw:
            {
                //vt*fw
                factorTotal *= lambdaW / (lambdaW + lambdaNw);
                this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                    viscosity_[nPhaseIdx], factorTotal, intersection);
                break;
            }
            case sn:
            {
                //vt*fn
                factorTotal *= lambdaNw / (lambdaW + lambdaNw);
                this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                    viscosity_[nPhaseIdx], factorTotal, intersection);
                break;
            }
            }
        }
        else
        {
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                                viscosity_[nPhaseIdx], factorNw, intersection, nPhaseIdx);
        }
    }
    switch (velocityType_)
    {
    case vt:
        update -= factorTotal / (volume * porosity); //-:v>0, if flow leaves the cell
        break;
    default:
        switch (saturationType_)
        {
        case sw:
            update -= factorW / (volume * porosity); //-:v>0, if flow leaves the cell
            break;
        case sn:
            update -= factorNw / (volume * porosity); //-:v>0, if flow leaves the cell
            break;
        }
        break;
    }
}

/*!
 * \brief Function which calculates the source update
 *
 *\copydetails FVTransport::getSource(Scalar&,const Element&,CellData&)
 *
 * Source of the fluid phase has to be defined as mass flux (\f$\text{kg}/(\text{m}^3 \text{s}\f$).
 */
template<class TypeTag>
void FVSaturation2P<TypeTag>::getSource(Scalar& update, const Element& element, CellData& cellDataI)
{
    // cell volume, assume linear map here
    Scalar volume = element.geometry().volume();

    using std::max;
    Scalar porosity = max(problem_.spatialParams().porosity(element), porosityThreshold_);

    if (compressibility_)
    {
        viscosity_[wPhaseIdx] = cellDataI.viscosity(wPhaseIdx);
        viscosity_[nPhaseIdx] = cellDataI.viscosity(nPhaseIdx);
    }

    PrimaryVariables sourceVec(0.0);
    problem_.source(sourceVec, element);

    if (compressibility_)
    {
        sourceVec[wPhaseIdx] /= cellDataI.density(wPhaseIdx);
        sourceVec[nPhaseIdx] /= cellDataI.density(nPhaseIdx);
    }
    else
    {
        sourceVec[wPhaseIdx] /= density_[wPhaseIdx];
        sourceVec[nPhaseIdx] /= density_[nPhaseIdx];
    }

    //get mobilities
    Scalar lambdaW = cellDataI.mobility(wPhaseIdx);
    Scalar lambdaNw = cellDataI.mobility(nPhaseIdx);
    if (compressibility_)
    {
        lambdaW /= cellDataI.density(wPhaseIdx);
        lambdaNw /= cellDataI.density(nPhaseIdx);
    }

    switch (saturationType_)
    {
    case sw:
    {
        if (sourceVec[wPhaseIdx] < 0 && cellDataI.saturation(wPhaseIdx) < threshold_)
            sourceVec[wPhaseIdx] = 0.0;

        update += sourceVec[wPhaseIdx] / porosity;
        break;
    }
    case sn:
    {
        if (sourceVec[nPhaseIdx] < 0 && cellDataI.saturation(nPhaseIdx) < threshold_)
            sourceVec[nPhaseIdx] = 0.0;

        update += sourceVec[nPhaseIdx] / porosity;
        break;
    }
    }

    switch (velocityType_)
    {
    case vt:
    {
        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx],
                (sourceVec[wPhaseIdx] + sourceVec[nPhaseIdx]) * -1 * volume, element);
        break;
    }
    default:
    {
        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                            viscosity_[nPhaseIdx], sourceVec[wPhaseIdx] * -1 * volume, element,
                wPhaseIdx);
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNw, viscosity_[wPhaseIdx],
                                            viscosity_[nPhaseIdx], sourceVec[nPhaseIdx] * -1 * volume, element,
                nPhaseIdx);
        break;
    }
    }
}

//! Sets the initial solution \f$ S_0 \f$.
template<class TypeTag>
void FVSaturation2P<TypeTag>::initialize()
{
    ParentType::initialize();

    if (!compressibility_)
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
    }

    // iterate through leaf grid an evaluate c0 at cell center
    for (const auto& element : elements(problem_.gridView()))
    {
        PrimaryVariables initSol(0.0);
        problem_.initial(initSol, element);

        int eIdxGlobal = problem_.variables().index(element);

        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        switch (saturationType_)
        {
        case sw:
        {
                cellData.setSaturation(wPhaseIdx, initSol[saturationIdx]);
                cellData.setSaturation(nPhaseIdx, 1 - initSol[saturationIdx]);
            break;
        }
        case sn:
        {
                cellData.setSaturation(wPhaseIdx,1 -initSol[saturationIdx]);
                cellData.setSaturation(nPhaseIdx, initSol[saturationIdx]);

            break;
        }
        }
    }

    velocity_->initialize();
    capillaryFlux_->initialize();
    gravityFlux_->initialize();

    return;
}

/*!
 * \brief Updates constitutive relations and stores them in the variable class
 *
 * Stores mobility, fractional flow function and capillary pressure for all grid cells.
 */
template<class TypeTag>
void FVSaturation2P<TypeTag>::updateMaterialLaws()
{
    // iterate through leaf grid an evaluate c0 at cell center
    for (const auto& element : elements(problem_.gridView()))
    {
        int eIdxGlobal = problem_.variables().index(element);

        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        //determine phase saturations from primary saturation variable
        Scalar satW = cellData.saturation(wPhaseIdx);

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(elementI.geometry().center());
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

        const Scalar pc = fluidMatrixInteraction.pc(satW);

        cellData.setCapillaryPressure(pc);

        // initialize mobilities
        const Scalar mobilityW = fluidMatrixInteraction.krw(satW) / viscosity_[wPhaseIdx];
        const Scalar mobilityNw = fluidMatrixInteraction.krn(satW) / viscosity_[nPhaseIdx];

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNw);

        //initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNw));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNw / (mobilityW + mobilityNw));
    }
    return;
}

} // end namespace Dumux
#endif
