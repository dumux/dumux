// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_FVSATURATION2P_HH
#define DUMUX_FVSATURATION2P_HH

#include <dumux/decoupled/2p/transport/transportproperties2p.hh>
#include "dumux/decoupled/common/fv/fvtransport.hh"

/**
 * @file
 * @brief  Finite Volume discretization of a saturation transport equation
 * @author Markus Wolff
 */

namespace Dumux
{
//! \ingroup FVSaturation2p
//! \brief The finite volume discretization of a saturation transport equation
/*! Provides a Finite Volume implementation for the evaluation
 *  of equations of the form
 *
 *  \f[
 *  \phi \frac{\partial (\rho_\alpha S_\alpha)}{\partial t} + \text{div}\, (\rho_\alpha \boldsymbol{v_\alpha}) = q_\alpha,
 *  \f]
 *
 *  where \f$ S_\alpha \f$ is the saturation of phase alpha (wetting (\f$ w \f$), non-wetting (\f$ n \f$)) and \f$ \boldsymbol v_\alpha \f$ is the phase velocity defined by the multi-phase Darcy equation.
 *  If a phase velocity is reconstructed from the pressure solution it can be directly inserted in the previous equation. In the incompressible case the equation is further divided by the phase density \f$ \rho_\alpha \f$. If a total velocity is reconstructed the saturation equation is reformulated into:
 *
 * \f[
 *  \phi \frac{\partial S_w}{\partial t} + f_w \text{div}\, \boldsymbol{v}_{t} + f_w \lambda_n \boldsymbol{K}\left(\text{grad}\, p_c + (\rho_n-\rho_w) g \text{grad} z \right)= q_\alpha,
 * \f]
 * to get a wetting phase saturation or
 * \f[
 * \phi \frac{\partial S_n}{\partial t} + f_n \text{div}\, \boldsymbol{v}_{t} - f_n \lambda_w \boldsymbol{K}\left(\text{grad}\, p_c + (\rho_n-\rho_w) g \text{grad} z \right)= q_\alpha,
 * \f]
 * if the non-wetting phase saturation is the primary transport variable.
 *
 *  The total velocity formulation is only implemented for incompressible fluids and \f$ f_\alpha \f$ is the fractional flow function, \f$ \lambda_\alpha \f$ is the mobility, \f$ \boldsymbol K \f$ the absolute permeability,
 *  \f$ p_c \f$ the capillary pressure, \f$ \rho \f$ the fluid density, \f$ g \f$ the gravity constant, and \f$ q \f$ the source term.
 *
 *
 *  In the IMPES models the default setting is:
 *
 *      - formulation: \f$ p_w-S_w \f$ (Property: <tt>Formulation</tt> defined as <tt>DecoupledTwoPCommonIndices::pwSw</tt>)
 *      - compressibility: disabled (Property: <tt>EnableCompressibility</tt> set to <tt>false</tt>)
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVSaturation2P: public FVTransport<TypeTag>
{
    typedef FVTransport<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Velocity) Velocity;
    typedef typename GET_PROP_TYPE(TypeTag, CapillaryFlux) CapillaryFlux;
    typedef typename GET_PROP_TYPE(TypeTag, GravityFlux) GravityFlux;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressEqIdx = Indices::pressEqIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

    Velocity& velocity()
    {
        return *velocity_;
    }

    Velocity& velocity() const
    {
        return *velocity_;
    }

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
    // Function which calculates the flux update
    void getFlux(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    // Function which calculates the boundary flux update
    void getFluxOnBoundary(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    // Function which calculates the source update
    void getSource(Scalar& update, const Element& element, CellData& cellDataI);

    // Sets the initial solution
    void initialize();

    // Update the values of the material laws and constitutive relations.
    void updateMaterialLaws();


    /* \brief Writes the current values of the primary transport variable into the <tt>transportedQuantity</tt>-vector (comes as function argument)
     *
     * \copydetails FVTransport::getTransportedQuantity(TransportSolutionType&)
     *
     */
    void getTransportedQuantity(TransportSolutionType& transportedQuantity)
    {
        int size = problem_.gridView().size(0);
        transportedQuantity.resize(size);
        for (int i = 0; i < size; i++)
        {
            switch (saturationType_)
            {
            case Sw:
            {
                transportedQuantity[i] = problem_.variables().cellData(i).saturation(wPhaseIdx);
                break;
            }
            case Sn:
            {
                transportedQuantity[i] = problem_.variables().cellData(i).saturation(nPhaseIdx);
                break;
            }
            }
        }
    }

    /*! \brief Updates the primary transport variable.
     *
     * \copydetails FVTransport::updateTransportedQuantity(TransportSolutionType&)
     */
    void updateTransportedQuantity(TransportSolutionType& updateVec)
    {
        updateSaturationSolution(updateVec);

//        std::cout<<"update = "<<updateVec<<"\n";
    }

    /*! \brief Globally updates the saturation solution
     *
     * \param updateVec Vector containing the global update.
     */
    void updateSaturationSolution(TransportSolutionType& updateVec)
    {
        Scalar dt = problem_.timeManager().timeStepSize();
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            updateSaturationSolution(i, updateVec[i][0], dt);
        }
    }

    /*! \brief Updates the saturation solution of a cell
     *
     * Calculates secondary saturation variables and stores saturations.
     *
     * \param globalIdx Global cell index
     * \param update Cell saturation update
     * \param dt Current time step
     */
    void updateSaturationSolution(int globalIdx, Scalar update, Scalar dt)
    {
        CellData& cellData = problem_.variables().cellData(globalIdx);

        switch (saturationType_)
        {
        case Sw:
        {
            Scalar sat = cellData.saturation(wPhaseIdx) + dt*update;

                cellData.setSaturation(wPhaseIdx, sat);
                cellData.setSaturation(nPhaseIdx, 1 - sat);
            break;
        }
        case Sn:
        {
            Scalar sat = cellData.saturation(nPhaseIdx) + dt*update;

                cellData.setSaturation(wPhaseIdx,1 -sat);
                cellData.setSaturation(nPhaseIdx, sat);
            break;
        }
        }
    }

    /*! \brief Adds saturation output to the output file
     *
     * Adds the phase saturation to the output. If the velocity is calculated in the transport model it is also added to the output.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        TransportSolutionType *saturationW = writer.allocateManagedBuffer(size);
        TransportSolutionType *saturationN = writer.allocateManagedBuffer(size);

        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            (*saturationW)[i] = cellData.saturation(wPhaseIdx);
            (*saturationN)[i] = cellData.saturation(nPhaseIdx);
        }

        writer.attachCellData(*saturationW, "wetting saturation");
        writer.attachCellData(*saturationN, "nonwetting saturation");

        if (velocity().calculateVelocityInTransport())
            velocity().addOutputVtkFields(writer);

        return;
    }

    /*! \brief  Function for serialization of the primary transport variable.
     *
     *\copydetails FVTransport::serializeEntity(std::ostream&,const Element&)
     *
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int globalIdx = problem_.variables().index(element);

        Scalar sat = 0.0;
        switch (saturationType_)
        {
        case Sw:
            sat = problem_.variables().cellData(globalIdx).saturation(wPhaseIdx);
        break;
        case Sn:
            sat = problem_.variables().cellData(globalIdx).saturation(nPhaseIdx);
            break;
        }

        outstream << sat;
    }

    /*! \brief  Function for deserialization of the primary transport variable.
     *
     *\copydetails FVTransport::deserializeEntity(std::istream&,const Element&)
     *
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int globalIdx = problem_.variables().index(element);

        Scalar sat = 0.;
        instream >> sat;

        switch (saturationType_)
        {
        case Sw:
            problem_.variables().cellData(globalIdx).setSaturation(wPhaseIdx, sat);
            problem_.variables().cellData(globalIdx).setSaturation(nPhaseIdx, 1-sat);
        break;
        case Sn:
            problem_.variables().cellData(globalIdx).setSaturation(nPhaseIdx, sat);
            problem_.variables().cellData(globalIdx).setSaturation(wPhaseIdx, 1-sat);
            break;
        }
    }


    /*! \brief Constructs a FVSaturation2P object
     *
     * \param problem A problem class object
     */
    FVSaturation2P(Problem& problem) :
            ParentType(problem), problem_(problem), threshold_(1e-6), switchNormals_(GET_PARAM(TypeTag, bool, SwitchNormals))
    {
        if (compressibility_ && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented,
                    "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (saturationType_ != Sw && saturationType_ != Sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (pressureType_ != pw && pressureType_ != pn && pressureType_ != pglobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }

        capillaryFlux_ = new CapillaryFlux(problem);
        gravityFlux_ = new GravityFlux(problem);
        velocity_ = new Velocity(problem);

        if (!compressibility_)
        {
            const Element& element = *(problem_.gridView().template begin<0> ());
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
    }

    //! Destructor
    ~FVSaturation2P()
    {
        delete capillaryFlux_;
        delete gravityFlux_;
        delete velocity_;
    }

private:
    Problem& problem_;
    Velocity* velocity_;
    CapillaryFlux* capillaryFlux_;
    GravityFlux* gravityFlux_;

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation);
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    const Scalar threshold_;
    const bool switchNormals_;
};

/*! \brief Function which calculates the flux update
 *
 * \copydetails FVTransport::getFlux(Scalar&,const Intersection&,CellData&)
 *
 * If a total velocity formulation is used this functions calculates not only the advective flux but also fluxes due to gravity and capillary diffusion.
 * These have to be defined separately as implementation of a DiffusivePart or ConvectivePart (e.g. GravityPart / CapillaryDiffusion ) and added to the property system via properties <tt>CapillaryFlux</tt> and <tt>GravityFlux</tt>.
 */
template<class TypeTag>
void FVSaturation2P<TypeTag>::getFlux(Scalar& update, const Intersection& intersection, CellData& cellDataI)
{
    ElementPointer elementI = intersection.inside();
    ElementPointer elementJ = intersection.outside();

    const CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(*elementJ));

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = elementI->geometry().center();
    const GlobalPosition& globalPosJ = elementJ->geometry().center();

    // cell volume, assume linear map here
    Scalar volume = elementI->geometry().volume();
    Scalar porosity = problem_.spatialParameters().porosity(*elementI);

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
    Scalar factorNW = (cellDataI.fluxData().velocity(nPhaseIdx, isIndex) * unitOuterNormal) * faceArea;
    Scalar factorTotal = factorW + factorNW;

    // distance vector between barycenters
    GlobalPosition distVec = globalPosJ - globalPosI;
    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    bool takeNeighbor = (elementI->level() < elementJ->level());
    //get phase potentials
    bool upwindWI =
            (takeNeighbor) ? !cellDataJ.fluxData().isUpwindCell(wPhaseIdx, intersection.indexInOutside()) :
                    cellDataI.fluxData().isUpwindCell(wPhaseIdx, isIndex);
    bool upwindNWI =
            (takeNeighbor) ? !cellDataJ.fluxData().isUpwindCell(nPhaseIdx, intersection.indexInOutside()) :
                    cellDataI.fluxData().isUpwindCell(nPhaseIdx, isIndex);

    Scalar lambdaW = 0;
    Scalar lambdaNW = 0;

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

    if (upwindNWI)
    {
        lambdaNW = cellDataI.mobility(nPhaseIdx);
        if (compressibility_)
        {
            lambdaNW /= cellDataI.density(nPhaseIdx);
        } //divide by density because lambda is saved as lambda*density
    }
    else
    {
        lambdaNW = cellDataJ.mobility(nPhaseIdx);
        if (compressibility_)
        {
            lambdaNW /= cellDataJ.density(nPhaseIdx);
        } //divide by density because lambda is saved as lambda*density
    }

    switch (velocityType_)
    {
    case vt:
    {
        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorTotal, intersection);

        //determine phase saturations from primary saturation variable
        Scalar satWI = cellDataI.saturation(wPhaseIdx);
        Scalar satWJ = cellDataJ.saturation(wPhaseIdx);

        Scalar pcI = cellDataI.capillaryPressure();
        Scalar pcJ = cellDataJ.capillaryPressure();

        // calculate the saturation gradient
        GlobalPosition pcGradient = unitOuterNormal;
        pcGradient *= (pcI - pcJ) / dist;

        // get the diffusive part
        FieldVector flux(0.);
        capillaryFlux().getFlux(flux, intersection, satWI, satWJ, pcGradient);
        Scalar capillaryFlux = (flux * unitOuterNormal * faceArea);

        flux = 0.0;
        gravityFlux().getFlux(flux, intersection, satWI, satWJ);
        Scalar gravityFlux =  (flux * unitOuterNormal * faceArea);

        switch (saturationType_)
        {
        case Sw:
        {
            //vt*fw
            factorTotal *= lambdaW / (lambdaW + lambdaNW);
            break;
        }
        case Sn:
        {
            //vt*fn
            factorTotal *= lambdaNW / (lambdaW + lambdaNW);
            capillaryFlux *= -1;
            gravityFlux *= -1;
            break;
        }
        }
        factorTotal -= capillaryFlux;
        factorTotal += gravityFlux;

        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], 10 * capillaryFlux, intersection);
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], 10 * gravityFlux, intersection);

        break;
    }
    default:
    {
        if (compressibility_)
        {
            factorW /= cellDataI.density(wPhaseIdx);
            factorNW /= cellDataI.density(nPhaseIdx);
        }

        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorNW, intersection, nPhaseIdx);

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
        case Sw:
            update -= factorW / (volume * porosity); //-:v>0, if flow leaves the cell
            break;
        case Sn:
            update -= factorNW / (volume * porosity); //-:v>0, if flow leaves the cell
            break;
        }
        break;
    }
}

/*! \brief Function which calculates the boundary flux update
 *
 * \copydetails FVTransport::getFluxOnBoundary(Scalar&,const Intersection&,CellData&)
 *
 * Dirichlet boundary condition is a phase saturation depending on the formulation (\f$ S_w \f$ (default) or \f$ S_n \f$),
 * Neumann boundary condition are phase mass fluxes (\f$ q_w \f$ (default) or \f$ q_n \f$  [\f$\text{kg}/(\text{m}^2 \text{s}\f$])
 */
template<class TypeTag>
void FVSaturation2P<TypeTag>::getFluxOnBoundary(Scalar& update, const Intersection& intersection, CellData& cellDataI)
{
    ElementPointer elementI = intersection.inside();

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = elementI->geometry().center();

    // center of face in global coordinates
    const GlobalPosition& globalPosJ = intersection.geometry().center();

    // cell volume, assume linear map here
    Scalar volume = elementI->geometry().volume();
    Scalar porosity = problem_.spatialParameters().porosity(*elementI);

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
    Scalar factorNW = (cellDataI.fluxData().velocity(nPhaseIdx, isIndex) * unitOuterNormal) * faceArea;
    Scalar factorTotal = factorW + factorNW;

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
        case Sw:
        {
            satWBound = satBound;
            break;
        }
        case Sn:
        {
            satWBound = 1 - satBound;
            break;
        }
        }

        Scalar pcBound = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*elementI), satWBound);

        Scalar lambdaW = 0;
        Scalar lambdaNW = 0;

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
                lambdaW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*elementI), satWBound)
                        / FluidSystem::viscosity(cellDataI.fluidState(), wPhaseIdx);
            }
            else
            {
                lambdaW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*elementI), satWBound)
                        / viscosity_[wPhaseIdx];
            }
        }

        if (cellDataI.fluxData().isUpwindCell(nPhaseIdx, isIndex))
        {
            lambdaNW = cellDataI.mobility(nPhaseIdx);
            if (compressibility_)
            {
                lambdaNW /= cellDataI.density(nPhaseIdx);
            } //divide by density because lambda is saved as lambda*density
        }
        else
        {
            if (compressibility_)
            {
                lambdaNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*elementI), satWBound)
                        / FluidSystem::viscosity(cellDataI.fluidState(), nPhaseIdx);
            }
            else
            {
                lambdaNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*elementI), satWBound)
                        / viscosity_[nPhaseIdx];
            }
        }

        switch (velocityType_)
        {
        case vt:
        {
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorTotal, intersection);

            Scalar pcI = cellDataI.capillaryPressure();

            // calculate the saturation gradient
            GlobalPosition pcGradient = unitOuterNormal;
            pcGradient *= (pcI - pcBound) / dist;

            // get the diffusive part -> give 1-sat because sat = S_n and lambda = lambda(S_w) and pc = pc(S_w)
            FieldVector flux(0.);
            capillaryFlux().getFlux(flux, intersection, satWI, satWBound, pcGradient);
            Scalar capillaryFlux = flux * unitOuterNormal * faceArea;

            flux = 0.0;
            gravityFlux().getFlux(flux, intersection, satWI, satWBound);
            Scalar gravityFlux = flux * unitOuterNormal * faceArea;

            switch (saturationType_)
            {
            case Sw:
            {
                //vt*fw
                factorTotal *= lambdaW / (lambdaW + lambdaNW);
                break;
            }
            case Sn:
            {
                //vt*fn
                factorTotal *= lambdaNW / (lambdaW + lambdaNW);
                capillaryFlux *= -1; //add cflFlux for time-stepping
                this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
                this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorNW, intersection, nPhaseIdx);
                gravityFlux *= -1;
                break;
            }
            }
            //vt*fw
            factorTotal -= capillaryFlux;
            factorTotal += gravityFlux;

            //add cflFlux for time-stepping
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], 10 * capillaryFlux, intersection);
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], 10 * gravityFlux, intersection);

            break;
        }
        default:
        {
            if (compressibility_)
            {
                factorW /= cellDataI.density(wPhaseIdx);
                factorNW /= cellDataI.density(nPhaseIdx);
            }

            //add cflFlux for time-stepping
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorNW, intersection, nPhaseIdx);

            break;
        }
        }
    } //end dirichlet boundary

    if (bcType.isNeumann(satEqIdx))
    {
        problem_.neumann(boundValues, intersection);
        factorW = boundValues[wPhaseIdx];
        factorNW = boundValues[nPhaseIdx];
        factorW *= faceArea;
        factorNW *= faceArea;

        //get mobilities
        Scalar lambdaW, lambdaNW;

        lambdaW = cellDataI.mobility(wPhaseIdx);
        lambdaNW = cellDataI.mobility(nPhaseIdx);
        if (compressibility_)
        {
            lambdaW /= cellDataI.density(wPhaseIdx);
            lambdaNW /= cellDataI.density(nPhaseIdx);
            factorW /= cellDataI.density(wPhaseIdx);
            factorNW /= cellDataI.density(nPhaseIdx);
        }
        else
        {
            factorW /= density_[wPhaseIdx];
            factorNW /= density_[nPhaseIdx];
        }

        switch (velocityType_)
        {
        case vt:
        {
            //add cflFlux for time-stepping
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorW + factorNW, intersection);
            break;
        }
        default:
        {
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorNW, intersection, nPhaseIdx);

            break;
        }
        }

    } //end neumann boundary
    if (bcType.isOutflow(satEqIdx))
    {
        //get mobilities
        Scalar lambdaW = cellDataI.mobility(wPhaseIdx);
        Scalar lambdaNW = cellDataI.mobility(nPhaseIdx);
        if (compressibility_)
        {
            lambdaW /= cellDataI.density(wPhaseIdx);
            lambdaNW /= cellDataI.density(nPhaseIdx);
        }

        if (velocityType_ == vt)
        {
            switch (saturationType_)
            {
            case Sw:
            {
                //vt*fw
                factorTotal *= lambdaW / (lambdaW + lambdaNW);
                this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorTotal, intersection);
                break;
            }
            case Sn:
            {
                //vt*fn
                factorTotal *= lambdaNW / (lambdaW + lambdaNW);
                this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorTotal, intersection);
                break;
            }
            }
        }
        else
        {
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorW, intersection, wPhaseIdx);
            this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], factorNW, intersection, nPhaseIdx);
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
        case Sw:
            update -= factorW / (volume * porosity); //-:v>0, if flow leaves the cell
            break;
        case Sn:
            update -= factorNW / (volume * porosity); //-:v>0, if flow leaves the cell
            break;
        }
        break;
    }
}

/*! \brief Function which calculates the source update
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

    Scalar porosity = problem_.spatialParameters().porosity(element);

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
    Scalar lambdaNW = cellDataI.mobility(nPhaseIdx);
    if (compressibility_)
    {
        lambdaW /= cellDataI.density(wPhaseIdx);
        lambdaNW /= cellDataI.density(nPhaseIdx);
    }

    switch (saturationType_)
    {
    case Sw:
    {
        if (sourceVec[wPhaseIdx] < 0 && cellDataI.saturation(wPhaseIdx) < threshold_)
            sourceVec[wPhaseIdx] = 0.0;

        update += sourceVec[wPhaseIdx] / porosity;
        break;
    }
    case Sn:
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
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx],
                (sourceVec[wPhaseIdx] + sourceVec[nPhaseIdx]) * volume, element);
        break;
    }
    default:
    {
        //add cflFlux for time-stepping
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], sourceVec[wPhaseIdx] * volume, element,
                wPhaseIdx);
        this->evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosity_[wPhaseIdx], viscosity_[nPhaseIdx], sourceVec[nPhaseIdx] * volume, element,
                nPhaseIdx);
        break;
    }
    }
}

//! Sets the initial solution \f$ S_0 \f$.
template<class TypeTag>
void FVSaturation2P<TypeTag>::initialize()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        PrimaryVariables initSol(0.0);
        problem_.initial(initSol, *eIt);

        int globalIdx = problem_.variables().index(*eIt);

        CellData& cellData = problem_.variables().cellData(globalIdx);

        switch (saturationType_)
        {
        case Sw:
        {
                cellData.setSaturation(wPhaseIdx, initSol[saturationIdx]);
                cellData.setSaturation(nPhaseIdx, 1 - initSol[saturationIdx]);
            break;
        }
        case Sn:
        {
                cellData.setSaturation(wPhaseIdx,1 -initSol[saturationIdx]);
                cellData.setSaturation(nPhaseIdx, initSol[saturationIdx]);

            break;
        }
        }
    }

    return;
}

/*! \brief Updates constitutive relations and stores them in the variable class
 *
 * Stores mobility, fractional flow function and capillary pressure for all grid cells.
 *
 */
template<class TypeTag>
void FVSaturation2P<TypeTag>::updateMaterialLaws()
{
    ElementIterator eItBegin = problem_.gridView().template begin<0>();
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem_.variables().index(*eIt);

        CellData& cellData = problem_.variables().cellData(globalIdx);

        //determine phase saturations from primary saturation variable
        Scalar satW = cellData.saturation(wPhaseIdx);
        Scalar satNW = cellData.saturation(nPhaseIdx);

        Scalar pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

        cellData.setSaturation(wPhaseIdx, satW);
        cellData.setSaturation(nPhaseIdx, satNW);

        cellData.setCapillaryPressure(pc);

        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosity_[wPhaseIdx];
        Scalar mobilityNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosity_[nPhaseIdx];

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNW);

        //initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNW));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNW / (mobilityW + mobilityNW));
    }
    return;
}

}
#endif
