// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
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

#include "dumux/decoupled/2p/transport/fv/diffusivepart.hh"
#include "dumux/decoupled/2p/transport/fv/convectivepart.hh"
#include <dumux/decoupled/2p/transport/transportproperties.hh>
#include <dumux/decoupled/2p/2pproperties.hh>
#include "evalcflflux_default.hh"

/**
 * @file
 * @brief  Finite Volume discretization of a saturation transport equation
 * @author Markus Wolff
 */

namespace Dumux
{
//! \ingroup Saturation2p
//! \brief The finite volume discretization of a saturation transport equation
/*! Provides a Finite Volume implementation for the evaluation
 *  of equations of the form
 *
 *  \f[\frac{\partial S_\alpha}{\partial t} + \text{div}\, \boldsymbol{v_\alpha} = q_\alpha,\f]
 *
 *  where \f$S_\alpha\f$ is the saturation of phase alpha (wetting (w), non-wetting (n)) and \f$\boldsymbol{v}_\alpha\f$ is the phase velocity calculated by the multi-phase Darcy equation,
 *  and of the form
 *
 * \f[\frac{\partial S_w}{\partial t} + f_w \text{div}\, \boldsymbol{v}_{t} + f_w \lambda_n \boldsymbol{K}\left(\text{grad}\, p_c + (\rho_n-\rho_w) g \text{grad} z \right)= q_\alpha,\f]
 *
 * \f[\frac{\partial S_n}{\partial t} + f_n \text{div}\, \boldsymbol{v}_{t} - f_n \lambda_w \boldsymbol{K}\left(\text{grad}\, p_c + (\rho_n-\rho_w) g \text{grad} z \right)= q_\alpha,\f]
 *
 *  where \f$f_\alpha\f$ is the fractional flow function, \f$\lambda_\alpha\f$ is the mobility, \f$\boldsymbol{K}\f$ the absolute permeability,
 *  \f$p_c\f$ the capillary pressure, \f$\rho\f$ the fluid density, \f$g\f$ the gravity constant, and \f$q\f$ the source term.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVSaturation2P
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElementContainer;
    typedef Dune::GenericReferenceElements<Scalar, dim - 1> ReferenceElementFaceContainer;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DiffusivePart)) DiffusivePart;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ConvectivePart)) ConvectivePart;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(EvalCflFluxFunction)) EvalCflFluxFunction;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Indices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NonwettingPhase)) NonwettingPhase;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressEqIdx = Indices::pressEqIdx,
        satEqIdx = Indices::satEqIdx
    };

    typedef typename SolutionTypes::ScalarSolution RepresentationType;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    DiffusivePart& diffusivePart()
    {
        return *diffusivePart_;
    }

    const DiffusivePart& diffusivePart() const
    {
        return *diffusivePart_;
    }

    ConvectivePart& convectivePart()
    {
        return *convectivePart_;
    }

    const ConvectivePart& convectivePart() const
    {
        return *convectivePart_;
    }

    EvalCflFluxFunction& evalCflFluxFunction()
    {
        return *evalCflFluxFunction_;
    }

    const EvalCflFluxFunction& evalCflFluxFunction() const
    {
        return *evalCflFluxFunction_;
    }

public:
    //! Calculate the update vector.
    /*!
     *  \param[in]  t         time
     *  \param[in] dt         time step size
     *  \param[in] updateVec  vector for the update values
     *  \param[in] impes      variable is true if an impes algorithm is used and false if the transport part is solved independently
     *
     *  This method calculates the update vector \f$\boldsymbol{u}\f$ of the discretized equation
     *
     *  \f[S_{new} = S_{old} + u,\f]
     *
     *  where \f$ u = \sum_{element faces} \boldsymbol{F}\f$ and \f$\boldsymbol{F}\f$ is the flux over an element face.
     *
     *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
     *  employing a CFL condition.
     */
    void update(const Scalar t, Scalar& dt, RepresentationType& updateVec, bool impes);

    //! Sets the initial solution \f$S_0\f$.
    void initialize();

    //! Update the values of the material laws and constitutive relations.
    /*!
     *  Constitutive relations like capillary pressure-saturation relationships, mobility-saturation relationships... are updated and stored in the variable class
     *  of type Dumux::VariableClass2P. The update has to be done when new saturation are available.
     */
    void updateMaterialLaws(RepresentationType& saturation, bool iterate);

    void updateSaturationSolution()
    {
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            Scalar sat = problem_.variables().primaryVariablesGlobal(satEqIdx)[i];
            updateSaturationSolution(i, sat);
        }
    }

    void updateSaturationSolution(int globalIdx, Scalar sat)
    {
        CellData& cellData = problem_.variables().cellData(globalIdx);
        switch (saturationType_)
        {
        case Sw:
        {
            if (compressibility_)
            {
                cellData.fluidState().setSaturation(wPhaseIdx, sat);
                cellData.fluidState().setSaturation(nPhaseIdx, 1 - sat);
            }
            else
            {
                cellData.setSaturation(wPhaseIdx, sat);
                cellData.setSaturation(nPhaseIdx, 1 - sat);
            }
            break;
        }
        case Sn:
        {
            if (compressibility_)
            {
                cellData.fluidState().setSaturation(wPhaseIdx, 1 - sat);
                cellData.fluidState().setSaturation(nPhaseIdx, sat);
            }
            else
            {
                cellData.setSaturation(wPhaseIdx,1 -sat);
                cellData.setSaturation(nPhaseIdx, sat);
            }
            break;
        }
        }
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        typename Variables::ScalarSolutionType *saturationW = writer.allocateManagedBuffer(size);
        typename Variables::ScalarSolutionType *saturationN = writer.allocateManagedBuffer(size);

        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            (*saturationW)[i] = cellData.saturation(wPhaseIdx);
            (*saturationN)[i] = cellData.saturation(nPhaseIdx);
        }

        writer.attachCellData(*saturationW, "wetting saturation");
        writer.attachCellData(*saturationN, "nonwetting saturation");

        return;
    }

    void indicatorSaturation(RepresentationType &indicator, Scalar &globalMin, Scalar &globalMax);

    // serialization methods
    //! Function needed for restart option.
    template<class Restarter>
    void serialize(Restarter &res)
    {
        problem_.variables().serialize<Restarter>(res);
    }

    //! Function needed for restart option.
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        problem_.variables().deserialize<Restarter>(res);
    }

    //! Constructs a FVSaturation2P object
    /**

     * \param problem a problem class object
     */

    FVSaturation2P(Problem& problem) :
            problem_(problem), switchNormals_(false), threshold_(1e-6)
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

        diffusivePart_ = new DiffusivePart(problem);
        convectivePart_ = new ConvectivePart(problem);
        evalCflFluxFunction_ = new EvalCflFluxFunction(problem);
    }

    ~FVSaturation2P()
    {
        delete diffusivePart_;
        delete convectivePart_;
        delete evalCflFluxFunction_;
    }

private:
    Problem& problem_;
    DiffusivePart* diffusivePart_;
    ConvectivePart* convectivePart_;
    EvalCflFluxFunction* evalCflFluxFunction_;

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility));
    ;
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation));
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation));
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation));
    bool switchNormals_;

    const Scalar threshold_;
};

/*!
 * \name Methods for adaptive refinement / coarsening of grids
 */
// \{
/*!
 * Indicator for grid refinement: saturation gradient
 */
template<class TypeTag>
void FVSaturation2P<TypeTag>::indicatorSaturation(RepresentationType &indicator, Scalar &globalMin, Scalar &globalMax)
{
    // 1) calculate Indicator -> min, maxvalues
    // Schleife über alle Leaf-Elemente
    for (ElementIterator it = problem_.gridView().template begin<0>(); it != problem_.gridView().template end<0>();
            ++it)
    {
        // Bestimme maximale und minimale Sättigung
        // Index des aktuellen Leaf-Elements
        int indexi = problem_.variables().index(*it);
        globalMin = std::min(problem_.variables().saturation()[indexi][0], globalMin);
        globalMax = std::max(problem_.variables().saturation()[indexi][0], globalMax);

        // Berechne Verfeinerungsindikator an allen Zellen
        IntersectionIterator isend = problem_.gridView().iend(*it);
        for (IntersectionIterator is = problem_.gridView().ibegin(*it); is != isend; ++is)
        {
            const typename IntersectionIterator::Intersection &intersection = *is;
            // Steige aus, falls es sich nicht um einen Nachbarn handelt
            if (!intersection.neighbor())
                continue;

            // Greife auf Nachbarn zu
            const Element &outside = *intersection.outside();
            int indexj = problem_.variables().index(outside);

            // Jede Intersection nur von einer Seite betrachten
            if (it.level() > outside.level() || (it.level() == outside.level() && indexi < indexj))
            {

                Scalar localdelta = std::abs(
                        problem_.variables().saturation()[indexi][0] - problem_.variables().saturation()[indexj][0]);
                indicator[indexi][0] = std::max(indicator[indexi][0], localdelta);
                indicator[indexj][0] = std::max(indicator[indexj][0], localdelta);
            }
        }
    }
}

template<class TypeTag>
void FVSaturation2P<TypeTag>::update(const Scalar t, Scalar& dt, RepresentationType& updateVec, bool impes = false)
{
    if (!impes)
    {
        updateMaterialLaws();
    }

    // initialize dt very large
    dt = 1E100;

    // set update vector to zero
    updateVec = 0;

    // some phase properties
    //const GlobalPosition& gravity = problem_.gravity();

    BoundaryTypes bcType;

    ElementIterator eItBegin = problem_.gridView().template begin<0>();
    Scalar densityW = 0;
    Scalar densityNW = 0;
    Scalar viscosityW = 0;
    Scalar viscosityNW = 0;

    if (!compressibility_)
    {
        Scalar temp = problem_.temperature(*eItBegin);
        Scalar pRef = problem_.referencePressure(*eItBegin);
        densityW = WettingPhase::density(temp, pRef);
        densityNW = NonwettingPhase::density(temp, pRef);
        viscosityW = WettingPhase::viscosity(temp, pRef);
        viscosityNW = NonwettingPhase::viscosity(temp, pRef);
    }

    // compute update vector
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
#if HAVE_MPI
        if (eIt->partitionType() != Dune::InteriorEntity)
        {
            continue;
        }
#endif
        //
        const GlobalPosition& globalPos = eIt->geometry().center();

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().volume();

        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        CellData& cellDataI = problem_.variables().cellData(globalIdxI);

        //for benchmark only!
        //        problem_.variables().storeSrn(residualSatNW, globalIdxI);
        //for benchmark only!

        Scalar porosity = problem_.spatialParameters().porosity(*eIt);

        Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
        Scalar lambdaNWI = cellDataI.mobility(nPhaseIdx);

        if (compressibility_)
        {
            viscosityW = cellDataI.viscosity(wPhaseIdx);
            viscosityNW = cellDataI.viscosity(nPhaseIdx);
        }

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            // local number of faces
            int isIndex = isIt->indexInInside();

            GlobalPosition unitOuterNormal = isIt->centerUnitOuterNormal();
            if (switchNormals_)
                unitOuterNormal *= -1.0;

            Scalar faceArea = isIt->geometry().volume();

            //get velocity*normalvector*facearea/(volume*porosity)
            Scalar factorW = (cellDataI.fluxData().velocity(wPhaseIdx, isIndex) * unitOuterNormal) * faceArea;
            Scalar factorNW = (cellDataI.fluxData().velocity(nPhaseIdx, isIndex) * unitOuterNormal) * faceArea;
            Scalar factorTotal = factorW + factorNW;

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem_.variables().index(*neighborPointer);

                CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                // distance vector between barycenters
                GlobalPosition distVec = globalPosNeighbor - globalPos;
                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                bool takeNeighbor = (eIt->level() < neighborPointer->level());
                //get phase potentials
                bool upwindWI =
                        (takeNeighbor) ? !cellDataJ.fluxData().isUpwindCell(wPhaseIdx, isIt->indexInOutside()) :
                                cellDataI.fluxData().isUpwindCell(wPhaseIdx, isIndex);
                bool upwindNWI =
                        (takeNeighbor) ? !cellDataJ.fluxData().isUpwindCell(nPhaseIdx, isIt->indexInOutside()) :
                                cellDataI.fluxData().isUpwindCell(nPhaseIdx, isIndex);

                Scalar lambdaW = 0;
                Scalar lambdaNW = 0;

                //upwinding of lambda dependend on the phase potential gradients
                if (upwindWI)
                {
                    lambdaW = lambdaWI;
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
                    lambdaNW = lambdaNWI;
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
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorTotal, *isIt);

                    //determine phase saturations from primary saturation variable
                    Scalar satWI = cellDataI.saturation(wPhaseIdx);
                    Scalar satWJ = cellDataJ.saturation(wPhaseIdx);

                    Scalar pcI = cellDataI.capillaryPressure();
                    Scalar pcJ = cellDataJ.capillaryPressure();

                    // calculate the saturation gradient
                    GlobalPosition pcGradient = unitOuterNormal;
                    pcGradient *= (pcI - pcJ) / dist;

                    // get the diffusive part
                    Scalar diffPart = diffusivePart()(*eIt, isIndex, satWI, satWJ, pcGradient) * unitOuterNormal
                            * faceArea;

                    Scalar convPart = convectivePart()(*eIt, isIndex, satWI, satWJ) * unitOuterNormal * faceArea;

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
                        diffPart *= -1;
                        convPart *= -1;
                        break;
                    }
                    }
                    factorTotal -= diffPart;
                    factorTotal += convPart;

                    //add cflFlux for time-stepping
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, 10 * diffPart, *isIt);
                    evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, 10 * convPart, *isIt);

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
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorW, *isIt,
                                wPhaseIdx);
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorNW, *isIt,
                                nPhaseIdx);

                    break;
                }
                }
            } //end intersection with neighbor element

            // handle boundary face
            if (isIt->boundary())
            {
                //get boundary type
                problem_.boundaryTypes(bcType, *isIt);
                PrimaryVariables boundValues(0.0);

                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().center();

                // distance vector between barycenters
                GlobalPosition distVec = globalPosFace - globalPos;
                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                if (bcType.isDirichlet(satEqIdx))
                {
                    problem_.dirichlet(boundValues, *isIt);

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

                    Scalar pcBound = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satBound);

                    Scalar lambdaW = 0;
                    Scalar lambdaNW = 0;

                    //upwinding of lambda dependend on the phase potential gradients
                    if (cellDataI.fluxData().isUpwindCell(wPhaseIdx, isIndex))
                    {
                        lambdaW = lambdaWI;
                        if (compressibility_)
                        {
                            lambdaW /= cellDataI.density(wPhaseIdx);
                        } //divide by density because lambda is saved as lambda*density
                    }
                    else
                    {
                        if (compressibility_)
                        {
                            lambdaW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satWBound)
                                    / FluidSystem::viscosity(cellDataI.fluidState(), wPhaseIdx);
                        }
                        else
                        {
                            lambdaW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satWBound)
                                    / viscosityW;
                        }
                    }

                    if (cellDataI.fluxData().isUpwindCell(nPhaseIdx, isIndex))
                    {
                        lambdaNW = lambdaNWI;
                        if (compressibility_)
                        {
                            lambdaNW /= cellDataI.density(nPhaseIdx);
                        } //divide by density because lambda is saved as lambda*density
                    }
                    else
                    {
                        if (compressibility_)
                        {
                            lambdaNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satWBound)
                                    / FluidSystem::viscosity(cellDataI.fluidState(), nPhaseIdx);
                        }
                        else
                        {
                            lambdaNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satWBound)
                                    / viscosityNW;
                        }
                    }
                    //                    std::cout<<lambdaW<<" "<<lambdaNW<<std::endl;

                    switch (velocityType_)
                    {
                    case vt:
                    {
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorTotal, *isIt);

                        Scalar pcI = cellDataI.capillaryPressure();

                        // calculate the saturation gradient
                        GlobalPosition pcGradient = unitOuterNormal;
                        pcGradient *= (pcI - pcBound) / dist;

                        // get the diffusive part -> give 1-sat because sat = S_n and lambda = lambda(S_w) and pc = pc(S_w)
                        Scalar diffPart = diffusivePart()(*eIt, isIndex, satWI, satWBound, pcGradient) * unitOuterNormal
                                * faceArea;

                        Scalar convPart = convectivePart()(*eIt, isIndex, satWI, satWBound) * unitOuterNormal
                                * faceArea;

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
                            diffPart *= -1;                        //add cflFlux for time-stepping
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorW, *isIt,
                                    wPhaseIdx);
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorNW,
                                    *isIt, nPhaseIdx);
                            convPart *= -1;
                            break;
                        }
                        }
                        //vt*fw
                        factorTotal -= diffPart;
                        factorTotal += convPart;

                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, 10 * diffPart, *isIt);
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, 10 * convPart, *isIt);

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
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorW, *isIt,
                                    wPhaseIdx);
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorNW, *isIt,
                                    nPhaseIdx);

                        break;
                    }
                    }
                } //end dirichlet boundary

                if (bcType.isNeumann(satEqIdx))
                {
                    problem_.neumann(boundValues, *isIt);
                    factorW = boundValues[wPhaseIdx];
                    factorNW = boundValues[nPhaseIdx];
                    factorW *= faceArea;
                    factorNW *= faceArea;

                    //get mobilities
                    Scalar lambdaW, lambdaNW;

                    lambdaW = lambdaWI;
                    lambdaNW = lambdaNWI;
                    if (compressibility_)
                    {
                        lambdaW /= cellDataI.density(wPhaseIdx);
                        lambdaNW /= cellDataI.density(nPhaseIdx);
                        factorW /= cellDataI.density(wPhaseIdx);
                        factorNW /= cellDataI.density(nPhaseIdx);
                    }
                    else
                    {
                        factorW /= densityW;
                        factorNW /= densityNW;
                    }

                    switch (velocityType_)
                    {
                    case vt:
                    {
                        //add cflFlux for time-stepping
                        evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorW + factorNW,
                                *isIt);
                        break;
                    }
                    default:
                    {
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorW, *isIt,
                                    wPhaseIdx);
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorNW, *isIt,
                                    nPhaseIdx);

                        break;
                    }
                    }

                } //end neumann boundary
                if (bcType.isOutflow(satEqIdx))
                {
                    //get mobilities
                    Scalar lambdaW = lambdaWI;
                    Scalar lambdaNW = lambdaNWI;
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
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorTotal,
                                    *isIt);
                            break;
                        }
                        case Sn:
                        {
                            //vt*fn
                            factorTotal *= lambdaNW / (lambdaW + lambdaNW);
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorTotal,
                                    *isIt);
                            break;
                        }
                        }
                    }
                    else
                    {
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorW, *isIt,
                                    wPhaseIdx);
                            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, factorNW, *isIt,
                                    nPhaseIdx);

                        break;
                    }

                } //end outflow boundary
            } //end boundary
              // add to update vector

            switch (velocityType_)
            {
            case vt:
                updateVec[globalIdxI] -= factorTotal / (volume * porosity); //-:v>0, if flow leaves the cell
                break;
            default:
                switch (saturationType_)
                {
                case Sw:
                    updateVec[globalIdxI] -= factorW / (volume * porosity); //-:v>0, if flow leaves the cell
                    break;
                case Sn:
                    updateVec[globalIdxI] -= factorNW / (volume * porosity); //-:v>0, if flow leaves the cell
                    break;
                }
                break;
            }
        } // end all intersections
        PrimaryVariables sourceVec(0.0);
        problem_.source(sourceVec, *eIt);

        if (compressibility_)
        {
            sourceVec[wPhaseIdx] /= cellDataI.density(wPhaseIdx);
            sourceVec[nPhaseIdx] /= cellDataI.density(nPhaseIdx);
        }
        else
        {
            sourceVec[wPhaseIdx] /= densityW;
            sourceVec[nPhaseIdx] /= densityNW;
        }

        //get mobilities
        Scalar lambdaW = lambdaWI;
        Scalar lambdaNW = lambdaNWI;
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

            updateVec[globalIdxI] += sourceVec[wPhaseIdx] / porosity;
            break;
        }
        case Sn:
        {
            if (sourceVec[nPhaseIdx] < 0 && cellDataI.saturation(nPhaseIdx) < threshold_)
                sourceVec[nPhaseIdx] = 0.0;

            updateVec[globalIdxI] += sourceVec[nPhaseIdx] / porosity;
            break;
        }
        }

        switch (velocityType_)
        {
        case vt:
        {
            //add cflFlux for time-stepping
            evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW,
                    (sourceVec[wPhaseIdx] + sourceVec[nPhaseIdx]) * volume, *eIt);
            break;
        }
        default:
        {
                //add cflFlux for time-stepping
                evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, sourceVec[wPhaseIdx] * volume,
                        *eIt, wPhaseIdx);
                evalCflFluxFunction().addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, sourceVec[nPhaseIdx] * volume,
                        *eIt, nPhaseIdx);
            break;
        }
        }

        //calculate time step
        dt = std::min(dt, evalCflFluxFunction().getCflFluxFunction(globalPos, *eIt) * (porosity * volume));

        cellDataI.setUpdate(updateVec[globalIdxI]);
    } // end grid traversal
}

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

        // initialize cell concentration
        problem_.variables().primaryVariablesGlobal(satEqIdx)[globalIdx] = initSol[saturationIdx];
    }

    return;
}

template<class TypeTag>
void FVSaturation2P<TypeTag>::updateMaterialLaws(RepresentationType& saturation = *(new RepresentationType(0)),
        bool iterate = false)
{
    ElementIterator eItBegin = problem_.gridView().template begin<0>();

    Scalar temp = problem_.temperature(*eItBegin);
    Scalar pRef = problem_.referencePressure(*eItBegin);
//    Scalar densityW = WettingPhase::density(temp, pRef);
//    Scalar densityNW = NonwettingPhase::density(temp, pRef);
    Scalar viscosityW = WettingPhase::viscosity(temp, pRef);
    Scalar viscosityNW = NonwettingPhase::viscosity(temp, pRef);

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem_.variables().index(*eIt);

        CellData& cellData = problem_.variables().cellData(globalIdx);

//        Scalar temperature = problem_.temperature(*eIt);

        //determine phase saturations from primary saturation variable
        Scalar satW = 0;
        Scalar satNW = 0;
        switch (saturationType_)
        {
        case Sw:
        {
            satW = problem_.variables().primaryVariablesGlobal(satEqIdx)[globalIdx];
            satNW = 1 - satW;
            break;
        }
        case Sn:
        {
            satNW = problem_.variables().primaryVariablesGlobal(satEqIdx)[globalIdx];
            satW = 1 - satNW;
            break;
        }
        }

        Scalar pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

        cellData.setSaturation(wPhaseIdx, satW);
        cellData.setSaturation(nPhaseIdx, satNW);
        cellData.setCapillaryPressure(pc);

        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosityW;
        Scalar mobilityNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosityNW;

        // initialize mobilities
        cellData.mobility(wPhaseIdx) = mobilityW;
        cellData.mobility(nPhaseIdx) = mobilityNW;

        //initialize fractional flow functions
        cellData.fracFlowFunc(wPhaseIdx) = mobilityW / (mobilityW + mobilityNW);
        cellData.fracFlowFunc(nPhaseIdx) = mobilityNW / (mobilityW + mobilityNW);
    }
    return;
}

}
#endif
