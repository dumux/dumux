// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_VARIABLECLASS2P_HH
#define DUMUX_VARIABLECLASS2P_HH

#include <dumux/decoupled/common/variableclass.hh>
#include <dumux/decoupled/2p/2pfluidstate.hh>
#include "2pproperties.hh"

/**
 * @file
 * @brief  Class including the variables and data of discretized data of the constitutive relations
 * @author Markus Wolff
 */

namespace Dumux
{
/*!
 * \ingroup DecoupledModel
 */
//! Class including the variables and data of discretized data of the constitutive relations.
/*! The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the decoupled two-phase flow is stored, as well as discretized data of constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure. Thus, they have to be callculated just once in every time step or every iteration step.
 *
 * Template parameters are:

 - GridView a DUNE gridview type
 - Scalar type used for scalar quantities
 */
template<class TypeTag>
class VariableClass2P: public VariableClass<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef VariableClass<TypeTag> ParentClass;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW, pn = Indices::pressureNW, pglobal = Indices::pressureGlobal,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation));

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;//!<type for vector of scalars
    typedef typename SolutionTypes::PhaseProperty PhasePropertyType;//!<type for vector of phase properties
    typedef typename SolutionTypes::FluidProperty FluidPropertyType;//!<type for vector of fluid properties
    typedef typename SolutionTypes::PhasePropertyElemFace PhasePropertyElemFaceType;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef typename SolutionTypes::DimVecElemFace DimVecElemFaceType;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars

public:
    const int codim_;

    ScalarSolutionType saturation_;
    PhasePropertyType mobility_;//store lambda for efficiency reasons
    PhasePropertyType fracFlowFunc_;
    ScalarSolutionType capillaryPressure_;
    DimVecElemFaceType velocitySecondPhase_;

    FluidPropertyType density_;
    FluidPropertyType viscosity_;

    ScalarSolutionType volumecorrection_;

public:

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */

    VariableClass2P(const GridView& gridView, Scalar& initialSat = *(new Scalar(1)),
            Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim>(0))) :
        VariableClass<TypeTag> (gridView, initialVel), codim_(0)
    {
        initializeGlobalVariablesDiffPart(initialVel);
        initializeGlobalVariablesTransPart(initialSat);
    }

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param codim codimension of the entity of which data has to be strored
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass2P(const GridView& gridView, int codim, Scalar& initialSat = *(new Scalar(1)), Dune::FieldVector<
            Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim>(0))) :
        VariableClass<TypeTag> (gridView, codim, initialVel), codim_(codim)
    {
        initializeGlobalVariablesDiffPart(initialVel);
        initializeGlobalVariablesTransPart(initialSat);
    }

    // serialization methods
    template<class Restarter>
    void serialize(Restarter &res)
    {
        res.template serializeEntities<0> (*this, this->gridView());
    }
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        res.template deserializeEntities<0> (*this, this->gridView());
    }

    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int globalIdx = this->elementMapper().map(element);
        outstream << this->pressure()[globalIdx] << "  " << saturation_[globalIdx];
    }
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int globalIdx = this->elementMapper().map(element);
        instream >> this->pressure()[globalIdx] >> saturation_[globalIdx];
    }

private:
    void initializeGlobalVariablesDiffPart(Dune::FieldVector<Scalar, dim>& initialVel)
    {
        //resize to grid size
        velocitySecondPhase_.resize(this->gridSize());//depends on pressure
        for (int i=0; i<2; i++) //for both phases
        {
        density_[i].resize(this->gridSize());//depends on pressure
        viscosity_[i].resize(this->gridSize());//depends on pressure
        }
        //initialise variables
        velocitySecondPhase_ = initialVel;
    }
    void initializeGlobalVariablesTransPart(int initialSat)
    {
        //resize to grid size
        saturation_.resize(this->gridSize());
        for (int i=0; i<2; i++) //for both phases
        {
        mobility_[i].resize(this->gridSize());//lambda is dependent on saturation! ->choose same size
        fracFlowFunc_[i].resize(this->gridSize());//depends on saturation
        }
        capillaryPressure_.resize(this->gridSize());//depends on saturation
        volumecorrection_.resize(this->gridSize());//dS/dt for correction of pressure equation

        //initialise variables
        saturation_ = initialSat;
    }

public:
    //Write saturation and pressure into file
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        if (codim_ == 0)
        {
            ScalarSolutionType *pressure = writer.template createField<Scalar, 1> (this->gridSize());
            ScalarSolutionType *saturation = writer.template createField<Scalar, 1> (this->gridSize());

            *pressure = this->pressure();
            *saturation = saturation_;

            writer.addCellData(pressure, "pressure");
            writer.addCellData(saturation, "saturation");

//            // exemplary output for phase-dependent stuff
//            ScalarSolutionType *fractionalflowW = writer.template createField<Scalar, 1> (this->gridSize());
//            *fractionalflowW = fracFlowFunc_[wPhaseIdx];
//            writer.addCellData(fractionalflowW, "fractional flow function w-phase");
        }
        if (codim_ == dim)
        {
            ScalarSolutionType *pressure = writer.template createField<Scalar, 1> (this->gridSize());
            ScalarSolutionType *saturation = writer.template createField<Scalar, 1> (this->gridSize());

            *pressure = this->pressure();
            *saturation = saturation_;

            writer.addVertexData(pressure, "pressure");
            writer.addVertexData(saturation, "saturation");
        }

        return;
    }

    //! Return the vector of the transported quantity, which is the saturation for an IMPES scheme
    ScalarSolutionType& transportedQuantity()
    {
        return saturation_;
    }

    //! Return saturation vector
    ScalarSolutionType& saturation()
    {
        return saturation_;
    }

    const ScalarSolutionType& saturation() const
    {
        return saturation_;
    }

    //! Return velocity vector
    DimVecElemFaceType& velocitySecondPhase()
    {
        return velocitySecondPhase_;
    }

    const DimVecElemFaceType& velocitySecondPhase() const
    {
        return velocitySecondPhase_;
    }

    //! Return vector of wetting phase potential gradients
    Scalar& potentialWetting(int Idx1, int Idx2)
    {
        return this->potential(Idx1, Idx2)[wPhaseIdx];
    }

    const Scalar& potentialWetting(int Idx1, int Idx2) const
    {
        return this->potential(Idx1, Idx2)[wPhaseIdx];
    }

    //! Return vector of non-wetting phase potential gradients
    Scalar& potentialNonwetting(int Idx1, int Idx2)
    {
        return this->potential(Idx1, Idx2)[nPhaseIdx];
    }

    const Scalar& potentialNonwetting(int Idx1, int Idx2) const
    {
        return this->potential(Idx1, Idx2)[nPhaseIdx];
    }

    //! Return vector of wetting phase mobilities
    Scalar& mobilityWetting(int Idx)
    {
        return mobility_[wPhaseIdx][Idx][0];
    }

    const Scalar& mobilityWetting(int Idx) const
    {
        return mobility_[wPhaseIdx][Idx][0];
    }

    //! Return vector of non-wetting phase mobilities
    Scalar& mobilityNonwetting(int Idx)
    {
        return mobility_[nPhaseIdx][Idx][0];
    }

    const Scalar& mobilityNonwetting(int Idx) const
    {
        return mobility_[nPhaseIdx][Idx][0];
    }

    //! Return vector of wetting phase fractional flow functions
    Scalar& fracFlowFuncWetting(int Idx)
    {
        return fracFlowFunc_[wPhaseIdx][Idx][0];
    }

    const Scalar& fracFlowFuncWetting(int Idx) const
    {
        return fracFlowFunc_[wPhaseIdx][Idx][0];
    }

    //! Return vector of non-wetting phase fractional flow functions
    Scalar& fracFlowFuncNonwetting(int Idx)
    {
        return fracFlowFunc_[nPhaseIdx][Idx][0];
    }

    const Scalar& fracFlowFuncNonwetting(int Idx) const
    {
        return fracFlowFunc_[nPhaseIdx][Idx][0];
    }

    //! Return capillary pressure vector
    Scalar& capillaryPressure(int Idx)
    {
        return capillaryPressure_[Idx][0];
    }

    const Scalar& capillaryPressure(int Idx) const
    {
        return capillaryPressure_[Idx][0];
    }

    //! Return density vector
    Scalar& densityWetting(int Idx)
    {
        return density_[wPhaseIdx][Idx][0];
    }
    const Scalar& densityWetting(int Idx) const
    {
        return density_[wPhaseIdx][Idx][0];
    }

    //! Return density vector
    Scalar& densityNonwetting(int Idx)
    {
        return density_[nPhaseIdx][Idx][0];
    }

    const Scalar& densityNonwetting(int Idx) const
    {
        return density_[nPhaseIdx][Idx][0];
    }

    //! Return density vector
    Scalar& viscosityWetting(int Idx)
    {
        return viscosity_[wPhaseIdx][Idx][0];
    }

    const Scalar& viscosityWetting(int Idx) const
    {
        return viscosity_[wPhaseIdx][Idx][0];
    }

    //! Return density vector
    Scalar& viscosityNonwetting(int Idx)
    {
        return viscosity_[nPhaseIdx][Idx][0];
    }

    const Scalar& viscosityNonwetting(int Idx) const
    {
        return viscosity_[nPhaseIdx][Idx][0];
    }

    Scalar& volumecorrection(int Idx)
    {
        return volumecorrection_[Idx][0];
    }

    const Scalar& volumecorrection(int Idx) const
    {
        return volumecorrection_[Idx][0];
    }

    //! Get saturation
    /*! evaluate saturation at given element
     @param element entity of codim 0
     \return value of saturation
     */
    Dune::FieldVector<Scalar, 1>& satElement(const Element& element)
    {
        return saturation_[this->elementMapper().map(element)];;
    }

    const Dune::FieldVector<Scalar, 1>& satElement(const Element& element) const
    {
        return saturation_[this->elementMapper().map(element)];;
    }

};
}
#endif
