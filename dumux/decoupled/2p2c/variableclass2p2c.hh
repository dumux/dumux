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
#ifndef DUMUX_VARIABLECLASS2P2C_HH
#define DUMUX_VARIABLECLASS2P2C_HH

#include <dumux/decoupled/common/variableclass.hh>
#include <dumux/decoupled/2p2c/dec2p2cfluidstate.hh>

/**
 * @file
 * @brief  Class including the variables and data of discretized data of the constitutive relations
 * @author Markus Wolff, Benjamin Faigle
 */

namespace Dumux
{
/*!
 * \ingroup IMPEC
 */
//! Class including the variables and data of discretized data of the constitutive relations.
/*! The variables of compositional two-phase flow, which are one pressure and one saturation are stored in this class.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class VariableClass2P2C: public VariableClass<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;
    typedef typename SolutionTypes::ElementMapper ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;

    typedef VariableClass<TypeTag> ParentClass;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW, pn = Indices::pressureNW, pglobal = Indices::pressureGlobal
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx
    };

    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;//!<type for vector of scalars
    typedef typename SolutionTypes::PhaseProperty PhasePropertyType;//!<type for vector of phase properties = FluidPropertyType
    typedef typename SolutionTypes::FluidProperty FluidPropertyType;//!<type for vector of fluid properties = PhasePropertyType -> decoupledproperties.hh
    typedef typename SolutionTypes::PhasePropertyElemFace PhasePropertyElemFaceType;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef typename SolutionTypes::DimVecElemFace DimVecElemFaceType;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars

protected:
    const int codim_;

    ScalarSolutionType saturation_;
    PhasePropertyType mobility_; //store lambda for efficiency reasons
    ScalarSolutionType capillaryPressure_;
    DimVecElemFaceType velocitySecondPhase_;    //necessary??

    FluidPropertyType density_;
    FluidPropertyType numericalDensity_;
    FluidPropertyType viscosity_;

    ScalarSolutionType volErr_;

    TransportSolutionType totalConcentration_;
    PhasePropertyType massfrac_;

    //vectors for partial derivatives
    PhasePropertyType dv_; // although it is a component property
    ScalarSolutionType dv_dp_;

    TransportSolutionType updateEstimate_;
    ScalarSolutionType errorCorrection_;
    ScalarSolutionType perimeter_;
    Dune::BlockVector<Dune::FieldVector<int,1> > subdomain_;

public:

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */

    VariableClass2P2C(const GridView& gridView, Scalar& initialSat = *(new Scalar(1)),
            Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim>(0))) :
        VariableClass<TypeTag> (gridView, initialVel), codim_(0)
    {
        initialize2p2cVariables(initialVel, initialSat);
    }

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param codim codimension of the entity of which data has to be strored
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass2P2C(const GridView& gridView, int codim, Scalar& initialSat = *(new Scalar(1)), Dune::FieldVector<
            Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim>(0))) :
        VariableClass<TypeTag> (gridView, codim, initialVel), codim_(codim)
    {
        initialize2p2cVariables(initialVel, initialSat);
    }

    //! serialization methods -- same as Dumux::VariableClass2P
    //@{
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
        outstream << this->pressure()[globalIdx] << "  " << totalConcentration(globalIdx, 0)
                  << "  " << totalConcentration(globalIdx,1);
    }
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int globalIdx = this->elementMapper().map(element);
        instream >> this->pressure()[globalIdx] >> totalConcentration(globalIdx, 0)
                 >> totalConcentration(globalIdx,1);
    }
    //@}


private:
    //! initializes the miscible 2p variables
    /*! Internal method that prepares the vectors holding the 2p2c variables for the current
     *  simulation.
     *
     *\param initialVel Vector containing the initial velocity
     *\param initialSat Initial value for the saturation
     */
    void initialize2p2cVariables(Dune::FieldVector<Scalar, dim>& initialVel,
                                 Scalar initialSat)
    {
        //resize to grid size
        int size_ = this->gridSize();
        // a) global variables
        velocitySecondPhase_.resize(size_);//depends on pressure
            velocitySecondPhase_ = initialVel;
        for (int i=0; i<2; i++)    //for both phases
        {
        density_[i].resize(size_);//depends on pressure
        numericalDensity_[i].resize(size_);//depends on pressure
        viscosity_[i].resize(size_);//depends on pressure
        };
        // b) transport variables
        saturation_.resize(size_);
            saturation_ = initialSat;
        capillaryPressure_.resize(size_);//depends on saturation
        mobility_[0].resize(size_);//lambda is dependent on saturation! ->choose same size
        mobility_[1].resize(size_);

        // c) compositional stuff
        totalConcentration_.resize(2);  // resize solution vector to 2 primary transport variables
        updateEstimate_.resize(2);
        for (int i=0; i<2; i++) //for both phases
        {
        totalConcentration_[i].resize(size_);   // resize vector to number of cells
        updateEstimate_[i].resize(size_);
        massfrac_[i].resize(size_);
        dv_[i].resize(size_);
        }
        dv_dp_.resize(size_);
        volErr_.resize(size_);
        errorCorrection_.resize(size_);
        perimeter_.resize(size_);
    }


public:

    //! Writes output into file
    /*!
     * Creates output for standard 2p2c models.
     * \param writer Applied VTK-writer
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size_ = this->gridSize();
        if (codim_ == 0)
        {
            // pressure & saturation
            ScalarSolutionType *pressure = writer.allocateManagedBuffer (size_);
            ScalarSolutionType *saturation = writer.allocateManagedBuffer (size_);
            *pressure = this->pressure();
            *saturation = saturation_;
            if (GET_PROP_VALUE(TypeTag, PressureFormulation) == GET_PROP_TYPE(TypeTag, TwoPTwoCIndices)::pressureW)
                writer.attachCellData(*pressure, "pressure w_phase");
            else
                writer.attachCellData(*pressure, "pressure nw_phase");
            writer.attachCellData(*saturation, "water saturation");
            if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
            {
                ScalarSolutionType *pc = writer.allocateManagedBuffer (size_);
                *pc = capillaryPressure_;
                writer.attachCellData(*pc, "capillary pressure");
            }

            //total Concentration
            ScalarSolutionType *totalConcentration1 = writer.allocateManagedBuffer (size_);
            ScalarSolutionType *totalConcentration2 = writer.allocateManagedBuffer (size_);
            *totalConcentration1 = totalConcentration_[0];
            *totalConcentration2 = totalConcentration_[1];
            writer.attachCellData(*totalConcentration1, "totalConcentration w_Comp");
            writer.attachCellData(*totalConcentration2, "totalConcentration n_Comp");

            // numerical stuff
            ScalarSolutionType *volErr = writer.allocateManagedBuffer (size_);
            *volErr = volErr_;
            writer.attachCellData(*volErr, "volume Error");

            // composition
            ScalarSolutionType *massfraction1W = writer.allocateManagedBuffer (size_);
            ScalarSolutionType *massfraction1NW = writer.allocateManagedBuffer (size_);
            *massfraction1W = massfrac_[0];
            *massfraction1NW = massfrac_[1];
            writer.attachCellData(*massfraction1W, "massfraction1 in w_phase");
            writer.attachCellData(*massfraction1NW, "massfraction1NW nw_phase");

            // phase properties
            ScalarSolutionType *mobilityW = writer.allocateManagedBuffer (size_);
            ScalarSolutionType *mobilityNW = writer.allocateManagedBuffer (size_);
            ScalarSolutionType *densityW = writer.allocateManagedBuffer (size_);
            ScalarSolutionType *densityNW = writer.allocateManagedBuffer (size_);
            *mobilityW = mobility_[0];
            *mobilityNW = mobility_[1];
            *densityW = density_[0];
            *densityNW = density_[1];
            writer.attachCellData(*mobilityW, "mobility w_phase");
            writer.attachCellData(*mobilityNW, "mobility nw_phase");
            writer.attachCellData(*densityW, "density w_phase");
            writer.attachCellData(*densityNW, "density nw_phase");
        }
        if (codim_ == dim)
        {
            ScalarSolutionType *pressure = writer.allocateManagedBuffer (size_);
            ScalarSolutionType *saturation = writer.allocateManagedBuffer (size_);

            *pressure = this->pressure();
            *saturation = saturation_;

            writer.attachVertexData(*pressure, "pressure");
            writer.attachVertexData(*saturation, "saturation");
        }

        return;
    }

    //! Communicate transported quantity to ghosts/overlap
    void communicateTransportedQuantity()
    {
#if HAVE_MPI
        ElementHandleAssign<typename TransportSolutionType::block_type::block_type, typename TransportSolutionType::block_type, ElementMapper> elementHandleW(totalConcentration_[wCompIdx], this->elementMapper());
        this->gridView().communicate(elementHandleW,
                Dune::InteriorBorder_All_Interface,
                Dune::ForwardCommunication);
        ElementHandleAssign<typename TransportSolutionType::block_type::block_type, typename TransportSolutionType::block_type, ElementMapper> elementHandleNW(totalConcentration_[nCompIdx], this->elementMapper());
        this->gridView().communicate(elementHandleNW,
                Dune::InteriorBorder_All_Interface,
                Dune::ForwardCommunication);
#endif
    }

    //! \name Access functions
    //@{
    //! Return the vector of the transported quantity
    /*! For an immiscible IMPES scheme, this is the saturation. For Miscible simulations, however,
     *  the total concentration of all components is transported.
     */
    TransportSolutionType& transportedQuantity()
    {
        return totalConcentration_;
    }

    //! Returs a reference to the total concentration vector
    const TransportSolutionType& totalConcentration() const
    {
        return totalConcentration_;
    }

    //! Returs a reference to the total concentration
    /*!\param Idx Element index
     * \param compIdx Index of the component */
    Scalar& totalConcentration(int Idx, int compIdx)
    {
        return totalConcentration_[compIdx][Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::totalConcentration()
    const Scalar& totalConcentration(int Idx, int compIdx) const
    {
        return totalConcentration_[compIdx][Idx][0];
    }

    //! Return saturation vector
    ScalarSolutionType& saturation()
    {
        return saturation_;
    }
    //! Return saturation vector
    const ScalarSolutionType& saturation() const
    {
        return saturation_;
    }
    //! \copydoc Dumux::VariableClass2P2C::saturation()
    Scalar& saturation(int Idx)
    {
        return saturation_[Idx][0];
    }

    //! Return vector of wetting phase mobilities
    /*! \param Idx Element index*/
    Scalar& mobilityWetting(int Idx)
    {
        return mobility_[wPhaseIdx][Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::mobilityWetting()
    const Scalar& mobilityWetting(int Idx) const
    {
        return mobility_[wPhaseIdx][Idx][0];
    }

    //! Return vector of non-wetting phase mobilities
    /*! \param Idx Element index*/
    Scalar& mobilityNonwetting(int Idx)
    {
        return mobility_[nPhaseIdx][Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::mobilityNonwetting()
    const Scalar& mobilityNonwetting(int Idx) const
    {
        return mobility_[nPhaseIdx][Idx][0];
    }

    //! Return capillary pressure vector
    const ScalarSolutionType& capillaryPressure() const
    {
        return capillaryPressure_;
    }
    //! Return capillary pressure of specific Element
    /*! \param Idx Element index*/
    Scalar& capillaryPressure(int Idx)
    {
        return capillaryPressure_[Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::capillaryPressure()
    const Scalar& capillaryPressure(int Idx) const
    {
        return capillaryPressure_[Idx][0];
    }

    //! Return wetting phase density
    /*! \param Idx Element index*/
    Scalar& densityWetting(int Idx)
    {
        return density_[wPhaseIdx][Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::densityWetting()
    const Scalar& densityWetting(int Idx) const
    {
        return density_[wPhaseIdx][Idx][0];
    }

    //! Return nonwetting phase density
    /*! \param Idx Element index*/
    Scalar& densityNonwetting(int Idx)
    {
        return density_[nPhaseIdx][Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::densityNonwetting()
    const Scalar& densityNonwetting(int Idx) const
    {
        return density_[nPhaseIdx][Idx][0];
    }

    /*! \brief Returns numerical phase density (current mass/volume, not from EOS)
     * \param Idx Element index
     * \param phaseIdx Index of the phase
     */
    Scalar& numericalDensity(int Idx, int phaseIdx)
    {
        return numericalDensity_[phaseIdx][Idx][0];
    }

    /*! \brief Returns numerical phase density (current mass/volume, not from EOS)
     * \param phaseIdx Index of the phase
     */
    const ScalarSolutionType& numericalDensity(int phaseIdx) const
    {
        return numericalDensity_[phaseIdx];
    }

    //! Returns correction term for volume mismatch \f$\varepsilon = \frac{v_t - \phi}{\Delta t}\f$
    /*! \param Idx Element index*/
    Scalar& errorCorrection(int Idx)
    {
        return errorCorrection_[Idx][0];
    }

    //! Returns correction term vector for volume mismatch \f$\varepsilon\f$
    const ScalarSolutionType& errorCorrection() const
    {
        return errorCorrection_;
    }

    //! Return viscosity of the wetting phase
    /*! \param Idx Element index*/
    Scalar& viscosityWetting(int Idx)
    {
        return viscosity_[wPhaseIdx][Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::viscosityWetting()
    const Scalar& viscosityWetting(int Idx) const
    {
        return viscosity_[wPhaseIdx][Idx][0];
    }

    //! Return viscosity of the nonwetting phase
    /*! \param Idx Element index*/
    Scalar& viscosityNonwetting(int Idx)
    {
        return viscosity_[nPhaseIdx][Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::viscosityNonwetting()
    const Scalar& viscosityNonwetting(int Idx) const
    {
        return viscosity_[nPhaseIdx][Idx][0];
    }

    //! Returns a reference to the vector holding the volume error
    ScalarSolutionType& volErr()
    {
        return volErr_;
    }
    //! Returns a reference to the volume error
    /*! \param Idx Element index*/
    const Scalar& volErr(int Idx) const
    {
        return volErr_[Idx][0];
    }


    //! Returns a reference to mass fraction of first component in wetting phase
    /*! \param Idx Element index*/
    Scalar& wet_X1(int Idx)
    {
        return massfrac_[wPhaseIdx][Idx][0];
    }
    //! Returns a reference to mass fraction of first component in nonwetting phase
    /*! \param Idx Element index*/
    Scalar& nonwet_X1(int Idx)
    {
        return massfrac_[nPhaseIdx][Idx][0];
    }
    //! Returns a reference to mass fractions of first component
    /*!\param Idx Element index
     * \param phaseIdx Index of the phase */
    const Scalar& massFracComp1(int Idx, int phaseIdx) const
    {
        return massfrac_[phaseIdx][Idx][0];
    }

    //! Return the estimated update vector for pressure calculation
    TransportSolutionType& updateEstimate()
    {
        return updateEstimate_;
    }
    //! Return the estimate of the next update
    /*!\param Idx Element index
     * \param compIdx Index of the component */
    Scalar& updateEstimate(int Idx, int compIdx)
    {
        return updateEstimate_[compIdx][Idx][0];
    }

    //! Communicate the estimate of the next update to ghost/overlap cells
    void communicateUpdateEstimate()
    {
#if HAVE_MPI
        ElementHandleAssign<typename TransportSolutionType::block_type::block_type, typename TransportSolutionType::block_type, ElementMapper> elementHandleW(updateEstimate_[wCompIdx], this->elementMapper());
        this->gridView().communicate(elementHandleW,
                Dune::InteriorBorder_All_Interface,
                Dune::ForwardCommunication);
        ElementHandleAssign<typename TransportSolutionType::block_type::block_type, typename TransportSolutionType::block_type, ElementMapper> elementHandleNW(updateEstimate_[nCompIdx], this->elementMapper());
        this->gridView().communicate(elementHandleNW,
                Dune::InteriorBorder_All_Interface,
                Dune::ForwardCommunication);
#endif
    }

    //! Returs a reference to the volume derivatives vector
    const PhasePropertyType& dv() const
    {
        return dv_;
    }
    //! Returs a reference to the volume derivatives
    /*!\param Idx Element index
     * \param compIdx Index of the component */
    Scalar& dv(int Idx, int compIdx)
    {
        return dv_[compIdx][Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::dv()
    const Scalar& dv(int Idx, int compIdx) const
    {
        return dv_[compIdx][Idx][0];
    }

    //! Returs a reference to the volume derivatives (wrt pressure) vector
    const ScalarSolutionType& dv_dp() const
    {
        return dv_dp_;
    }
    //! Returs a reference to the volume derivatives (wrt pressure)
    /*!\param Idx Element index*/
    Scalar& dv_dp(int Idx)
    {
        return dv_dp_[Idx][0];
    }
    //! \copydoc Dumux::VariableClass2P2C::dv_dp()
    const Scalar& dv_dp(int Idx) const
    {
        return dv_dp_[Idx][0];
    }

    //! Return the vector holding subdomain information
    Dune::BlockVector<Dune::FieldVector<int,1> >& subdomain()
    {
        return subdomain_;
    }
    //! Return specific subdomain information
    /*!\param Idx Element index */
    int& subdomain(int Idx)
    {
        return subdomain_[Idx][0];
    }

    //! Return the perimeter of each element
    /*!\param Idx Element index */
    double& perimeter(int Idx)
    {
        return perimeter_[Idx][0];
    }
    //! Returs a reference to the perimeter vector
    const ScalarSolutionType& perimeter() const
    {
        return perimeter_;
    }
    //@}
};
}
#endif
