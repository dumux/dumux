// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
#ifndef DUMUX_FVTRANSPORT_HH
#define DUMUX_FVTRANSPORT_HH

#include <dune/grid/common/gridenums.hh>
#include <dumux/decoupled/common/transportproperties.hh>
#include <dumux/decoupled/common/decoupledproperties.hh>

/**
 * @file
 * @brief  Finite Volume discretization of a  transport equation
 * @author Markus Wolff
 */

namespace Dumux
{
//! \ingroup IMPET
//! \brief The finite volume discretization of a transport equation
/*!
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport
{
    typedef typename GET_PROP_TYPE(TypeTag, TransportModel) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GET_PROP_TYPE(TypeTag, EvalCflFluxFunction) EvalCflFluxFunction;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dim> FieldVector;

protected:
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
     *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
     *  employing a CFL condition.
     */
    void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impes);

    void getFlux(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    void getFluxOnBoundary(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    void getSource(Scalar& update, const Element& element, CellData& cellDataI);

    //! Sets the initial solution \f$S_0\f$.
    void initialize();

    //! Update the values of the material laws and constitutive relations.
    /*!
     *  Constitutive relations like capillary pressure-saturation relationships, mobility-saturation relationships... are updated and stored in the variable class
     *  of type Dumux::VariableClass2P. The update has to be done when new saturation are available.
     */
    void updateMaterialLaws();

    void getTransportedQuantity(TransportSolutionType& transportedQuantity);

    void updateTransportedQuantity(TransportSolutionType& updateVec);

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {}

//    void indicator(TransportSolutionType &indicator, Scalar &globalMin, Scalar &globalMax);


    /*! \name general methods for serialization, output */
    //@{
    // serialization methods
    //! Function needed for restart option.
    void serializeEntity(std::ostream &outstream, const Element &element)
    {}

    void deserializeEntity(std::istream &instream, const Element &element)
    {}
    //@}

    //! Constructs a FVTransport object
    /**

     * \param problem a problem class object
     */

    FVTransport(Problem& problem) :
            problem_(problem), switchNormals_(GET_PARAM(TypeTag, bool, SwitchNormals))
    {
        evalCflFluxFunction_ = new EvalCflFluxFunction(problem);
    }

    ~FVTransport()
    {
        delete evalCflFluxFunction_;
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    Problem& problem_;
    bool switchNormals_;

    EvalCflFluxFunction* evalCflFluxFunction_;
};

template<class TypeTag>
void FVTransport<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impes = false)
{
    if (!impes)
    {
        asImp_().updateMaterialLaws();
    }

    // initialize dt very large
    dt = 1E100;

    // resize update vector and set to zero
    updateVec.resize(problem_.gridView().size(0));
    updateVec = 0;

    // compute update vector
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
#if HAVE_MPI
        if (eIt->partitionType() != Dune::InteriorEntity)
        {
            continue;
        }
#endif

        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        CellData& cellDataI = problem_.variables().cellData(globalIdxI);

        Scalar update = 0;
        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            FieldVector unitOuterNormal = isIt->centerUnitOuterNormal();
            if (switchNormals_)
                unitOuterNormal *= -1.0;

            // handle interior face
            if (isIt->neighbor())
            {
                //add flux to update
                asImp_().getFlux(update, *isIt, cellDataI);
            } //end intersection with neighbor element
            // handle boundary face
            else if (isIt->boundary())
            {
                //add boundary flux to update
                asImp_().getFluxOnBoundary(update, *isIt, cellDataI);
            } //end boundary
        } // end all intersections

        //add flux update to global update vector
        updateVec[globalIdxI] += update;

        //add source to global update vector
        Scalar source = 0.;
        asImp_().getSource(source,*eIt, cellDataI);
        updateVec[globalIdxI] += source;

        //calculate time step
        dt = std::min(dt, evalCflFluxFunction().getDt(*eIt));

        //store update
        cellDataI.setUpdate(updateVec[globalIdxI]);
    } // end grid traversal
}

}
#endif
