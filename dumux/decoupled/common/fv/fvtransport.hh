// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
/*!\brief The finite volume discretization of a transport equation
 *  Base class for finite volume (FV) implementations of an explicitly treated transport equation.
 *  The class provides a method to calculate the explicit update to get a new solution of the transported quantity:
 *  \f[
 *      u_{new} = u_{old} + \Delta t \Delta u_{update}
 *  \f]
 *  A certain transport equation defined in a implementation of this base class must be splitted into a flux term and a source term.
 *  Corresponding functions (<tt>getSource()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt>) have to be defined in the implementation.
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
    //! \cond \private
    EvalCflFluxFunction& evalCflFluxFunction()
    {
        return *evalCflFluxFunction_;
    }

    const EvalCflFluxFunction& evalCflFluxFunction() const
    {
        return *evalCflFluxFunction_;
    }
    //! \endcond

public:

    // Calculate the update vector.
    void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet);

    /*! \brief Function which calculates the flux update
     *
     * Function computes the inter-cell flux term and adds it to the update.
     *
     * \param update The cell update
     * \param intersection Intersection of two grid elements
     * \param cellDataI Object containing all model relevant cell data
     */
    void getFlux(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    /*! \brief Function which calculates the boundary flux update
     *
     * Function computes the boundary-flux term and  adds it to the update.
     *
     * \param update The cell update
     * \param intersection Intersection of two grid elements
     * \param cellDataI Object containing all model relevant cell data
     */
    void getFluxOnBoundary(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    /*! \brief Function which calculates the source update
     *
     * Function computes the source term and adds it to the update.
     *
     * \param update The cell update
     * \param element Grid element
     * \param cellDataI Object containing all model relevant cell data
     */
    void getSource(Scalar& update, const Element& element, CellData& cellDataI);

    //! Sets the initial solution \f$ S_0 \f$.
    void initialize();


    /*! \brief Updates constitutive relations and stores them in the variable class*/
    void updateMaterialLaws();

    /*! \brief Writes the current values of the primary transport variable into the <tt>transportedQuantity</tt>-vector (comes as function argument)
     *
     * \param transportedQuantity Vector of the size of global numbers of degrees of freedom of the primary transport variable.
     */
    void getTransportedQuantity(TransportSolutionType& transportedQuantity);

    /*! \brief Updates the primary transport variable.
     *
     * \param updateVec Vector containing the global update.
     */
    void updateTransportedQuantity(TransportSolutionType& updateVec);

    /*! \brief Adds transport output to the output file
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {}

    /*! \brief  Function for serialization of the primary transport variable.
     *
     *  Function needed for restart option. Writes the primary transport variable of a grid element to a restart file.
     *
     *  \param outstream Stream into the restart file.
     *  \param element Grid element
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {}

    /*! \brief  Function for deserialization of the primary transport variable.
     *
     *  Function needed for restart option. Reads the the primary transport variable of a grid element from a restart file.
     *
     *  \param instream Stream from the restart file.
     *  \param element Grid element
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {}


    //! Constructs a FVTransport object
    /**

     * \param problem A problem class object
     */
    FVTransport(Problem& problem) :
            problem_(problem), switchNormals_(GET_PARAM(TypeTag, bool, SwitchNormals))
    {
        evalCflFluxFunction_ = new EvalCflFluxFunction(problem);
    }

    //! Destructor
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

/*! \brief Calculate the update vector.
 *  \param t         current time
 *  \param dt        time step size
 *  \param updateVec  vector containing the update values
 *  \param impet      variable should be true if an impet algorithm is used and false if the transport part is solved independently
 *
 *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
 *  employing a CFL condition.
 */
template<class TypeTag>
void FVTransport<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false)
{
    if (!impet)
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
