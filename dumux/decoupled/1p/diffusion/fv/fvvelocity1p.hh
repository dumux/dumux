// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
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
#ifndef DUMUX_FVVELOCITY1P_HH
#define DUMUX_FVVELOCITY1P_HH

/**
 * @file
 * @brief  Single Phase Finite Volume Model
 * @author Markus Wolff
 */


namespace Dumux
{
//! \ingroup OnePhase
//! \brief Single Phase Finite Volume Model
/*! Calculates velocities from a known pressure field in context of a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v} = q.\f]
 * The pressure has to be given as piecewise constant cell values.
 * The velocity is calculated following  Darcy's law as
 * \f[\boldsymbol{v} = -\frac{1}{\mu} \boldsymbol{K} \left(\text{grad}\, p + \rho g  \text{grad}\, z\right),\f]
 * where, \f$p\f$ is the pressure, \f$\boldsymbol{K}\f$ the absolute permeability, \f$\mu\f$ the viscosity, \f$\rho\f$ the density and \f$g\f$ the gravity constant.
 *
 * @tparam TypeTag The Type Tag
 */

template<class TypeTag>
class FVVelocity1P
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
     typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

     typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

     typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
     typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;

     typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
     typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        pressEqIdx = Indices::pressEqIdx // only one equation!
    };

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    //! The Constructor
    /**
     * \param problem a problem class object
     */
    FVVelocity1P(Problem& problem)
    : problem_(problem), gravity_(problem.gravity())
      {
        const Element& element = *(problem_.gridView().template begin<0> ());
        Scalar temperature = problem_.temperature(element);
        Scalar referencePress = problem_.referencePressure(element);

        density_ = Fluid::density(temperature, referencePress);
        viscosity_ = Fluid::viscosity(temperature, referencePress);
      }


    //! Calculate the velocity.
    /*!
     *
     *  Given the piecewise constant pressure \f$p\f$,
     *  this method calculates the velocity field
     */
    void calculateVelocity(const Intersection&, CellData&);

    void calculateVelocityOnBoundary(const Intersection&, CellData&);

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        Dune::BlockVector<Dune::FieldVector<Scalar, dim> > &velocity = *(writer.template allocateManagedBuffer<Scalar, dim> (
                problem_.gridView().size(0)));

        // compute update vector
        ElementIterator eItEnd = problem_.gridView().template end<0>();
        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // cell index
            int globalIdx = problem_.variables().index(*eIt);

            CellData& cellData = problem_.variables().cellData(globalIdx);

            Dune::FieldVector<Scalar, 2*dim> flux(0);
            // run through all intersections with neighbors and boundary
            IntersectionIterator
            isItEnd = problem_.gridView().iend(*eIt);
            for (IntersectionIterator
                    isIt = problem_.gridView().ibegin(*eIt); isIt
                    !=isItEnd; ++isIt)
            {
                int isIndex = isIt->indexInInside();

                flux[isIndex] = isIt->geometry().volume()
                        * (isIt->centerUnitOuterNormal() * cellData.fluxData().velocity(isIndex));
            }

            Dune::FieldVector<Scalar, dim> refVelocity(0);
            refVelocity[0] = 0.5 * (flux[1] - flux[0]);
            refVelocity[1] = 0.5 * (flux[3] - flux[2]);

            typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
            const Dune::FieldVector<Scalar, dim>& localPos = ReferenceElements::general(eIt->geometry().type()).position(0,
                    0);

            // get the transposed Jacobian of the element mapping
            const FieldMatrix& jacobianInv = eIt->geometry().jacobianInverseTransposed(localPos);
            FieldMatrix jacobianT(jacobianInv);
            jacobianT.invert();

            // calculate the element velocity by the Piola transformation
            Dune::FieldVector<Scalar, dim> elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= eIt->geometry().integrationElement(localPos);

            velocity[globalIdx] = elementVelocity;
        }

        writer.attachCellData(velocity, "velocity", dim);

        return;
    }
private:
    Problem &problem_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant
    Scalar density_;
    Scalar viscosity_;
};
template<class TypeTag>
void FVVelocity1P<TypeTag>::calculateVelocity(const Intersection& intersection, CellData& cellDataI)
{
    ElementPointer elementI = intersection.inside();
    ElementPointer elementJ = intersection.outside();

    int globalIdxJ = problem_.variables().index(*elementJ);

    CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = elementI->geometry().center();
    const GlobalPosition& globalPosJ = elementJ->geometry().center();

    //get face index
    int isIndexI = intersection.indexInInside();
    int isIndexJ = intersection.indexInOutside();

    //get face normal
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosJ - globalPosI;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    // compute vectorized permeabilities
    FieldMatrix meanPermeability(0);

    problem_.spatialParameters().meanK(meanPermeability, problem_.spatialParameters().intrinsicPermeability(*elementI),
            problem_.spatialParameters().intrinsicPermeability(*elementJ));

    Dune::FieldVector < Scalar, dim > permeability(0);
    meanPermeability.mv(unitOuterNormal, permeability);

    permeability /= viscosity_;

    //calculate potential gradients
    Scalar potential = (cellDataI.pressure() - cellDataJ.pressure()) / dist;

    potential += density_ * (unitOuterNormal * gravity_);

    //store potentials for further calculations (velocity, saturation, ...)
    cellDataI.fluxData().setPotential(isIndexI, potential);
    cellDataJ.fluxData().setPotential(isIndexJ, -potential);

    //calculate the gravity term
    GlobalPosition velocity(permeability);
    velocity *= (cellDataI.pressure() - cellDataJ.pressure()) / dist;

    GlobalPosition gravityTerm(unitOuterNormal);
    gravityTerm *= (gravity_ * permeability) * density_;

    velocity += gravityTerm;

    //store velocities
    cellDataI.fluxData().setVelocity(isIndexI, velocity);
    cellDataI.fluxData().setVelocityMarker(isIndexI);

    cellDataJ.fluxData().setVelocity(isIndexJ, velocity);
    cellDataJ.fluxData().setVelocityMarker(isIndexJ);
    return;
}

template<class TypeTag>
void FVVelocity1P<TypeTag>::calculateVelocityOnBoundary(const Intersection& intersection, CellData& cellData)
{
    ElementPointer element = intersection.inside();

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
        const GlobalPosition& globalPosI = element->geometry().center();

        // center of face in global coordinates
        const GlobalPosition& globalPosJ = intersection.geometry().center();

        // distance vector between barycenters
        GlobalPosition distVec = globalPosJ - globalPosI;

        // compute distance between cell centers
        Scalar dist = distVec.two_norm();

        //permeability vector at boundary
        // compute vectorized permeabilities
        FieldMatrix meanPermeability(0);

        problem_.spatialParameters().meanK(meanPermeability, problem_.spatialParameters().intrinsicPermeability(*element));

        //multiply with normal vector at the boundary
        Dune::FieldVector < Scalar, dim > permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);
        permeability /= viscosity_;

        Scalar pressBound = boundValues;

        //calculate potential gradients
        Scalar potential = (cellData.pressure() - pressBound) / dist;

        potential += density_ * (unitOuterNormal * gravity_);

        //store potentials for further calculations (velocity, saturation, ...)
        cellData.fluxData().setPotential(isIndex, potential);

        //calculate the gravity term
        GlobalPosition velocity(permeability);
        velocity *= (cellData.pressure() - pressBound) / dist;

        GlobalPosition gravityTerm(unitOuterNormal);
        gravityTerm *= (gravity_ * permeability) * density_;

        velocity += gravityTerm;

        //store velocities
        cellData.fluxData().setVelocity(isIndex, velocity);
        cellData.fluxData().setVelocityMarker(isIndex);

    } //end dirichlet boundary

    else
    {
        problem_.neumann(boundValues, intersection);
        GlobalPosition velocity(unitOuterNormal);

        velocity *= boundValues[pressEqIdx] / density_;

        //store potential gradients for further calculations
        cellData.fluxData().setPotential(isIndex, boundValues[pressEqIdx]);

        //store velocity
        cellData.fluxData().setVelocity(isIndex, velocity);
        cellData.fluxData().setVelocityMarker(isIndex);
    } //end neumann boundary
    return;
}
}
#endif
