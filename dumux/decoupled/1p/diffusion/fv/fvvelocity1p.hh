// $Id$
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

#include <dumux/decoupled/1p/diffusion/fv/fvpressure1p.hh>

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
class FVVelocity1P: public FVPressure1P<TypeTag>
{
    typedef FVVelocity1P<TypeTag> ThisType;
    typedef FVPressure1P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Fluid)) Fluid;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    //! The Constructor
    /**
     * \param problem a problem class object
     */
    FVVelocity1P(Problem& problem)
    : FVPressure1P<TypeTag>(problem)
      {   }


    //! Calculate the velocity.
    /*!
     *
     *  Given the piecewise constant pressure \f$p\f$,
     *  this method calculates the velocity field
     */
    void calculateVelocity();


    void initialize()
    {
        ParentType::initialize();

        calculateVelocity();

        return;
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);

        typename Variables::ScalarSolutionType *velocityX = writer.template createField<Scalar, 1> (
                this->problem().gridView().size(0));
        typename Variables::ScalarSolutionType *velocityY = writer.template createField<Scalar, 1> (
                this->problem().gridView().size(0));

        // compute update vector
        ElementIterator eItEnd = this->problem().gridView().template end<0>();
        for (ElementIterator eIt = this->problem().gridView().template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // cell index
            int globalIdx = this->problem().variables().index(*eIt);


            Dune::FieldVector<Scalar, 2*dim> flux(0);
            // run through all intersections with neighbors and boundary
            IntersectionIterator
            isItEnd = this->problem().gridView().iend(*eIt);
            for (IntersectionIterator
                    isIt = this->problem().gridView().ibegin(*eIt); isIt
                    !=isItEnd; ++isIt)
            {
                int isIndex = isIt->indexInInside();

                flux[isIndex] = isIt->geometry().volume() * (isIt->centerUnitOuterNormal() * this->problem().variables().velocityElementFace(*eIt, isIndex));
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

            (*velocityX)[globalIdx] = elementVelocity[0];
            (*velocityY)[globalIdx] = elementVelocity[1];
        }

        writer.addCellData(velocityX, "x-velocity");
        writer.addCellData(velocityY, "y-velocity");

        return;
    }

};
template<class TypeTag>
void FVVelocity1P<TypeTag>::calculateVelocity()
{
    // compute update vector
    ElementIterator eItEnd = this->problem().gridView().template end<0>();
    for (ElementIterator eIt = this->problem().gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        //
        GlobalPosition globalPos = eIt->geometry().center();

        // cell index
        int globalIdxI = this->problem().variables().index(*eIt);

        Scalar pressI = this->problem().variables().pressure()[globalIdxI];

        Scalar temperatureI = this->problem().temperature(globalPos, *eIt);
        Scalar referencePressI = this->problem().referencePressure(globalPos, *eIt);

        Scalar densityI = Fluid::density(temperatureI, referencePressI);

        // run through all intersections with neighbors and boundary
        IntersectionIterator
        isItEnd = this->problem().gridView().iend(*eIt);
        for (IntersectionIterator
                isIt = this->problem().gridView().ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // local number of facet
            int isIndex = isIt->indexInInside();

            Dune::FieldVector<Scalar,dimWorld> unitOuterNormal = isIt->centerUnitOuterNormal();

            // get absolute permeability
            FieldMatrix permeabilityI(this->problem().spatialParameters().intrinsicPermeability(globalPos, *eIt));

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->problem().variables().index(*neighborPointer);


                // cell center in global coordinates
                const GlobalPosition& globalPos = eIt->geometry().center();

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                // get absolute permeability
                FieldMatrix permeabilityJ(this->problem().spatialParameters().intrinsicPermeability(globalPosNeighbor, *neighborPointer));

                // compute vectorized permeabilities
                FieldMatrix meanPermeability(0);

//                // harmonic mean of permeability
                for (int x = 0; x < dim; x++)
                {
                    meanPermeability[x][x] = 2 * permeabilityI[x][x] * permeabilityJ[x][x] / (permeabilityI[x][x]
                            + permeabilityJ[x][x]);
                    for (int y = 0; y < dim; y++)
                    {
                        if (x != y)
                        {//use arithmetic mean for the off-diagonal entries to keep the tensor property!
                            meanPermeability[x][y] = 0.5 * (permeabilityI[x][y] + permeabilityJ[x][y]);
                        }
                    }
                }

                Dune::FieldVector<Scalar,dim> permeability(0);
                meanPermeability.mv(unitOuterNormal,permeability);

                Scalar pressJ = this->problem().variables().pressure()[globalIdxJ];

                Scalar temperatureJ = this->problem().temperature(globalPos, *eIt);
                Scalar referencePressJ = this->problem().referencePressure(globalPos, *eIt);

                Scalar densityJ = Fluid::density(temperatureJ, referencePressJ);

                //calculate potential gradients
                Scalar potential = this->problem().variables().potential(globalIdxI, isIndex);

                Scalar density = (potential> 0.) ? densityI : densityJ;

                density= (potential == 0.) ? 0.5 * (densityI + densityJ) : density;


                potential = (pressI - pressJ) / dist;


                potential += density * (unitOuterNormal * this->gravity);

                //store potentials for further calculations (velocity, saturation, ...)
                this->problem().variables().potential(globalIdxI, isIndex) = potential;

                //do the upwinding depending on the potentials
                density = (potential> 0.) ? densityI : densityJ;
                density = (potential == 0.) ? 0.5 * (densityI + densityJ) : density;

                //calculate the gravity term
                Dune::FieldVector<Scalar,dimWorld> velocity(permeability);
                velocity *= (pressI - pressJ)/dist;

                Dune::FieldVector<Scalar,dimWorld> gravityTerm(unitOuterNormal);
                gravityTerm *= (this->gravity*permeability)*density;

                //store velocities
                this->problem().variables().velocity()[globalIdxI][isIndex] = (velocity + gravityTerm);

            }//end intersection with neighbor

            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition globalPosFace = isIt->geometry().center();

                //get boundary type
                BoundaryConditions::Flags bcType = this->problem().bctype(globalPosFace, *isIt);

                // cell center in global coordinates
                GlobalPosition globalPos = eIt->geometry().center();

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPosFace - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                //multiply with normal vector at the boundary
                Dune::FieldVector<Scalar,dim> permeability(0);
                permeabilityI.mv(unitOuterNormal, permeability);

                if (bcType == BoundaryConditions::dirichlet)
                {
                    Scalar pressBound = this->problem().dirichlet(globalPosFace, *isIt);

                    //calculate the gravity term
                    Dune::FieldVector<Scalar,dimWorld> velocity(permeability);
                    velocity *= (pressI - pressBound)/dist;

                    Dune::FieldVector<Scalar,dimWorld> gravityTerm(unitOuterNormal);
                    gravityTerm *= (this->gravity*permeability)*densityI;

                    this->problem().variables().velocity()[globalIdxI][isIndex] = (velocity + gravityTerm);

                  }//end dirichlet boundary

                else
                {
                   Scalar J = this->problem().neumann(globalPosFace, *isIt);
                    Dune::FieldVector<Scalar,dimWorld> velocity(unitOuterNormal);

                    velocity *= J/= densityI;

                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocity;

                }//end neumann boundary
            }//end boundary
        }// end all intersections
    }// end grid traversal
//                        printvector(std::cout, this->problem().variables().velocity(), "velocity", "row", 4, 1, 3);
    return;
}
}
#endif
