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
 * \ingroup SequentialOnePModel
 * \brief  Single phase finite volume velocity reconstruction
 */

#ifndef DUMUX_FVVELOCITY1P_HH
#define DUMUX_FVVELOCITY1P_HH

#include <dumux/porousmediumflow/1p/sequential/properties.hh>

namespace Dumux {
/*!
 * \ingroup SequentialOnePModel
 * \brief Single phase finite volume velocity reconstruction
 *
 * Calculates velocities from a known pressure field applying a finite volume discretization.
 * The pressure has to be given as piecewise constant cell values.
 * The velocity is calculated following  Darcy's law as
 * \f[
 * \boldsymbol v = -\frac{1}{\mu} \boldsymbol K \left(\textbf{grad}\, p + \rho g  \textbf{grad}\, z\right),
 * \f]
 * where, \f$ p \f$ is the pressure, \f$ \boldsymbol K \f$ the absolute permeability,
 *        \f$ \mu \f$ the viscosity, \f$ \rho \f$ the density and \f$ g \f$ the gravity constant.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVVelocity1P
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        pressEqIdx = Indices::pressureEqIdx // only one equation!
    };

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    /*!
     * \brief Constructs a FVVelocity1P object
     *
     * \param problem A problem class object
     */
    FVVelocity1P(Problem& problem)
    : problem_(problem), gravity_(problem.gravity())
      {
        const auto element = *problem_.gridView().template begin<0>();
        Scalar temperature = problem_.temperature(element);
        Scalar referencePress = problem_.referencePressure(element);

        density_ = FluidSystem::density(temperature, referencePress);
        viscosity_ = FluidSystem::viscosity(temperature, referencePress);
      }

    // Calculates the velocity at a cell-cell interface.
    void calculateVelocity(const Intersection&, CellData&);

    // Calculates the velocity at a boundary.
    void calculateVelocityOnBoundary(const Intersection&, CellData&);

    /*!
     * \brief Adds velocity output to the output file
     *
     * Adds the velocities to the output.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        Dune::BlockVector<Dune::FieldVector<Scalar, dim> > &velocity = *(writer.template allocateManagedBuffer<Scalar, dim> (
                problem_.gridView().size(0)));

        // compute update vector
        for (const auto& element : elements(problem_.gridView()))
        {
            // cell index
            int eIdxGlobal = problem_.variables().index(element);

            CellData& cellData = problem_.variables().cellData(eIdxGlobal);

            const typename Element::Geometry& geometry = element.geometry();

            // get corresponding reference element
            const auto refElement = referenceElement(geometry);

            const int numberOfFaces = refElement.size(1);
            std::vector<Scalar> flux(numberOfFaces,0);

            // run through all intersections with neighbors and boundary
            for (const auto& intersection : intersections(problem_.gridView(), element))
            {
                int isIndex = intersection.indexInInside();

                flux[isIndex] = intersection.geometry().volume()
                        * (intersection.centerUnitOuterNormal() * cellData.fluxData().velocity(isIndex));
            }

            // calculate velocity on reference element as the Raviart-Thomas-0
            // interpolant of the fluxes
            Dune::FieldVector<Scalar, dim> refVelocity;
            // simplices
            if (refElement.type().isSimplex()) {
                for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                {
                    refVelocity[dimIdx] = -flux[dim - 1 - dimIdx];
                    for (int fIdx = 0; fIdx < dim + 1; fIdx++)
                    {
                        refVelocity[dimIdx] += flux[fIdx]/(dim + 1);
                    }
                }
            }
            // cubes
            else if (refElement.type().isCube()){
                for (int i = 0; i < dim; i++)
                    refVelocity[i] = 0.5 * (flux[2*i + 1] - flux[2*i]);
            }
            // 3D prism and pyramids
            else {
                DUNE_THROW(Dune::NotImplemented, "velocity output for prism/pyramid not implemented");
            }

            const Dune::FieldVector<Scalar, dim>& localPos
              = refElement.position(0, 0);

            // get the transposed Jacobian of the element mapping
            const typename Element::Geometry::JacobianTransposed& jacobianT =
                geometry.jacobianTransposed(localPos);

            // calculate the element velocity by the Piola transformation
            Dune::FieldVector<Scalar, dim> elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= geometry.integrationElement(localPos);

            velocity[eIdxGlobal] = elementVelocity;
        }

        writer.attachCellData(velocity, "velocity", dim);

        return;
    }
private:
    Problem &problem_;
    const GravityVector& gravity_; //!< vector including the gravity constant
    Scalar density_;
    Scalar viscosity_;
};

/*!
 * \ingroup SequentialOnePModel
 * \brief Calculates the velocity at a cell-cell interface.
 *
 * Calculates the velocity at a cell-cell interface from a given pressure field.
 *
 * \param intersection Intersection of two grid cells
 * \param cellData Object containing all model relevant cell data
 */
template<class TypeTag>
void FVVelocity1P<TypeTag>::calculateVelocity(const Intersection& intersection, CellData& cellData)
{
    auto elementI = intersection.inside();
    auto elementJ = intersection.outside();

    int eIdxGlobalJ = problem_.variables().index(elementJ);

    CellData& cellDataJ = problem_.variables().cellData(eIdxGlobalJ);

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = elementI.geometry().center();
    const GlobalPosition& globalPosJ = elementJ.geometry().center();

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
    DimMatrix meanPermeability(0);

    problem_.spatialParams().meanK(meanPermeability, problem_.spatialParams().intrinsicPermeability(elementI),
            problem_.spatialParams().intrinsicPermeability(elementJ));

    Dune::FieldVector < Scalar, dim > permeability(0);
    meanPermeability.mv(unitOuterNormal, permeability);

    permeability /= viscosity_;

    // calculate potential gradients
    Scalar potential = (cellData.pressure() - cellDataJ.pressure()) / dist;

    potential += density_ * (unitOuterNormal * gravity_);

    // store potentials for further calculations (velocity, saturation, ...)
    cellData.fluxData().setPotential(isIndexI, potential);
    cellDataJ.fluxData().setPotential(isIndexJ, -potential);

    // calculate the gravity term
    VelocityVector velocity(permeability);
    velocity *= (cellData.pressure() - cellDataJ.pressure()) / dist;

    GravityVector gravityTerm(unitOuterNormal);
    gravityTerm *= (gravity_ * permeability) * density_;

    velocity += gravityTerm;

    // store velocities
    cellData.fluxData().setVelocity(isIndexI, velocity);
    cellData.fluxData().setVelocityMarker(isIndexI);

    cellDataJ.fluxData().setVelocity(isIndexJ, velocity);
    cellDataJ.fluxData().setVelocityMarker(isIndexJ);
    return;
}

/*!
 * \ingroup SequentialOnePModel
 * \brief Calculates the velocity at a boundary.
 *
 * Calculates the velocity at a boundary from a given pressure field.
 *
 * \param intersection Boundary intersection
 * \param cellData Object containing all model relevant cell data
 */
template<class TypeTag>
void FVVelocity1P<TypeTag>::calculateVelocityOnBoundary(const Intersection& intersection, CellData& cellData)
{
    auto element = intersection.inside();

    // get face index
    int isIndex = intersection.indexInInside();

    // get face normal
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

    BoundaryTypes bcType;
    // get boundary type
    problem_.boundaryTypes(bcType, intersection);
    PrimaryVariables boundValues(0.0);

    if (bcType.isDirichlet(pressEqIdx))
    {
        problem_.dirichlet(boundValues, intersection);

        // get global coordinates of cell centers
        const GlobalPosition& globalPosI = element.geometry().center();

        // center of face in global coordinates
        const GlobalPosition& globalPosJ = intersection.geometry().center();

        // distance vector between barycenters
        GlobalPosition distVec = globalPosJ - globalPosI;

        // compute distance between cell centers
        Scalar dist = distVec.two_norm();

        // permeability vector at boundary
        // compute vectorized permeabilities
        DimMatrix meanPermeability(0);

        problem_.spatialParams().meanK(meanPermeability, problem_.spatialParams().intrinsicPermeability(element));

        // multiply with normal vector at the boundary
        Dune::FieldVector < Scalar, dim > permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);
        permeability /= viscosity_;

        Scalar pressBound = boundValues;

        // calculate potential gradients
        Scalar potential = (cellData.pressure() - pressBound) / dist;

        potential += density_ * (unitOuterNormal * gravity_);

        // store potentials for further calculations (velocity, saturation, ...)
        cellData.fluxData().setPotential(isIndex, potential);

        // calculate the gravity term
        VelocityVector velocity(permeability);
        velocity *= (cellData.pressure() - pressBound) / dist;

        GravityVector gravityTerm(unitOuterNormal);
        gravityTerm *= (gravity_ * permeability) * density_;

        velocity += gravityTerm;

        // store velocities
        cellData.fluxData().setVelocity(isIndex, velocity);
        cellData.fluxData().setVelocityMarker(isIndex);

    } // end Dirichlet boundary

    else
    {
        problem_.neumann(boundValues, intersection);
        VelocityVector velocity(unitOuterNormal);

        velocity *= boundValues[pressEqIdx] / density_;

        // store potential gradients for further calculations
        cellData.fluxData().setPotential(isIndex, boundValues[pressEqIdx]);

        //store velocity
        cellData.fluxData().setVelocity(isIndex, velocity);
        cellData.fluxData().setVelocityMarker(isIndex);
    } // end Neumann boundary
    return;
}
} // end namespace Dumux
#endif
