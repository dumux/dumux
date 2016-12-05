// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Base class for all fully implicit problems
 */
#ifndef DUMUX_IMPLICIT_STAGGERED_PROBLEM_HH
#define DUMUX_IMPLICIT_STAGGERED_PROBLEM_HH


#include <dumux/implicit/problem.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitBaseProblems
 * \brief Base class for all fully implicit problems
 *
 * \note All quantities are specified assuming a threedimensional
 *       world. Problems discretized using 2D grids are assumed to be
 *       extruded by \f$1 m\f$ and 1D grids are assumed to have a
 *       cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class StaggeredProblem : public ImplicitProblem<TypeTag>
{
private:
    using ParentType = ImplicitProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using VtkMultiWriter = typename GET_PROP_TYPE(TypeTag, VtkMultiWriter);
    using NewtonMethod = typename GET_PROP_TYPE(TypeTag, NewtonMethod);
    using NewtonController = typename GET_PROP_TYPE(TypeTag, NewtonController);
    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using PointSourceHelper = typename GET_PROP_TYPE(TypeTag, PointSourceHelper);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);


    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using CoordScalar = typename GridView::Grid::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    enum { adaptiveGrid = GET_PROP_VALUE(TypeTag, AdaptiveGrid) };

    using GridAdaptModel = ImplicitGridAdapt<TypeTag, adaptiveGrid>;
    using BoundingBoxTree = Dumux::BoundingBoxTree<GridView>;

    // copying a problem is not a good idea
    StaggeredProblem(const StaggeredProblem &);

public:

    StaggeredProblem(TimeManager &timeManager, const GridView &gridView) : ParentType(timeManager, gridView)
    {}

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param scvFace the sub control volume face
     *
     * The method returns the boundary types information.
     */
    CellCenterPrimaryVariables ccDirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        // forward it to the method which only takes the global coordinate
        if (isBox)
        {
            DUNE_THROW(Dune::InvalidStateException, "dirichlet(scvf) called for box method.");
        }
        else
            return asImp_().dirichletAtPos(scvf.center());
    }

    CellCenterPrimaryVariables ccDirichlet(const Element &element, const SubControlVolume &scv) const
    {
        // forward it to the method which only takes the global coordinate
        if (!isBox)
        {
            DUNE_THROW(Dune::InvalidStateException, "dirichlet(scv) called for cell-centered method.");
        }
        else
            return asImp_().dirichletAtPos(scv.dofPosition());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the dirichlet condition ought to be
     *            set in global coordinates
     *
     * For this method, the \a values parameter stores primary variables.
     */
    CellCenterPrimaryVariables ccDirichletAtPos(const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are dirichlet, but does not provide "
                   "a dirichlet() method.");
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the dirichlet condition ought to be
     *            set in global coordinates
     *
     * For this method, the \a values parameter stores primary variables.
     */
    FacePrimaryVariables faceDirichletAtPos(const GlobalPosition &globalPos, const int direction) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are dirichlet, but does not provide "
                   "a dirichlet() method.");
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local subcontrolvolume index
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    CellCenterPrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolvars,
                             const SubControlVolumeFace& scvFace) const
    {
        // forward it to the interface with only the global position
        return asImp_().neumannAtPos(scvFace.center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param globalPos The position of the boundary face's integration point in global coordinates
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    CellCenterPrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Neumann conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are neumann, but does not provide "
                   "a neumannAtPos() method.");
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The subcontrolvolume
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    CellCenterPrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        // forward to solution independent, fully-implicit specific interface
        return asImp_().sourceAtPos(scv.dofPosition());
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param globalPos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    CellCenterPrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a sourceAtPos() method.");
    }


    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position of the center of the finite volume
     *            for which the initial values ought to be
     *            set (in global coordinates)
     *
     * For this method, the \a values parameter stores primary variables.
     */
    CellCenterPrimaryVariables initialCCValuesAtPos(const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for initial values)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a initialAtPos() method.");
    }


    /*!
     * \brief Evaluate the initial value for a facet.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position of the center of the finite volume
     *            for which the initial values ought to be
     *            set (in global coordinates)
     *
     * For this method, the \a values parameter stores primary variables.
     */
    FacePrimaryVariables initialFaceValueAtPos(const GlobalPosition &globalPos, const int direction) const
    {
        // Throw an exception (there is no reasonable default value
        // for initial values)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a initialAtPos() method.");
    }

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }


};
} // namespace Dumux

#include <dumux/implicit/adaptive/gridadaptpropertydefaults.hh>

#endif
