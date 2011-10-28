/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_HH

#include "MpNcfluxvariables.hh"
#include "diffusion/diffusion.hh"
#include "energy/MpNclocalresidualenergy.hh"
#include "mass/MpNclocalresidualmass.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dumux/common/math.hh>


namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \ingroup BoxLocalResidual
 * \brief two-phase, N-component specific details needed to
 *        approximately calculate the local defect in the BOX scheme.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the
 * two-phase, N-component twophase flow.
 */
template<class TypeTag>
class MPNCLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;

protected:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) Implementation;
    typedef MPNCLocalResidual<TypeTag> ThisType;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        enableEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)),
        enableKineticEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableKineticEnergy)),

        enableDiffusion = GET_PROP_VALUE(TypeTag, PTAG(EnableDiffusion)),
        enableKinetic = GET_PROP_VALUE(TypeTag, PTAG(EnableKinetic)),
        enableSmoothUpwinding = GET_PROP_VALUE(TypeTag, PTAG(EnableSmoothUpwinding)),

        phase0NcpIdx = Indices::phase0NcpIdx
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GridView::CollectiveCommunication CollectiveCommunication;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    typedef MPNCLocalResidualEnergy<TypeTag, enableEnergy, enableKineticEnergy> EnergyResid;
    typedef MPNCLocalResidualMass<TypeTag, enableKinetic> MassResid;

public:
    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(PrimaryVariables &storage, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        storage =0;

        // compute mass and energy storage terms
        MassResid::computeStorage(storage, volVars);
        Valgrind::CheckDefined(storage);
        EnergyResid::computeStorage(storage, volVars);
        Valgrind::CheckDefined(storage);
    }

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within all sub-control volumes of an
     *        element.
     */
    void addPhaseStorage(PrimaryVariables &storage,
                         const Element &element,
                         int phaseIdx) const
    {
        // create a finite volume element geometry
        FVElementGeometry fvElemGeom;
        fvElemGeom.update(this->gridView_(), element);

        // calculate volume variables
        ElementVolumeVariables volVars;
        this->model_().setHints(element, volVars);
        volVars.update(this->problem_(),
                       element,
                       fvElemGeom,
                       /*useOldSolution=*/false);

        // calculate the phase storage for all sub-control volumes
        for (int scvIdx=0;
             scvIdx < fvElemGeom.numVertices;
             scvIdx++)
        {
            PrimaryVariables tmp(0.0);

            // compute mass and energy storage terms in terms of
            // averaged quantities
            MassResid::addPhaseStorage(tmp,
                                       volVars[scvIdx],
                                       phaseIdx);
            EnergyResid::addPhaseStorage(tmp,
                                         volVars[scvIdx],
                                         phaseIdx);

            // multiply with volume of sub-control volume
            tmp *= fvElemGeom.subContVol[scvIdx].volume;

            // Add the storage of the current SCV to the total storage
            storage += tmp;
        }
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(PrimaryVariables &source,
                       int scvIdx)
     {
        Valgrind::SetUndefined(source);
        this->problem_().boxSDSource(source,
                                     this->elem_(),
                                     this->fvElemGeom_(),
                                     scvIdx,
                                     this->curVolVars_() );
        const VolumeVariables &volVars = this->curVolVars_(scvIdx);

        PrimaryVariables tmp(0);
        MassResid::computeSource(tmp, volVars);
        source += tmp;
        Valgrind::CheckDefined(source);

        /*
         *      EnergyResid also called in the MassResid
         *      1) Makes some sense because energy is also carried by mass
         *      2) The mass transfer between the phases is needed.
         */
//        tmp = 0.;
//        EnergyResid::computeSource(tmp, volVars);
//        source += tmp;
//        Valgrind::CheckDefined(source);
     };


    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     */
    void computeFlux(PrimaryVariables &flux, int faceIdx) const
    {
        FluxVariables fluxVars(this->problem_(),
                               this->elem_(),
                               this->fvElemGeom_(),
                               faceIdx,
                               this->curVolVars_());

        flux = 0.0;
        MassResid::computeFlux(flux, fluxVars, this->curVolVars_() );
        Valgrind::CheckDefined(flux);
/*
 *      EnergyResid also called in the MassResid
 *      1) Makes some sense because energy is also carried by mass
 *      2) The component-wise mass flux in each phase is needed.
 */
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     */
    void eval(const Element &element)
    { ParentType::eval(element); }

    /*!
     * \brief Evaluate the local residual.
     */
    void eval(const Element &element,
              const FVElementGeometry &fvGeom,
              const ElementVolumeVariables &prevVolVars,
              const ElementVolumeVariables &curVolVars,
              const ElementBoundaryTypes &bcType)
    {
        ParentType::eval(element,
                         fvGeom,
                         prevVolVars,
                         curVolVars,
                         bcType);

        for (int i = 0; i < this->fvElemGeom_().numVertices; ++i) {
            // add the two auxiliary equations, make sure that the
            // dirichlet boundary condition is conserved
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                if (!bcType[i].isDirichlet(phase0NcpIdx + phaseIdx))
                {
                    this->residual_[i][phase0NcpIdx + phaseIdx] =
                        this->curVolVars_(i).phaseNcp(phaseIdx);
                }
            }
        }
    }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namepace

#endif
