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
 *
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters for pore network models.
 */
#ifndef DUMUX_DROP_SOLVER_HH
#define DUMUX_DROP_SOLVER_HH

#include <algorithm>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <memory>

#include <dumux/common/exceptions.hh>
#include <dumux/nonlinear/findscalarroot.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/facetensoraverage.hh>
#include <dumux/flux/upwindscheme.hh>
#include <dumux/porousmediumflow/droplet/timeloopdroplet.hh>

#include "dropintersection.hh"
#include "droplet.hh"

namespace Dumux
{

/**
 * \brief The base class for Drop models.
 */
template<class Imp, class TypeTag, bool IsCoupled>
class DropletSolverBase
{
    using Implementation = Imp;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridView = typename GridGeometry::GridView;
public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Drop = Dumux::Droplet<GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridIndex = typename IndexTraits<GridView>::GridIndex;
private:
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using TimeLoop= TimeLoopDroplet<GetPropType<TypeTag, Properties::Scalar>>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ThisType = DropletSolverBase<Imp, TypeTag, IsCoupled>;


public:
    DropletSolverBase(const Problem& problem)
    : problem_(problem)
    {
        initialContactAngle_ = getParamFromGroup<Scalar>(problem_.paramGroup(), "Drop.InitialContactAngle", 90); //shoudl be converted to radian
        initialContactAngle_ = M_PI * initialContactAngle_ / 180;

        if (hasParamInGroup("", "Drop.InitialContactRadius"))
        {
            initialContactRadius_ = getParamFromGroup<Scalar>(problem_.paramGroup(), "Drop.InitialContactRadius", 1e-3); //TODO
            initialVolume_=  M_PI / 3 * initialContactRadius_ * initialContactRadius_ * initialContactRadius_ * (1 - cos(initialContactAngle_)) * (1 - cos(initialContactAngle_)) * (2 + cos(initialContactAngle_)) / (sin(initialContactAngle_) * sin(initialContactAngle_) * sin(initialContactAngle_));
        }
        else if (hasParamInGroup("", "Drop.InitialVolume"))
        {
            initialVolume_ = getParamFromGroup<Scalar>(problem_.paramGroup(), "Drop.InitialVolume", 1e-3);
            initialContactRadius_ = computeContactRadius_(initialVolume_, initialContactAngle_);
        }

        initialCenter_ = getParamFromGroup<GlobalPosition>(problem_.paramGroup(), "Drop.InitialCenter", GlobalPosition(0.0));
        distanceBetweenDroplets_ = getParamFromGroup<Scalar>(problem_.paramGroup(), "Drop.DistanceBetweenDroplets", 0.0);

        surfaceTension_ = getParamFromGroup<Scalar>(problem_.paramGroup(),"SpatialParameters.SurfaceTension", 0.0725);

        initialRadius_ = initialContactRadius_;
        initialRadius_ /= sin(initialContactAngle_);

        initialHeight_ = initialRadius_ - initialRadius_ * cos(initialContactAngle_);


        interfaceElementIndex();

        initDroplet_(); // TODO
    }

    void interfaceElementIndex()
    {
        const auto& gridGeometry = problem_.gridGeometry();
        auto fvGeometry = localView(gridGeometry);

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (isAtInterface(element, scv)) //Todo
               {
                    interfaceScvs_.push_back(scv);
                    interfaceElementsIndex_.push_back(scv.elementIndex());
               }
            }
        }
    }

    void update()
    {
        for (auto& droplet : droplets_)
        {
            updateDropGeometry_(droplet);
        }

    }

    void dispenseDroplet()
    {
        if (ifDispenseDroplet_())
            initDroplet_();
    }

    bool checkDropletsGeometry() const
    {
        for (const auto& droplet : droplets_)
        {
            if (!checkDropletGeometry_(droplet))
                return false;
        }

        return true;
    }


    Scalar Pc(const Element& element, const Drop& droplet) const
    {
        const auto radius = droplet.radius();
        if (radius == 0.0)
            return 0.0;

        return 2 * surfaceTension_ / radius;
    }

    Scalar Pc(const Drop& droplet) const
    {
        const auto radius = droplet.radius();
        if (radius == 0.0)
            return 0.0;

        return 2 * surfaceTension_ / radius;
    }

    Scalar dropRadius(const Element& element, const Drop& droplet) const
    {
        return droplet.radius();
    }


    Scalar dropContactAngle(const Element& element, const Drop& droplet) const
    {
        return 180* droplet.contactAngle() / M_PI;
    }

    Scalar dropVolume(const Element& element, const Drop& droplet) const
    {
        return droplet.volume();
    }

    Scalar dropStickRadius(const Element& element, const Drop& droplet) const
    {
        return droplet.stickRadius();
    }

    Scalar dropHeight(const Element& element, const Drop& droplet) const
    {
        return droplet.height();
    }

    auto dropCenter(const Element& element, const Drop& droplet) const
    {
        return droplet.center();
    }


    const Problem& problem() const
    {
        return problem_;
    }

    auto droplets() const
    {
        return droplets_;
    }

    auto droplet(const Element& element, const SubControlVolume& scv) const
    {
        auto dofIdxGlobal = scv.dofIndex();
        return droplets_.at(dofIdxGlobal);
    }

    auto droplet() const
    {
        return droplet_;
    }

    auto timeLoop()
    { return timeLoop_;}

    Scalar surfaceTension() const
    {
        return surfaceTension_;
    }

    bool hasDroplet(const Element& element, const SubControlVolume& scv) const
    {
        if (isAtInterface(element, scv))
            return droplets_.count(scv.dofIndex());
        return false;
    }

    bool isCoupledWithDroplet(const GlobalPosition &globalPos) const
    {
        if (problem().onUpperBoundary(globalPos))
        {
            for (const Drop& droplet : droplets_)
            {
                const auto& dropletVolume = droplet.volume();
                if (dropletVolume < 0.0)
                    continue;
                const auto& dropletCenter = droplet.center();
                const auto& dropletContactRadius = droplet.contactRadius();

                const auto distance = DropIntersection<Scalar, GlobalPosition>::distancePointToPoint(globalPos, dropletCenter);
                if (distance < dropletContactRadius)
                {
                    droplet_ = droplet;
                    return true;
                }
            }
        }
        return false;
    }

    bool isAtInterface(const Element& element, const SubControlVolume& scv) const //Todo
    {
        auto position = scv.dofPosition();
        return problem().onUpperBoundary(position);
    }

    void setTimeLoop(std::shared_ptr<TimeLoop> timeLoop)
    { timeLoop_ = timeLoop; }

    Scalar suggestTimeStepSize()
    {
        Scalar suggestedTimeStepSize = timeLoop_->maxTimeStepSize();

        for (auto& droplet : droplets_)
        {
            Scalar flux = -volumeFlux_(droplet);
            suggestedTimeStepSize = std::min(suggestedTimeStepSize, droplet.volume()/flux);
        }

        return suggestedTimeStepSize;
    }

private:

    void updateDropGeometry_(Drop& droplet)
    {

        Scalar radius = 0.0;
        Scalar contactAngle = 0.0;
        Scalar volume = 0.0;
        Scalar height = 0.0;
        Scalar contactRadius = droplet.contactRadius();
        const Scalar dt = timeLoop_->timeStepSize();
        Scalar flux = volumeFlux_(droplet);

        volume = computeDropVolume_(droplet.volume(), dt, flux);
        std::cout<<" ------------volume------------- "<<volume<<std::endl;
        if (volume < 1e-6*droplet.initialVolume())
        {
            auto it = std::find(droplets_.begin(), droplets_.end(), droplet);
            droplets_.erase(it);
            return;
        }
        contactAngle = computeContactAngle_(volume, contactRadius);

        radius = contactRadius;
        radius /= sin(contactAngle);

        height = radius - radius * cos(contactAngle);

        droplet.update(volume,
                       radius,
                       height,
                       contactAngle);
    }

    bool checkDropletGeometry_(const Drop& droplet) const
    {

        Scalar contactRadius = droplet.contactRadius();
        const Scalar dt = timeLoop_->timeStepSize();
        Scalar flux = volumeFlux_(droplet);

        Scalar volume = computeDropVolume_(droplet.volume(), dt, flux);

        if (volume < - (1e-6 * droplet.initialVolume()))
        {
            std::cout<<" ------------Checkvolume------------- "<<volume<<std::endl;
            return false;
        }

        return true;
    }

    void initDroplet_()
    {
        dropletCounter_++;

        auto dropletCenter = initialCenter_;
        dropletCenter[0] += (dropletCounter_ - 1) * distanceBetweenDroplets_;  //TODO for now it just change the position in x-direction. Extend it!

        Drop droplet(initialVolume_,
                     initialRadius_,
                     initialHeight_,
                     initialContactRadius_,
                     initialContactAngle_,
                     dropletCenter);

std::cout<<" ------------initialVolume------------- "<<initialVolume_<<std::endl;
        std::vector<GridIndex> dropletDoFs;
        std::vector<GlobalPosition> dropletDoFPositions;
        std::vector<GridIndex> dropletElems;

        for (const auto& scv : interfaceScvs_)
        {

            const auto dofIdxGlobal = scv.dofIndex();
            const auto dofPosition = scv.dofPosition();
            const auto elementIdx = scv.elementIndex();
            const auto distance = DropIntersection<Scalar, GlobalPosition>::distancePointToPoint(dofPosition, dropletCenter);

            if ( !problem().onUpperBoundary(dofPosition) || distance > initialContactRadius_ )
                continue;

            if (std::none_of(dropletDoFs.cbegin(), dropletDoFs.cend(), [&](GridIndex i){return i == dofIdxGlobal;}))
            {
                dropletDoFs.push_back(dofIdxGlobal);
                dropletElems.push_back(elementIdx);
                dropletDoFPositions.push_back(dofPosition);
            }

        }

        for (const auto dropdof : dropletDoFs)
        std::cout<<"   dropletDoFs  "<<dropdof<<std::endl;

        droplet.couplingData( dropletElems,
                              dropletDoFs,
                              dropletDoFPositions);

       droplets_.push_back(droplet);
    }

    Scalar volumeFlux_(const Drop& droplet) const
    {
        Scalar flux = 0.0;
        const auto& gridGeometry = problem_.gridGeometry();
        auto fvGeometry = localView(gridGeometry);

        for (const auto& elementIdx : droplet.elementIndices())
        {
            const auto& element = gridGeometry.element(elementIdx);
            fvGeometry.bindElement(element);
            flux += asImp_().totalVolumeFlux(element, fvGeometry);
        }

        return flux;
    }

     static Scalar computeDropVolume_(Scalar prevDropVolume, Scalar dt, Scalar flux)
    {
        Scalar curDropVol = prevDropVolume;
        curDropVol += (dt * flux);

        return curDropVol;
    }

    static Scalar computeContactAngle_(const Scalar dropVolume, const Scalar contactRadius) //Todo
    {
        auto evalContactAngle = [&](Scalar contactAngle)
        {
            Scalar res = (sin(contactAngle) * sin(contactAngle) * sin(contactAngle)) * dropVolume;
            res -= M_PI / 3 * contactRadius * contactRadius * contactRadius * (1 - cos(contactAngle)) * (1 - cos(contactAngle)) * (2 + cos(contactAngle));
            return res;
        };

        if (dropVolume < 1e-30)
            return 0.0;

        Scalar lowerLimit = 0;
        Scalar upperLimit = M_PI / 2 + 1;

        Scalar Theta = findScalarRootBrent(lowerLimit, upperLimit, evalContactAngle, 1e-6, 500);
        return Theta;
    }

    static Scalar computeContactRadius_(const Scalar initialDropVolume, const Scalar initialContactAngle)
    {
        Scalar contactRadius = 3 * initialDropVolume / (M_PI * (1 - cos(initialContactAngle)) * (1 - cos(initialContactAngle)) * (2 + cos(initialContactAngle)) / (sin(initialContactAngle) * sin(initialContactAngle) * sin(initialContactAngle)));
        contactRadius = std::cbrt(contactRadius);

        return contactRadius;
    }

    bool ifDispenseDroplet_()
    {
        Scalar temp = timeLoop_->time() / timeLoop_->dropletDispenseTimeInterval();

        return (std::abs(std::round(temp) - temp) / temp < 1e-10);
    }


    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this);}

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this);}

    std::size_t numDofs_;
    std::vector<Drop> droplets_{};
    mutable Drop droplet_;

    Scalar initialContactAngle_;
    Scalar initialVolume_;
    Scalar initialContactRadius_;
    Scalar initialRadius_;
    Scalar initialHeight_;
    GlobalPosition initialCenter_;
    Scalar distanceBetweenDroplets_;

    Scalar surfaceTension_;

    const Problem& problem_;
    std::shared_ptr<TimeLoop> timeLoop_;
    std::vector<int> interfaceElementsIndex_;
    std::vector<SubControlVolume> interfaceScvs_;
    int dropletCounter_ = 0;
};



template<class TypeTag, bool IsCoupled>
class DropletSolverTwoP : public DropletSolverBase <DropletSolverTwoP <TypeTag, IsCoupled>, TypeTag, IsCoupled>
{
    using ParentType = DropletSolverBase <DropletSolverTwoP<TypeTag, IsCoupled>, TypeTag, IsCoupled>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Drop = Dumux::Droplet<GridView>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Extrusion = Extrusion_t<GridGeometry>;
    using UpScheme = UpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>>;


    using ThisType = DropletSolverTwoP<TypeTag, IsCoupled>;
    // using DropDetachment =  DropletDetachment<ThisType, Scalar, IsCoupled>;

public:

    DropletSolverTwoP(const Problem& problem, const SolutionVector& sol)
    : ParentType(problem)
    , sol_(sol)
    {}


    auto totalVolumeFlux(const Element& element, const FVElementGeometry& fvGeometry) const
    {
        Scalar flux(0.0);

        flux += volumeFlux_(element, fvGeometry);

        return flux;
    }

    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    Scalar dropletMassFlux(const Element& element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf,
                      const ElementFluxVarsCache& elemFluxVarCache) const
    {
        auto phaseIdx = 0; //dropletPhaseIdx();

        Scalar fluxBeforeUpwinding = flux_(this->problem(), element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarCache);



        return fluxBeforeUpwinding * 1e6; // todo replace 1e6 with mobility * density
        // auto upwindTerm = [phaseIdx](const auto& volVars) { return volVars.density(phaseIdx) * volVars.mobility(phaseIdx); };

        // return UpScheme::apply(elemVolVars, scvf, upwindTerm, fluxBeforeUpwinding, phaseIdx);
    }

    void setSolutionVector(const SolutionVector& sol)
    { sol_ = sol; }

    void setGridVariables(std::shared_ptr<GridVariables> gridVariables)
    { gridVariables_ = gridVariables; }


private:

    //! Evaluates the flux coming from the pore to the droplet  //TODO
    Scalar volumeFlux_(const Element& element,
                const FVElementGeometry& fvGeometry) const
    {
        Scalar flux = 0.0;
        auto phaseIdx = 0; //dropletPhaseIdx();

        auto elemVolVars = localView(gridVariables_->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, sol_);
        auto elemFluxVarsCache = localView(gridVariables_->gridFluxVarsCache());
        elemFluxVarsCache.bindElement(element, fvGeometry, elemVolVars);

        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(this->problem(), element, fvGeometry);
        FluxVariables fluxVars;
        // The flux must be substracted:
        // On an inlet boundary, the flux part of the local residual will be positive, since all fluxes will leave the SCV towards to interior domain.
        // For the domain itself, however, the sign has to be negative, since mass is entering the system.
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary())
                continue;
            Scalar fluxBeforeUpwinding = flux_(this->problem(), element, fvGeometry, elemVolVars, scvf, 0, elemFluxVarsCache);

            // fluxVars.init(this->problem(), element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

            auto upwindTerm = [phaseIdx](const auto& volVars) { return volVars.mobility(phaseIdx); };
            // const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            // const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
            // // std::cout<<" ------- insideScv.dofIndex() ------"<<insideScv.dofIndex()<<std::endl;
            // // std::cout<<" ------- outsideScv.dofIndex() ------"<<outsideScv.dofIndex()<<std::endl;
            // // std::cout<<" ------- scvf.center() ------"<<scvf.center()<<std::endl;
            // // std::cout<<" ------- scvf.area() ------"<<scvf.area()<<std::endl;
            // // std::cout<<" ----------------- flux ---------"<<flux<<std::endl;

            flux += UpScheme::apply(elemVolVars, scvf, upwindTerm, fluxBeforeUpwinding, phaseIdx);
            // // std::cout<<" ----------------- flux ---------"<<flux<<std::endl;
        }

        return flux;

    }

    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    Scalar flux_(const Problem& problem,
                 const Element& element,
                 const FVElementGeometry& fvGeometry,
                 const ElementVolumeVariables& elemVolVars,
                 const SubControlVolumeFace& scvf,
                 const int phaseIdx,
                 const ElementFluxVarsCache& elemFluxVarCache) const
    {
        const auto& fluxVarCache = elemFluxVarCache[scvf];
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        auto insideK = insideVolVars.permeability();
        auto outsideK = outsideVolVars.permeability();

        // scale with correct extrusion factor
        insideK *= insideVolVars.extrusionFactor();
        outsideK *= outsideVolVars.extrusionFactor();

        const auto K = faceTensorAverage(insideK, outsideK, scvf.unitOuterNormal());
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        const auto& shapeValues = fluxVarCache.shapeValues();

        // evaluate gradP - rho*g at integration point
        Dune::FieldVector<Scalar, dimWorld> gradP(0.0);
        Scalar rho(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            if (enableGravity)
                rho += volVars.density(phaseIdx)*shapeValues[scv.indexInElement()][0];

            Scalar pressure = volVars.pressure(phaseIdx);
            if (this->isCoupledWithDroplet(scv.dofPosition()))
            {
                const auto& droplet = this->droplet();
                pressure = 1e5 + this->Pc(droplet);
            }
            // the global shape function gradient
            gradP.axpy(pressure, fluxVarCache.gradN(scv.indexInElement()));
            // std::cout<<"-----scv.indexInElement()---"<<scv.indexInElement()<<"  ---------   "<<scv.dofIndex()<<"  ---  "<<scv.dofPosition()<<"  ---- "<<pressure<<std::endl;
            // std::cout<<"-----fluxVarCache.gradN(scv.indexInElement())---"<<fluxVarCache.gradN(scv.indexInElement())<<std::endl;

        }

        // std::cout<<" ------------gradP------------- "<<gradP<<std::endl;
        // std::cout<<" ***************************************************************************************** "<<gradP<<std::endl;

        if (enableGravity)
            gradP.axpy(-rho, problem.spatialParams().gravity(scvf.center()));

        // apply the permeability and return the flux
        return -1.0*vtmv(scvf.unitOuterNormal(), K, gradP)*Extrusion::area(fvGeometry, scvf);
    }
    const SolutionVector& sol_;
    std::shared_ptr<GridVariables> gridVariables_;
};

} // namespace Dumux

#endif
