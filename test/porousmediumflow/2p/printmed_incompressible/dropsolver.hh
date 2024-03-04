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

#include <dumux/common/exceptions.hh>
#include <dumux/nonlinear/findscalarroot.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/indextraits.hh>

#include <dumux/porenetwork/common/poreproperties.hh>
#include "dropcalculations.hh"
#include "droplet.hh"
#include "dropintersection.hh"

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
    using TimeLoopDrop = TimeLoop<GetPropType<TypeTag, Properties::Scalar>>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ThisType = DropletSolverBase<Imp, TypeTag, IsCoupled>;


public:
    DropletSolverBase(const Problem& problem, const SolutionVector& sol, std::shared_ptr<TimeLoopDrop> timeLoop)
    : problem_(problem)
    , timeLoop_(timeLoop)
    , sol_(sol)
    {
        initialContactAngle_ = getParamFromGroup<Scalar>(problem_.paramGroup(), "Drop.InitialContactAngle", 90); //shoudl be converted to radian
        initialContactAngle_ = M_PI*initialContactAngle_/180;
        surfaceTension_ = getParamFromGroup<Scalar>(problem_.paramGroup(),"SpatialParameters.SurfaceTension", 0.0725);

        previousDropletsNew_ = {};
        droplets_ = {};

        initDroplets_();

        interfaceElementIndex();
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

                if (isAtInterface(element, scv))
                    interfaceElementsIndex_.push_back(scv.elementIndex());
            }
        }
    }

    void update()
    {
        for (const auto& droplet : droplets_)
        {
            updateDropGeometry_(droplet);
        }
    }


    Scalar Pc(const Element& element, const Drop& droplet) const
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

    auto prevDroplet(const Element& element, const SubControlVolume& scv) const
    {
        const auto dofIdxGlobal = scv.dofIndex();
        if (previousDropletsNew_.count(dofIdxGlobal))
            return previousDropletsNew_.at(dofIdxGlobal);

            Drop drop;
            return drop;
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

    bool isAtInterface(const Element& element, const SubControlVolume& scv) const
    {
        auto position = scv.dofPosition();

        return problem().onInlet(position);
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
        auto flux = 0.0;

        const auto& gridGeometry = problem_.gridGeometry();
        auto fvGeometry = localView(gridGeometry);

        for (const auto& elementIdx : droplet.elementIndices())
        {
            const auto& element = gridGeometry.element(elementIdx);
            fvGeometry.bindElement(element);

            flux += asImp_().totalFlux(element, fvGeometry, droplet);
        }

        volume = computeDropVolume_(droplet.volume(), dt, flux);

        contactAngle = computeContactAngle_(volume, contactRadius);

        radius = contactRadius;
        radius /= sin(M_PI - contactAngle);

        height = radius - radius * cos(contactAngle);

        droplet.update(volume,
                       radius,
                       height,
                       contactAngle);
    }

    void initDroplets_()
    {

        Scalar initialContactAngle = 0.0; //Todo
        Scalar initialVolume = 0.0; //Todo
        GlobalPosition initialCenter = 0.0; //Todo
        Scalar initialContactRadius = computeContactRadius_(initialVolume, initialContactAngle);

        Scalar initialRadius = initialContactRadius;
        initialRadius /= sin(M_PI - initialContactAngle);

        Scalar initialHeight = initialRadius - initialRadius * cos(initialContactAngle);


        const auto& gridGeometry = problem_.gridGeometry();
        auto fvGeometry = localView(gridGeometry);


        std::vector<GridIndexType> dropletDoFs;
        std::vector<GlobalPosition> dropletDoFPositions;
        std::vector<GridIndexType> dropletElems;
        for (const auto& elementIdx : interfaceElementsIndex_)
        {
            const auto& element = gridGeometry.element(elementIdx);
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (!isAtInterface(element, scv))
                    continue;

                const auto dofIdxGlobal = scv.dofIndex();
                const auto dofPosition = scv.dofPosition();
                const auto distance = DropIntersection<Scalar, GlobalPosition>::distancePointToPoint(dofPosition, initialCenter);

                if (distance > initialContactRadius)
                    continue;

                if (std::none_of(dropletDoFs.cbegin(), dropletDoFs.cend(), dofIdxGlobal))
                {
                    dropletDoFs.push_back(dofIdxGlobal);
                    dropletElems.push_back(elementIdx);
                    dropletDoFPositions.push_back(dofPosition);
                }

            }
        }

        Drop droplet(initialVolume,
                     initialRadius,
                     initialHeight,
                     initialContactRadius,
                     initialContactAngle,
                     initialCenter,
                     dropletElems,
                     dropletDoFs,
                     dropletDoFPositions);

        droplets_.push_back(droplet);
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
            Scalar upperLimit = M_PI / 2;

        Scalar Theta = findScalarRootBrent(lowerLimit, upperLimit, evalContactAngle, 1e-6, 500);
        return Theta;
    }

    static Scalar computeContactRadius_(const Scalar initialDropVolume, const Scalar initialContactAngle)
    {
        Scalar contactRadius = 3 * initialDropVolume / (M_PI * (1 - cos(initialContactAngle)) * (1 - cos(initialContactAngle)) * (2 + cos(initialContactAngle)) / (sin(initialContactAngle) * sin(initialContactAngle) * sin(initialContactAngle)));
        contactRadius = std::cbrt(stickRadius);

        return contactRadius;
    }


    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this);}

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this);}

    std::size_t numDofs_;
    std::vector<Drop> droplets_;
    mutable Drop droplet_;

    Scalar initialContactAngle_;
    Scalar surfaceTension_;

    const Problem& problem_;
    std::shared_ptr<TimeLoopDrop> timeLoop_;
    const SolutionVector& sol_;
    std::vector<int> interfaceElementsIndex_;
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
    using TimeLoopDrop = TimeLoop<GetPropType<TypeTag, Properties::Scalar>>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Drop = Dumux::PoreNetwork::Droplet<GridView>;


    using ThisType = DropletSolverTwoP<TypeTag, IsCoupled>;
    // using DropDetachment =  DropletDetachment<ThisType, Scalar, IsCoupled>;

public:

    DropletSolverTwoP(const Problem& problem, const GridVariables& gridVariables, const SolutionVector& sol, std::shared_ptr<TimeLoopDrop> timeLoop)
    : ParentType(problem, sol, timeLoop)
    , gridVariables_(gridVariables)
    , sol_{sol}
    {}

    bool update()
    {
        return ParentType::update();
    }

    auto totalFlux(const Element& element, const FVElementGeometry& fvGeometry, const SubControlVolume& scv) const
    {
        using std::isnan;
        Scalar flux(0.0);

        auto elemVolVars = localView(gridVariables_.curGridVolVars());
        elemVolVars.bind(element, fvGeometry, sol_);

        auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());
        elemFluxVarsCache.bindElement(element, fvGeometry, elemVolVars);
        const auto& scvf = fvGeometry.scvf(0);
        flux += flux_(element, fvGeometry, elemVolVars, scv);

        return flux;
    }

    auto totalFlux(const Element& element, const FVElementGeometry& fvGeometry, const Drop& droplet) const
    {
        Scalar flux(0.0);

        const auto dofIdxGlobal = droplet.dofIndex();
        auto&& scv = fvGeometry.scv(droplet.LocalDofIdx());

        return totalFlux(element, fvGeometry, scv);
    }

private:

        //! Evaluates the flux coming from the pore to the droplet
    Scalar flux_(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVolumeVariables& elemVolVars,
                    const SubControlVolume& scv) const
    {
        Scalar flux = 0.0;
        auto phaseIdx = dropletPhaseIdx();
        const auto dofIdxGlobal = scv.dofIndex();
        auto elemFluxVarsCache = localView(this->couplingManager().gridVariables(dropIdx).gridFluxVarsCache());
        elemFluxVarsCache.bindElement(element, fvGeometry, elemVolVars);

        ElementBoundaryTypes<dropIdx> elemBcTypes;
        elemBcTypes.update(this->couplingManager().problem(dropIdx), element, fvGeometry);
        FluxVariables<dropIdx> fluxVars;
        // The flux must be substracted:
        // On an inlet boundary, the flux part of the local residual will be positive, since all fluxes will leave the SCV towards to interior domain.
        // For the domain itself, however, the sign has to be negative, since mass is entering the system.
        const auto& scvf = fvGeometry.scvf(0);
        fluxVars.init(this->couplingManager().problem(dropIdx), element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        auto upwindTerm = [phaseIdx](const auto& volVars) { return volVars.mobility(phaseIdx); };
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

        if(dofIdxGlobal==insideScv.dofIndex())
            flux -= (fluxVars.advectiveFlux(phaseIdx, upwindTerm));
        else if(dofIdxGlobal==outsideScv.dofIndex())
            flux += (fluxVars.advectiveFlux(phaseIdx, upwindTerm));

        return flux;

    }
    const SolutionVector& sol_;
    const GridVariables& gridVariables_;
};

} // namespace Dumux

#endif
