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
 * \ingroup StokesDropsDarcyCoupling
 * \copydoc Dumux::StokesDropsDarcyDropManager
 */

#ifndef DUMUX_STOKES_DROPS_DARCY_DROPMANAGER_HH
#define DUMUX_STOKES_DROPS_DARCY_DROPMANAGER_HH

#include <dumux/discretization/localview.hh> // TODO

namespace Dumux {

/*!
 * \ingroup StokesDropsDarcyCoupling
 * \brief Drop manager to evaluate, update and store drop information.
 */
template<class MDTraits>
class DropManager
{
    using Scalar = typename MDTraits::Scalar;

public:
    static constexpr auto stokesIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto interfaceIdx = typename MDTraits::template SubDomain<2>::Index();
    static constexpr auto darcyIdx = typename MDTraits::template SubDomain<3>::Index();

private:
    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using InterfaceTypeTag = typename MDTraits::template SubDomain<2>::TypeTag;
    using DarcyTypeTag = typename MDTraits::template SubDomain<3>::TypeTag;

    using InterfaceProblem = GetPropType<InterfaceTypeTag, Properties::Problem>;
    using InterfaceFVGridGeometry = GetPropType<InterfaceTypeTag, Properties::FVGridGeometry>;
    using StokesIndices = typename GetPropType<StokesTypeTag, Properties::ModelTraits>::Indices;
    using DarcyIndices = typename GetPropType<DarcyTypeTag, Properties::ModelTraits>::Indices;

    using SolutionVector = typename MDTraits::SolutionVector;

    struct DropDomainInfo
    {
        std::vector<int> elementIndices;
        Scalar dropVolume;
        bool dropFormation;
        Scalar domainArea; // sum of all elements/scvs in drop domain
        Scalar dropRadius;
        Scalar dropSurfaceArea;
        Scalar dropContactArea;
        Scalar aDrop;
        Scalar aFree;
    };

    struct PoreClass {Scalar meanPoreRadius; Scalar percentage; };

public:
    //! Constructor
    DropManager()
    { }

    //! Initializes the drop manager
    // TODO move to constructor?
    void init(std::shared_ptr<const InterfaceProblem> interfaceProblem)
    {
        // get drop domain information from interface problem
        numberOfDropDomains_ = interfaceProblem->getNumberOfDropDomains();
        dropDomains_.resize(numberOfDropDomains_);
        computeDropDomains(interfaceProblem->fvGridGeometry()); // TODO
        contactAngle_ = interfaceProblem->spatialParams().contactAngle(); // TODO
        surfaceTension_ = interfaceProblem->spatialParams().surfaceTension(); // TODO
        porosity_ = interfaceProblem->spatialParams().porosityInterface(); // TODO heterogeneous porosity
//        poreSizeDistribution_ = interfaceProblem->spatialParams().poreSizeDistribution();
        // TODO work around, needs to be fixed!
        auto poreSizeDistribution = interfaceProblem->spatialParams().poreSizeDistribution();
        for (auto poreClass : poreSizeDistribution)
            poreSizeDistribution_.push_back({poreClass.meanPoreRadius, poreClass.percentage});
    }

    //! Assigns each element to a 'global' drop domain
    // TODO implement for 2D interface domain
    void computeDropDomains(const InterfaceFVGridGeometry& fvGridGeometry)
    {
        const int numberOfElements = fvGridGeometry.numScv();
        const int elementsPerDropDomain = numberOfElements/numberOfDropDomains_;

        int domainCounter = 0, elementCounter = 0;
        for (auto scv : scvs(localView(fvGridGeometry)))
        {
            if(elementCounter >= elementsPerDropDomain)
            {
                elementCounter = 0;
                domainCounter++;
            }
            dropDomains_[domainCounter].elementIndices.push_back(scv.elementIndex());
            elementCounter++;
        }
    }

    //! Checks if a drop forms in the given element, computes drop volume and area fractions
    void evaluateDropFormation(const SolutionVector& sol, const InterfaceFVGridGeometry& fvGridGeometry) const
    {
        for (auto dropDomain : dropDomains_)
        {
            if (dropDomain.dropVolume == 0)
            {
                for (auto elementIdx : dropDomain.elementIndices)
                {
                    const Scalar pressureFF = sol[stokesIdx][elementIdx][StokesIndices::pressureIdx];
                    const Scalar pressurePM = sol[darcyIdx][elementIdx][DarcyIndices::pressureIdx];
                    const Scalar criticalPoreRadius = 2 * surfaceTension_ / (pressurePM - pressureFF);

                    for (auto poreClass : poreSizeDistribution_)
                    {
                        if (poreClass.meanPoreRadius >= criticalPoreRadius)
                        {
                            const Scalar poreArea = 2 * poreClass.meanPoreRadius; // 2D
                            const Scalar interfaceArea = localView(fvGridGeometry).scv(elementIdx).volume();
                            dropDomain.domainArea += interfaceArea;
                            // add contribution to global drop volume (sum of all caps for current pore size)
                            dropDomain.dropVolume += interfaceArea * porosity_ * poreClass.percentage / poreArea;
                        }
                    }
                }
            }
            dropDomain.dropVolume > 0 ? dropDomain.dropFormation = true : dropDomain.dropFormation = false;

            // compute global drop radius
            // TODO 2D (should all be multiplied with width b = 1m)
            dropDomain.dropRadius = std::sqrt(2 * dropDomain.dropVolume / (2 * contactAngle_ - std::sin(contactAngle_)));
            // TODO necessary to store surface and contact area? needed elsewhere??
            dropDomain.dropSurfaceArea = 2 * contactAngle_ * dropDomain.dropRadius;
            dropDomain.dropContactArea = 2 * dropDomain.dropRadius * std::sin(contactAngle_);
            // compute drop-covered and drop-free fractions of the interface
            dropDomain.aDrop = dropDomain.dropContactArea / dropDomain.domainArea;
            dropDomain.aFree = 1.0 - dropDomain.aDrop;
            std::cout << "** dropmanager: drop volume = " << dropDomain.dropVolume << std::endl;
        }
    }

    //! Checks if a drop detaches in the given drop domain
    void checkDropDetachment() const
    {

    }

    //! Updates the drop volume in the given drop domain
    void updateDropVolume()
    {

    }

private:
    int numberOfDropDomains_;
    std::vector<DropDomainInfo> dropDomains_;
    std::vector<PoreClass> poreSizeDistribution_;
    Scalar porosity_;
    Scalar contactAngle_; // TODO obtain from ??
    Scalar surfaceTension_; // TODO obtain from fluidsystem! depends on fluid properties!

};


} // end namespace Dumux

#endif
