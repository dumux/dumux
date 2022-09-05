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
 * \brief Types required for solving a poisson problem using TPFA.
 */

#include <memory>
#include <cmath>

#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/gridvariables.hh>

namespace Dumux {


// Custom problem to define the parameters and BCs of our poisson problem
template<typename GridGeometry>
class PoissonProblem
{
public:
    explicit PoissonProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    {}

    template<typename GlobalPosition>
    double sourceAtPos(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;
        return -sin(globalPos[0])*cos(globalPos[1])
            - sin(globalPos[0])*cos(globalPos[1]);
    }

    template<typename GlobalPosition>
    double exactAtPos(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;
        return sin(globalPos[0])*cos(globalPos[1]);
    }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
};


// Custom grid variables to define a local view
template<typename GridGeometry,
         typename DofVector,
         typename IndexStrategy>
class PoissonGridVariables
: public GridVariables<GridGeometry, DofVector, IndexStrategy>
{
    using ParentType = GridVariables<GridGeometry, DofVector, IndexStrategy>;
    using Problem = PoissonProblem<GridGeometry>;

public:
    class LocalView
    {
    public:
        using GridVariables = PoissonGridVariables;

        explicit LocalView(const GridVariables& gridVars)
        : gridVars_(gridVars)
        {}

        template<typename GGLocalView>
        void bind(const GGLocalView) &
        {}

        template<typename GGLocalView>
        LocalView bind(const GGLocalView) &&
        { return std::move(*this); }

        const GridVariables& gridVariables() const
        { return gridVars_; }

    private:
        const GridVariables& gridVars_;
    };

    template<typename Initializer>
    PoissonGridVariables(std::shared_ptr<const GridGeometry> gridGeometry,
                         std::shared_ptr<const IndexStrategy> indexStrategy,
                         const Initializer& initializer)
    : ParentType(gridGeometry, indexStrategy, initializer)
    , problem_(std::make_shared<const Problem>(gridGeometry))
    {}

    const Problem& problem() const
    { return *problem_; }

private:
    std::shared_ptr<const Problem> problem_;
};


// custom local assembler for our poisson problem
template<typename GG, typename GV>
class PoissonLocalAssembler
{
    using SubControlVolume = typename GG::SubControlVolume;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;

    using LocalGridGeometry = typename GG::LocalView;
    using LocalGridVariables = typename GV::LocalView;

public:
    using GridGeometry = GG;
    using GridVariables = GV;

    PoissonLocalAssembler(const LocalGridGeometry& localGeom,
                          const LocalGridVariables& localVars)
    : localGeom_(localGeom)
    , localVars_(localVars)
    {}

    template<typename Matrix>
    void addJacobianEntries(Matrix& jacobian) const
    { addJacobianEntries_(jacobian); }

    template<typename Vector>
    void addResidualEntries(Vector& residual) const
    { addRHSEntries_(residual); }

    template<typename Matrix, typename Vector>
    void addJacobianAndResidualEntries(Matrix& jacobian, Vector& residual) const
    {
        addJacobianEntries(jacobian);
        addResidualEntries(residual);
    }

private:
    template<typename MultiIndex>
    decltype(auto) getDofIndex_(const MultiIndex& i) const
    { return localVars_.gridVariables().getDofIndex(i); }

    const auto& problem_() const
    { return localVars_.gridVariables().problem(); }

    template<typename Vector>
    void addRHSEntries_(Vector& residual) const
    {
        for (const auto& scv : scvs(localGeom_))
            LinearSystem::add(
                residual,
                getDofIndex_(MultiIndex{scv.dofIndex(), 0}),
                -1.0*problem_().sourceAtPos(scv.dofPosition())*scv.volume()
            );

        for (const auto& scvf : scvfs(localGeom_))
            LinearSystem::add(
                residual,
                getDofIndex_(MultiIndex{scvf.insideScvIdx(), 0}),
                computeFaceFlux_(scvf)
            );
    }

    template<typename Matrix>
    void addJacobianEntries_(Matrix& jacobian) const
    {
        for (const auto& scvf : scvfs(localGeom_))
        {
            const auto deriv = computeFaceTransmissibility_(scvf);
            LinearSystem::add(
                jacobian,
                getDofIndex_(MultiIndex{scvf.insideScvIdx(), 0}),
                getDofIndex_(MultiIndex{scvf.insideScvIdx(), 0}),
                -1.0*deriv
            );

            if (!scvf.boundary())
                LinearSystem::add(
                    jacobian,
                    getDofIndex_(MultiIndex{scvf.insideScvIdx(), 0}),
                    getDofIndex_(MultiIndex{scvf.outsideScvIdx(), 0}),
                    deriv
                );
        }
    }

    double computeFaceFlux_(const SubControlVolumeFace& scvf) const
    {
        const auto tij = computeFaceTransmissibility_(scvf);
        const auto ui = getInsideDofValue_(scvf);
        const auto uj = scvf.boundary() ? problem_().exactAtPos(scvf.ipGlobal())
                                        : getOutsideDofValue_(scvf);
        return tij*(uj - ui);
    }

    double computeFaceTransmissibility_(const SubControlVolumeFace& scvf) const
    {
        const auto& insideScv = localGeom_.scv(scvf.insideScvIdx());
        const auto insideT = Dumux::computeTpfaTransmissibility(scvf, insideScv, 1.0, 1.0)*scvf.area();
        if (scvf.boundary())
            return insideT;

        const auto& outsideScv = localGeom_.scv(scvf.outsideScvIdx());
        const auto outsideT = Dumux::computeTpfaTransmissibility(scvf, outsideScv, 1.0, 1.0)*scvf.area();
        return insideT*outsideT/(outsideT - insideT);
    }

    double getInsideDofValue_(const SubControlVolumeFace& scvf) const
    { return getDofValue_(localGeom_.scv(scvf.insideScvIdx())); }

    double getOutsideDofValue_(const SubControlVolumeFace& scvf) const
    { return getDofValue_(localGeom_.scv(scvf.outsideScvIdx())); }

    double getDofValue_(const SubControlVolume& scv) const
    {
        return LinearSystem::get(
            localVars_.gridVariables().dofs(),
            getDofIndex_(MultiIndex{scv.dofIndex(), 0})
        );
    }

    const LocalGridGeometry& localGeom_;
    const LocalGridVariables& localVars_;
};

} // namespace Dumux
