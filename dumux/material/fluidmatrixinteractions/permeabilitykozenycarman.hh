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
 * \brief The Kozeny-Carman relationship for the calculation of a porosity-dependent permeability.
 */
#ifndef DUMUX_PERMEABILITY_KOZENY_CARMAN_HH
#define DUMUX_PERMEABILITY_KOZENY_CARMAN_HH

namespace Dumux
{

/*!
 * \ingroup fluidmatrixinteractionslaws
 */

/**
 * \brief The Kozeny-Carman relationship for the calculation of a porosity-dependent permeability.
 *        When the porosity is implemented as solution-independent, using this relationship for the
 *        permeability leads to unnecessary overhead.
 */
template<class TypeTag, class PermType>
class PermeabilityKozenyCarman
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using InitPoroField = std::function<Scalar(const Element&, const SubControlVolume&)>;
    using InitPermField = std::function<PermType(const Element&, const SubControlVolume&)>;

public:
    PermeabilityKozenyCarman(const Problem& problem) : problemPtr_(&problem) {}

    // the initial parameter distribution
    void init(InitPoroField&& poro,
              InitPermField&& perm)
    {
        initPoro_ = poro;
        initPerm_ = perm;
    }

    // calculate permeability for a given scv
    PermType evaluatePermeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        auto initPoro = initPoro_(element, scv);
        auto poro = problem_().spatialParams().porosity(element, scv, elemSol);
        auto factor = std::pow((1.0 - initPoro)/(1.0 - poro), 2) * std::pow(poro/initPoro, 3);
        return applyFactorToPermeability_(initPerm_(element, scv), factor);
    }

private:
    const Problem& problem_() const
    {return *problemPtr_; }

    Scalar applyFactorToPermeability_(Scalar k, Scalar factor) const
    { return k*factor; }

    Tensor applyFactorToPermeability_(Tensor K, Scalar factor) const
    {
        Tensor result(K);
        for (int i = 0; i < dim; ++i)
            result[i][i] *= factor;
        return result;
    }

    const Problem* problemPtr_;
    InitPoroField initPoro_;
    InitPermField initPerm_;
};

} // namespace Dumux

#endif
