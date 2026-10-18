// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Contact
 * \brief Frictionless penalty contact between the parts of a contact surface
 */
#ifndef DUMUX_SOLIDMECHANICS_CONTACT_PENALTY_CONTACT_HH
#define DUMUX_SOLIDMECHANICS_CONTACT_PENALTY_CONTACT_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/solidmechanics/contact/closestpoint.hh>
#include <dumux/solidmechanics/contact/penalty.hh>
#include <dumux/solidmechanics/contact/search.hh>

namespace Dumux::Contact {

/*!
 * \ingroup Contact
 * \brief Frictionless penalty contact between the parts of a contact surface (self-contact and
 *        contact between bodies of one mesh)
 *
 * Every surface node \f$ s \f$ is tested against its candidate facets \f$ f \f$, the facets
 * opposite to it. For each, \f$ g_f = (\mathbf{x}_s - \mathbf{c}_0)\cdot\mathbf{n}_f \f$ is the
 * signed distance to the plane of the facet with the unit outward normal \f$ \mathbf{n}_f \f$,
 * negative if the node is behind the facet, \f$ d_f \f$ the distance to the facet itself with the
 * closest point \f$ \mathbf{y}_f = \sum_j N_j \mathbf{c}_j \f$, and
 * \f$ q_f = d_f^2 - g_f^2 \f$ the squared distance of the foot point
 * \f$ \mathbf{x}_s - g_f \mathbf{n}_f = \sum_j \xi_j \mathbf{c}_j \f$ from the facet. The energy
 * of the node is the weighted mean of the facet energies,
 * \f[
 *     E_s = \frac{1}{2} A_s \sum_f w_f \lambda(q_f) \psi(g_f), \quad
 *     w_f = \frac{\exp(-\tilde{d}_f/\varepsilon)}{\sum_k \exp(-\tilde{d}_k/\varepsilon)}, \quad
 *     \tilde{d}_f = \sqrt{d_f^2 + (\varepsilon/10)^2},
 * \f]
 * with the reference dual area \f$ A_s \f$, the energy density \f$ \psi \f$ of the penalty law,
 * the smoothing length \f$ \varepsilon \f$ and the lateral cutoff
 * \f$ \lambda(q) = (1 - q/L^2)^2 \f$ for \f$ q < L^2 \f$, zero beyond, with the search radius
 * \f$ L \f$. The nearest facet dominates where it is nearer than the others by several
 * \f$ \varepsilon \f$, facets about equally near are blended, and a facet whose foot point lies
 * outside of it fades out over the distance \f$ L \f$. All parts are continuously
 * differentiable, so the contact force changes continuously, also across concave corners and
 * at the edges of the opposite surface. The factor 1/2 accounts for every pair of opposite
 * surface parts being visited twice, once from each side.
 *
 * The derivatives follow from \f$ \partial g_f/\partial\mathbf{x}_s = \mathbf{n}_f \f$,
 * \f$ \partial g_f/\partial\mathbf{c}_j = -\xi_j \mathbf{n}_f \f$,
 * \f$ \partial d_f^2/\partial\mathbf{x}_s = 2(\mathbf{x}_s - \mathbf{y}_f) \f$ and
 * \f$ \partial d_f^2/\partial\mathbf{c}_j = -2N_j(\mathbf{x}_s - \mathbf{y}_f) \f$, which hold for
 * a plane and for the closest point on a convex set. Each facet contributes forces with zero sum
 * and zero moment. The Jacobian is computed by central differences of the residual.
 *
 * The candidate facets of every node are searched once for a given configuration
 * (updateCandidates) within the search radius \f$ \delta + m \f$ (offset of the law plus search
 * margin) and are kept fixed while the forces follow the displacement; a step from that
 * configuration is covered by them if no surface node moves by more than \f$ m/2 \f$
 * (candidatesCover). A facet is a candidate for a node if it
 * - has no corner within ExclusionRings - 1 surface edges of the node (one ring excludes the
 *   facets of the node itself),
 * - faces the node: the node pseudo-normal and the facet normal point against each other,
 * - has the node in front of it: \f$ |g_f| \geq \f$ MinAlignment \f$ d_f \f$, which excludes
 *   facets beside the node, such as the facets of its own surface around a tight fold,
 * - is penetrated by the node by at most MaxPenetration.
 *
 * Parameters (in the group Contact): PenaltyStiffness, Offset, SearchMargin (default: Offset),
 * SmoothingLength (default: Offset/10), ExclusionRings (default 2), MinAlignment (default 0.5),
 * MaxPenetration (default: Offset). Parts of the solid thinner than the search radius are not
 * resolved.
 */
template<class Surface, class Law = PenaltyLaw<typename Surface::GlobalPosition::value_type>>
class PenaltyContact
{
    using GlobalPosition = typename Surface::GlobalPosition;
    using Scalar = typename GlobalPosition::value_type;
    static constexpr int dim = Surface::dimension;
    using Normals = typename Surface::Normals;
    using Corners = typename Surface::Corners;

    // the geometry of a node relative to a candidate facet
    struct Term
    {
        std::size_t facet;
        std::array<Scalar, dim> closestPointWeights;
        std::array<Scalar, dim> footPointWeights;
        GlobalPosition fromClosestPoint;
        GlobalPosition normal;
        Scalar gap;
        Scalar squaredDistance;
    };

    // the energy of a node and the coefficients of its gradient: sum over the facets of
    // alpha*grad(d_f^2) + beta*grad(g_f), without the factor A_s/2
    struct Blend
    {
        Scalar energy = 0.0;
        Scalar pressure = 0.0;
        std::vector<Scalar> share, alpha, beta;
    };

    struct NodeTerms
    {
        std::size_t node;
        std::vector<Term> terms;
    };

public:
    using JacobianBlock = Dune::FieldMatrix<Scalar, dim, dim>;

    /*!
     * \brief A node in contact range: the facet with the largest share, its closest point, the
     *        mean plane distance and facet normal of the facets and the nodal contact pressure
     *        (the shares and cutoffs times the pressures of the facets)
     */
    struct Pair
    {
        std::size_t node;
        std::size_t facet;
        std::array<Scalar, dim> weights;
        Scalar gap;
        GlobalPosition normal;
        Scalar pressure;
    };

    explicit PenaltyContact(std::shared_ptr<const Surface> surface, const std::string& paramGroup = "")
    : surface_(surface)
    , law_(getParamFromGroup<Scalar>(paramGroup, "Contact.PenaltyStiffness"),
           getParamFromGroup<Scalar>(paramGroup, "Contact.Offset"))
    , searchMargin_(getParamFromGroup<Scalar>(paramGroup, "Contact.SearchMargin", law_.offset()))
    , smoothingLength_(getParamFromGroup<Scalar>(paramGroup, "Contact.SmoothingLength", 0.1*law_.offset()))
    , minAlignment_(getParamFromGroup<Scalar>(paramGroup, "Contact.MinAlignment", 0.5))
    , maxPenetration_(getParamFromGroup<Scalar>(paramGroup, "Contact.MaxPenetration", law_.offset()))
    , candidates_(surface->numNodes())
    {
        const int rings = std::max(1, getParamFromGroup<int>(paramGroup, "Contact.ExclusionRings", 2));
        computeExcludedNodes_(rings);
    }

    const Surface& surface() const
    { return *surface_; }

    const Law& law() const
    { return law_; }

    Scalar searchRadius() const
    { return law_.offset() + searchMargin_; }

    /*!
     * \brief Whether the candidates searched at displacement uOld cover all contacts at u
     *
     * This holds if no surface node moved by more than half the search margin, since then no
     * node came closer to a facet by more than the margin.
     */
    template<class DisplacementVector>
    bool candidatesCover(const DisplacementVector& uOld, const DisplacementVector& u) const
    {
        Scalar maxSquaredIncrement = 0.0;
        for (std::size_t s = 0; s < surface_->numNodes(); ++s)
        {
            const auto dofIdx = surface_->dofIndex(s);
            Scalar squaredIncrement = 0.0;
            for (int i = 0; i < dim; ++i)
                squaredIncrement += (u[dofIdx][i] - uOld[dofIdx][i])*(u[dofIdx][i] - uOld[dofIdx][i]);
            maxSquaredIncrement = std::max(maxSquaredIncrement, squaredIncrement);
        }
        return maxSquaredIncrement <= 0.25*searchMargin_*searchMargin_;
    }

    //! Search the candidate facets of every node in the configuration with displacement u
    template<class DisplacementVector>
    void updateCandidates(const DisplacementVector& u)
    {
        const auto x = surface_->positions(u);
        const auto n = surface_->normals(x);
        candidates_ = findCandidateFacets(*surface_, x, searchRadius(), candidateFilter(x, n));
    }

    /*!
     * \brief Add the facets in contact range at displacement u that are not candidates yet
     * \return whether facets were added, i.e. whether a solution computed with the previous
     *         candidates missed contact
     */
    template<class DisplacementVector>
    bool addMissingCandidates(const DisplacementVector& u)
    {
        const auto x = surface_->positions(u);
        const auto n = surface_->normals(x);
        const auto found = findCandidateFacets(*surface_, x, std::max(law_.offset(), maxPenetration_), candidateFilter(x, n));
        bool added = false;
        for (std::size_t s = 0; s < found.size(); ++s)
        {
            auto& candidates = candidates_[s];
            for (const auto f : found[s])
            {
                const auto pos = std::lower_bound(candidates.begin(), candidates.end(), f);
                if (pos == candidates.end() || *pos != f)
                {
                    candidates.insert(pos, f);
                    added = true;
                }
            }
        }
        return added;
    }

    const std::vector<std::size_t>& candidates(std::size_t nodeIdx) const
    { return candidates_[nodeIdx]; }

    //! The predicate deciding whether a facet within the search radius is a candidate for a node
    auto candidateFilter(const std::vector<GlobalPosition>& x, const Normals& n) const
    {
        return [this, &x, &n](std::size_t s, std::size_t f)
        {
            const auto& excluded = excludedNodes_[s];
            for (const auto node : surface_->facet(f))
                if (std::binary_search(excluded.begin(), excluded.end(), node))
                    return false;

            if (n.node[s]*n.facet[f] >= 0.0)
                return false;

            const auto term = makeTerm_(f, x[s], surface_->corners(f, x));
            return term.gap*term.gap >= minAlignment_*minAlignment_*term.squaredDistance && term.gap >= -maxPenetration_;
        };
    }

    //! The nodes in contact range of a candidate facet
    template<class DisplacementVector>
    std::vector<Pair> pairs(const DisplacementVector& u) const
    {
        const auto x = surface_->positions(u);
        std::vector<Pair> result;
        for (const auto& nodeTerms : nodeTerms_(x))
        {
            const auto blend = blend_(nodeTerms.terms);
            std::size_t dominant = 0;
            Scalar maxWeight = -1.0, weightSum = 0.0, gap = 0.0;
            GlobalPosition normal(0.0);
            for (std::size_t t = 0; t < nodeTerms.terms.size(); ++t)
            {
                const auto& term = nodeTerms.terms[t];
                const Scalar weight = blend.share[t]*cutoff_(term);
                weightSum += weight;
                gap += weight*term.gap;
                normal.axpy(weight, term.normal);
                if (weight > maxWeight)
                {
                    maxWeight = weight;
                    dominant = t;
                }
            }
            if (weightSum <= 0.0)
                continue;
            gap /= weightSum;
            normal /= normal.two_norm();
            const auto& term = nodeTerms.terms[dominant];
            result.push_back(Pair{nodeTerms.node, term.facet, term.closestPointWeights, gap, normal, blend.pressure});
        }
        return result;
    }

    //! The pairs with a positive contact pressure
    template<class DisplacementVector>
    std::vector<Pair> activePairs(const DisplacementVector& u) const
    {
        auto result = pairs(u);
        std::erase_if(result, [](const auto& p){ return p.pressure <= 0.0; });
        return result;
    }

    template<class DisplacementVector>
    Scalar energy(const DisplacementVector& u) const
    {
        const auto x = surface_->positions(u);
        Scalar e = 0.0;
        for (const auto& nodeTerms : nodeTerms_(x))
            e += 0.5*surface_->referenceArea(nodeTerms.node)*blend_(nodeTerms.terms).energy;
        return e;
    }

    //! Add the gradient of the contact energy to the residual, indexed by displacement dofs
    template<class Residual, class DisplacementVector>
    void addResidual(Residual& r, const DisplacementVector& u) const
    {
        const auto x = surface_->positions(u);
        for (const auto& nodeTerms : nodeTerms_(x))
        {
            const auto stencil = makeStencil_(nodeTerms);
            const auto forces = nodeForces_(nodeTerms, stencil, stencilPositions_(stencil, x));
            for (std::size_t i = 0; i < stencil.nodes.size(); ++i)
                r[surface_->dofIndex(stencil.nodes[i])] += forces[i];
        }
    }

    /*!
     * \brief Add the Jacobian of the contact residual
     * \param u the displacement
     * \param addBlock called as addBlock(rowDof, colDof, block) with the derivative of the residual
     *        at rowDof with respect to the displacement at colDof
     */
    template<class DisplacementVector, class AddBlock>
    void addJacobian(const DisplacementVector& u, const AddBlock& addBlock) const
    {
        const auto x = surface_->positions(u);
        const Scalar eps = 1e-4*std::min(law_.offset(), smoothingLength_);
        for (const auto& nodeTerms : nodeTerms_(x))
        {
            const auto stencil = makeStencil_(nodeTerms);
            const auto size = stencil.nodes.size();
            auto pos = stencilPositions_(stencil, x);
            std::vector<JacobianBlock> blocks(size*size, JacobianBlock(0.0));
            for (std::size_t col = 0; col < size; ++col)
            {
                for (int k = 0; k < dim; ++k)
                {
                    const Scalar original = pos[col][k];
                    pos[col][k] = original + eps;
                    const auto plus = nodeForces_(nodeTerms, stencil, pos);
                    pos[col][k] = original - eps;
                    const auto minus = nodeForces_(nodeTerms, stencil, pos);
                    pos[col][k] = original;
                    for (std::size_t row = 0; row < size; ++row)
                        for (int i = 0; i < dim; ++i)
                            blocks[row*size + col][i][k] = (plus[row][i] - minus[row][i])/(2.0*eps);
                }
            }
            for (std::size_t row = 0; row < size; ++row)
                for (std::size_t col = 0; col < size; ++col)
                    addBlock(surface_->dofIndex(stencil.nodes[row]), surface_->dofIndex(stencil.nodes[col]), blocks[row*size + col]);
        }
    }

    //! Call addEntry(rowDof, colDof) for every Jacobian entry the current candidates may need
    template<class AddEntry>
    void jacobianPattern(const AddEntry& addEntry) const
    {
        std::vector<std::size_t> dofs;
        for (std::size_t s = 0; s < surface_->numNodes(); ++s)
        {
            if (candidates_[s].empty())
                continue;
            dofs.assign(1, surface_->dofIndex(s));
            for (const auto f : candidates_[s])
                for (const auto node : surface_->facet(f))
                    dofs.push_back(surface_->dofIndex(node));
            std::sort(dofs.begin(), dofs.end());
            dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
            for (const auto row : dofs)
                for (const auto col : dofs)
                    addEntry(row, col);
        }
    }

private:
    // the surface nodes the energy of a node depends on (the node itself first) and the local
    // corner indices of its facets
    struct Stencil
    {
        std::vector<std::size_t> nodes;
        std::vector<std::array<std::size_t, dim>> corners;
    };

    // the terms of the nodes with candidates within the offset of the law, without the facets
    // of negligible share
    std::vector<NodeTerms> nodeTerms_(const std::vector<GlobalPosition>& x) const
    {
        std::vector<NodeTerms> result;
        for (std::size_t s = 0; s < surface_->numNodes(); ++s)
        {
            if (candidates_[s].empty())
                continue;
            NodeTerms nodeTerms{s, {}};
            bool inRange = false;
            for (const auto f : candidates_[s])
            {
                nodeTerms.terms.push_back(makeTerm_(f, x[s], surface_->corners(f, x)));
                const auto& term = nodeTerms.terms.back();
                inRange = inRange || (term.gap < law_.offset() && cutoff_(term) > 0.0);
            }
            if (!inRange)
                continue;
            const auto blend = blend_(nodeTerms.terms);
            std::vector<Term> kept;
            for (std::size_t t = 0; t < nodeTerms.terms.size(); ++t)
                if (blend.share[t] >= 1e-12)
                    kept.push_back(nodeTerms.terms[t]);
            nodeTerms.terms = std::move(kept);
            result.push_back(std::move(nodeTerms));
        }
        return result;
    }

    static Term makeTerm_(std::size_t f, const GlobalPosition& xs, const Corners& c)
    {
        const auto cp = closestPointOnSimplex(xs, c);
        auto normal = Surface::normal(c);
        normal /= normal.two_norm();
        const Scalar gap = (xs - c[0])*normal;
        auto footPoint = xs;
        footPoint.axpy(-gap, normal);
        return Term{f, cp.weights, barycentricWeights_(footPoint, c), xs - cp.point, normal, gap, cp.squaredDistance};
    }

    // barycentric coordinates of a point in the plane of a segment or triangle
    static std::array<Scalar, dim> barycentricWeights_(const GlobalPosition& p, const Corners& c)
    {
        if constexpr (dim == 2)
        {
            const auto e = c[1] - c[0];
            const Scalar t = (p - c[0])*e/e.two_norm2();
            return {1.0 - t, t};
        }
        else
        {
            const auto e0 = c[1] - c[0];
            const auto e1 = c[2] - c[0];
            const auto v = p - c[0];
            const Scalar d00 = e0*e0, d01 = e0*e1, d11 = e1*e1, d20 = v*e0, d21 = v*e1;
            const Scalar denominator = d00*d11 - d01*d01;
            const Scalar xi1 = (d11*d20 - d01*d21)/denominator;
            const Scalar xi2 = (d00*d21 - d01*d20)/denominator;
            return {1.0 - xi1 - xi2, xi1, xi2};
        }
    }

    Scalar cutoff_(const Term& term) const
    {
        const Scalar q = std::max(term.squaredDistance - term.gap*term.gap, Scalar(0.0));
        const Scalar range2 = searchRadius()*searchRadius();
        return q < range2 ? (1.0 - q/range2)*(1.0 - q/range2) : 0.0;
    }

    Blend blend_(const std::vector<Term>& terms) const
    {
        const auto numTerms = terms.size();
        const Scalar eps = smoothingLength_;
        const Scalar range2 = searchRadius()*searchRadius();
        std::vector<Scalar> distance(numTerms);
        for (std::size_t t = 0; t < numTerms; ++t)
            distance[t] = std::sqrt(terms[t].squaredDistance + 0.01*eps*eps);
        const Scalar minDistance = *std::min_element(distance.begin(), distance.end());

        Blend blend;
        blend.share.resize(numTerms);
        blend.alpha.resize(numTerms);
        blend.beta.resize(numTerms);
        Scalar sum = 0.0;
        for (std::size_t t = 0; t < numTerms; ++t)
        {
            blend.share[t] = std::exp(-(distance[t] - minDistance)/eps);
            sum += blend.share[t];
        }

        std::vector<Scalar> facetEnergy(numTerms), cutoff(numTerms), cutoffDerivative(numTerms), psi(numTerms);
        for (std::size_t t = 0; t < numTerms; ++t)
        {
            blend.share[t] /= sum;
            const Scalar q = std::max(terms[t].squaredDistance - terms[t].gap*terms[t].gap, Scalar(0.0));
            cutoff[t] = q < range2 ? (1.0 - q/range2)*(1.0 - q/range2) : 0.0;
            cutoffDerivative[t] = q < range2 ? -2.0*(1.0 - q/range2)/range2 : 0.0;
            psi[t] = law_.energyDensity(terms[t].gap);
            facetEnergy[t] = cutoff[t]*psi[t];
            blend.energy += blend.share[t]*facetEnergy[t];
            blend.pressure += blend.share[t]*cutoff[t]*law_.pressure(terms[t].gap);
        }

        for (std::size_t t = 0; t < numTerms; ++t)
        {
            const Scalar w = blend.share[t];
            blend.alpha[t] = w*cutoffDerivative[t]*psi[t] - w*(facetEnergy[t] - blend.energy)/(2.0*eps*distance[t]);
            blend.beta[t] = w*(-2.0*cutoffDerivative[t]*psi[t]*terms[t].gap - cutoff[t]*law_.pressure(terms[t].gap));
        }
        return blend;
    }

    Stencil makeStencil_(const NodeTerms& nodeTerms) const
    {
        Stencil stencil{{nodeTerms.node}, {}};
        for (const auto& term : nodeTerms.terms)
        {
            std::array<std::size_t, dim> corners;
            for (int j = 0; j < dim; ++j)
            {
                const auto node = surface_->facet(term.facet)[j];
                const auto it = std::find(stencil.nodes.begin(), stencil.nodes.end(), node);
                corners[j] = std::distance(stencil.nodes.begin(), it);
                if (it == stencil.nodes.end())
                    stencil.nodes.push_back(node);
            }
            stencil.corners.push_back(corners);
        }
        return stencil;
    }

    std::vector<GlobalPosition> stencilPositions_(const Stencil& stencil, const std::vector<GlobalPosition>& x) const
    {
        std::vector<GlobalPosition> pos(stencil.nodes.size());
        for (std::size_t i = 0; i < pos.size(); ++i)
            pos[i] = x[stencil.nodes[i]];
        return pos;
    }

    // residual contributions on the stencil nodes at their positions pos
    std::vector<GlobalPosition> nodeForces_(const NodeTerms& nodeTerms, const Stencil& stencil, const std::vector<GlobalPosition>& pos) const
    {
        std::vector<Term> terms;
        for (std::size_t t = 0; t < nodeTerms.terms.size(); ++t)
        {
            Corners corners;
            for (int j = 0; j < dim; ++j)
                corners[j] = pos[stencil.corners[t][j]];
            terms.push_back(makeTerm_(nodeTerms.terms[t].facet, pos[0], corners));
        }
        const auto blend = blend_(terms);
        const Scalar factor = 0.5*surface_->referenceArea(nodeTerms.node);

        std::vector<GlobalPosition> forces(pos.size(), GlobalPosition(0.0));
        for (std::size_t t = 0; t < terms.size(); ++t)
        {
            const Scalar alpha = 2.0*factor*blend.alpha[t];
            const Scalar beta = factor*blend.beta[t];
            forces[0].axpy(alpha, terms[t].fromClosestPoint);
            forces[0].axpy(beta, terms[t].normal);
            for (int j = 0; j < dim; ++j)
            {
                auto& f = forces[stencil.corners[t][j]];
                f.axpy(-alpha*terms[t].closestPointWeights[j], terms[t].fromClosestPoint);
                f.axpy(-beta*terms[t].footPointWeights[j], terms[t].normal);
            }
        }
        return forces;
    }

    void computeExcludedNodes_(int rings)
    {
        const auto numNodes = surface_->numNodes();
        std::vector<std::vector<std::size_t>> neighbors(numNodes);
        for (std::size_t s = 0; s < numNodes; ++s)
        {
            for (const auto f : surface_->nodeFacets(s))
                for (const auto node : surface_->facet(f))
                    if (node != s)
                        neighbors[s].push_back(node);
            std::sort(neighbors[s].begin(), neighbors[s].end());
            neighbors[s].erase(std::unique(neighbors[s].begin(), neighbors[s].end()), neighbors[s].end());
        }

        excludedNodes_.assign(numNodes, {});
        for (std::size_t s = 0; s < numNodes; ++s)
        {
            auto& excluded = excludedNodes_[s];
            excluded.push_back(s);
            std::vector<std::size_t> front{s};
            for (int ring = 1; ring < rings; ++ring)
            {
                std::vector<std::size_t> next;
                for (const auto node : front)
                    for (const auto neighbor : neighbors[node])
                        if (std::find(excluded.begin(), excluded.end(), neighbor) == excluded.end())
                        {
                            excluded.push_back(neighbor);
                            next.push_back(neighbor);
                        }
                front = std::move(next);
            }
            std::sort(excluded.begin(), excluded.end());
        }
    }

    std::shared_ptr<const Surface> surface_;
    Law law_;
    Scalar searchMargin_;
    Scalar smoothingLength_;
    Scalar minAlignment_;
    Scalar maxPenetration_;
    std::vector<std::vector<std::size_t>> candidates_;
    std::vector<std::vector<std::size_t>> excludedNodes_;
};

} // end namespace Dumux::Contact

#endif
