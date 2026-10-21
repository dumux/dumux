#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Generate solution.hh: exact solution and source terms of a plane-strain manufactured
solution for the large-deformation poroelastic model on the unit square.

The displacement d = t phi(x,y) (a, b) has a rank-one gradient, so J = 1 + t (a phi_x + b phi_y)
is linear in time, and so is the fluid pressure p_f = t psi(x,y). The storage term
alpha_B J + S_p p_f is therefore linear in time and implicit Euler integrates it without time
discretization error. Displacement and fluid pressure vanish on the boundary for all t.

The source terms are expressed in the derivatives of phi and psi, in J and in I1 = tr(F^T F).
"""

import numpy as np
import sympy as sp

x, y, t = sp.symbols("x y t", real=True)
G, K, alpha, Sp, kappa0, mu = sp.symbols("G K alpha Sp kappa0 mu", positive=True)
a, b, J, I1 = sp.symbols("a b J I1", real=True)

phiExact = sp.sin(sp.pi*x)*sp.sin(sp.pi*y)*sp.exp(x/2 + y/4)
psiExact = sp.sin(sp.pi*x)*sp.sin(sp.pi*y)*sp.exp(-y/2)

# amplitude such that the largest |1 - J| is 0.45 at t = 1; J spans [0.55, 1.29]
X, Y = np.meshgrid(np.linspace(0, 1, 401), np.linspace(0, 1, 401))
slope = (2*sp.diff(phiExact, x) + sp.diff(phiExact, y))/sp.sqrt(5)
slope = sp.lambdify((x, y), slope, "numpy")(X, Y)
amplitude = sp.Rational(45, 100)/sp.nsimplify(round(float(np.abs(slope).max()), 4))
aValue, bValue = float(2*amplitude/sp.sqrt(5)), float(amplitude/sp.sqrt(5))

phi, psi = sp.Function("phi")(x, y), sp.Function("psi")(x, y)
Jf, I1f = sp.Function("J")(x, y), sp.Function("I1")(x, y)
F = sp.eye(2) + t*sp.Matrix([a, b])*sp.Matrix([[sp.diff(phi, x), sp.diff(phi, y)]])
cofF = F.adjugate().T
JExact = F.det()
I1Exact = (F.T*F).trace()
assert sp.diff(JExact, t, 2) == 0

pf = t*psi
ps = K*(Jf**2 - 1)/(2*Jf)
P = G*Jf**sp.Rational(-2, 3)*F - G*(I1f + 1)/3*Jf**sp.Rational(-5, 3)*cofF + (ps - alpha*pf)*cofF
divP = sp.Matrix([sp.diff(P[i, 0], x) + sp.diff(P[i, 1], y) for i in range(2)]).subs({
    sp.Derivative(Jf, v): sp.diff(JExact, v) for v in (x, y)} | {
    sp.Derivative(I1f, v): sp.diff(I1Exact, v) for v in (x, y)})

# J^2 C^{-1} = cof(F)^T cof(F)
q = kappa0/mu*cofF.T*cofF*sp.Matrix([sp.diff(pf, x), sp.diff(pf, y)])
fluidSource = sp.diff(alpha*JExact + Sp*pf, t) - (sp.diff(q[0], x) + sp.diff(q[1], y))

derivatives = [(), (x,), (y,), (x, x), (x, y), (y, y)]
suffixes = ["", "X", "Y", "XX", "XY", "YY"]
symbols = {name: [sp.Symbol(name + s) for s in suffixes] for name in ("phi", "psi")}
placeholders = {Jf: J, I1f: I1}
for f, name in [(phi, "phi"), (psi, "psi")]:
    for v, s in zip(derivatives, symbols[name]):
        placeholders[sp.Derivative(f, *v) if v else f] = s


def toSymbols(e):
    return sp.simplify(e.doit().xreplace(placeholders))


definitions = {J: toSymbols(JExact), I1: toSymbols(I1Exact)}

fJ = JExact.doit().subs({phi: phiExact, a: aValue, b: bValue}).doit()
fJ = sp.lambdify((x, y, t), fJ, "numpy")
print(f"J at t=1: min {fJ(X, Y, 1.0).min():.4f}, max {fJ(X, Y, 1.0).max():.4f}")


def code(e):
    return sp.cxxcode(e, standard="c++17")


def wrap(text, indent):
    lines = [""]
    for term in text.replace(" + ", "\0+ ").replace(" - ", "\0- ").split("\0"):
        if lines[-1] and indent + len(lines[-1]) + len(term) > 100:
            lines.append("")
        lines[-1] += (" " if lines[-1] else "") + term
    return ("\n" + " "*indent).join(lines)


def cxxFunction(name, args, exprs):
    used = set().union(*(e.free_symbols for e in exprs))
    for s in (I1, J):
        if s in used:
            used |= definitions[s].free_symbols
    body = []
    for f in ("phi", "psi"):
        if used & set(symbols[f]):
            body.append(f"const auto [{', '.join(map(str, symbols[f]))}] = {f}Derivatives(x, y);")
    body += [f"const Scalar {s} = {wrap(code(definitions[s]), 8)};" for s in (J, I1) if s in used]
    if len(exprs) == 1:
        ret = "Scalar"
        body.append(f"return {wrap(code(exprs[0]), 8)};")
    else:
        ret = f"std::array<Scalar, {len(exprs)}>"
        oneLine = "return {" + ", ".join(code(e) for e in exprs) + "};"
        multiLine = ",\n        ".join(wrap(code(e), 8) for e in exprs)
        body.append(oneLine if len(oneLine) < 96 else f"return {{\n        {multiLine}\n    }};")
    signature = ", ".join(f"Scalar {arg}" for arg in args)
    body = "".join(f"    {line}\n" for line in body)
    return f"template<class Scalar>\n{ret} {name}({signature})\n{{\n{body}}}\n"


def derivativesFunction(name, exact):
    exprs = [sp.simplify(sp.diff(exact, *v)) if v else exact for v in derivatives]
    names = ", ".join(f"{name}_{s.lower()}" for s in suffixes[1:])
    values = ",\n        ".join(wrap(code(e), 8) for e in exprs)
    return (f"//! {name} and its derivatives {names}\n"
            f"template<class Scalar>\nstd::array<Scalar, 6> {name}Derivatives(Scalar x, Scalar y)\n"
            f"{{\n    return {{\n        {values}\n    }};\n}}\n")


functions = [
    "//! amplitude of the displacement d = t phi (a, b)\n"
    f"inline constexpr double a = {aValue!r};\ninline constexpr double b = {bValue!r};\n",
    derivativesFunction("phi", phiExact),
    derivativesFunction("psi", psiExact),
    cxxFunction("displacement", ["x", "y", "t"], [a*t*symbols["phi"][0], b*t*symbols["phi"][0]]),
    cxxFunction("fluidPressure", ["x", "y", "t"], [t*symbols["psi"][0]]),
    cxxFunction("solidPressure", ["x", "y", "t", "K"], [toSymbols(ps)]),
    cxxFunction("momentumSource", ["x", "y", "t", "G", "K", "alpha"],
                [toSymbols(-divP[0]), toSymbols(-divP[1])]),
    cxxFunction("fluidSource", ["x", "y", "t", "alpha", "Sp", "kappa0", "mu"],
                [toSymbols(fluidSource)]),
]

header = """// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \\file
 * \\ingroup PoromechanicsTests
 * \\brief Exact solution and source terms of the plane-strain manufactured solution
 *        for the large-deformation poroelastic model (generated by generate_solution.py).
 */
#ifndef DUMUX_TEST_POROELASTIC_LARGE_DEF_MMS_SOLUTION_HH
#define DUMUX_TEST_POROELASTIC_LARGE_DEF_MMS_SOLUTION_HH

#include <array>
#include <cmath>

namespace Dumux::PoroElasticLargeDefMMS {

"""
footer = """
} // end namespace Dumux::PoroElasticLargeDefMMS

#endif
"""
with open("solution.hh", "w") as out:
    out.write(header + "\n".join(functions) + footer)
print("wrote solution.hh")
