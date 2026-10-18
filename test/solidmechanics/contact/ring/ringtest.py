#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Runs the ring collapse test and checks it against the inextensible elastica.

The reference is the elastica ring of bending stiffness D under uniform external pressure p,
symmetric about both axes, whose midline touches the opposite one at the half wall thickness
h = t/2 (Flaherty, Keller, Rubinow, SIAM J. Appl. Math. 23, 1972, for h = 0). Its values for
h/R = 0.02, from shooting on D theta' = C + p r^2/2 - Q x/2 over the quarter ring: the contact
force is nearly linear in p for 6 <= pR^3/D <= 9, with slope dQ/dp = 1.566 R and zero at
p = 5.137 D/R^3; line contact (zero curvature at the contact point) starts at p = 10.07 D/R^3.
Checked: the midline stays outside h on the minor axis, the slope within 1 %, the zero of the
linear fit within 3 %, and the largest contact pressure leaves the minor axis between the two
listed pressures around 10.07.

usage: ringtest.py <executable> <parameter file> [arguments]
"""
import subprocess
import sys

executable, params = sys.argv[1:3]
output = subprocess.run(["./" + executable, params] + sys.argv[3:], check=True, capture_output=True, text=True).stdout
print(output)

rows = [[float(v) for v in line.split()] for line in output.splitlines() if line[:1].isdigit() and len(line.split()) == 7]
if not rows:
    sys.exit("no results")
p, semiMinorAxis, semiMajorAxis, area, force, peakPosition, numPairs = map(list, zip(*rows))

h, firstContact, slope, lineContact = 0.02, 5.137, 1.566, 10.07
errors = []
if min(semiMinorAxis) < h:
    errors.append(f"walls crossed: semi-minor axis {min(semiMinorAxis)} below {h}")
fitted = [(pk, qk) for pk, qk in zip(p, force) if 6.0 <= pk <= 9.0]
if len(fitted) < 2 or min(q for _, q in fitted) <= 0.0:
    errors.append("no contact force between p = 6 and 9")
else:
    n = len(fitted)
    meanP, meanQ = sum(x for x, _ in fitted)/n, sum(y for _, y in fitted)/n
    fittedSlope = sum((x - meanP)*(y - meanQ) for x, y in fitted)/sum((x - meanP)**2 for x, _ in fitted)
    fittedFirstContact = meanP - meanQ/fittedSlope
    print(f"slope {fittedSlope:.4f} (elastica {slope}), first contact {fittedFirstContact:.4f} (elastica {firstContact})")
    if abs(fittedSlope/slope - 1.0) > 0.01:
        errors.append(f"slope dQ/dp {fittedSlope} differs from {slope}")
    if abs(fittedFirstContact/firstContact - 1.0) > 0.03:
        errors.append(f"first contact at p = {fittedFirstContact} differs from {firstContact}")
below = max(k for k in range(len(p)) if p[k] < lineContact)
above = min(k for k in range(len(p)) if p[k] > lineContact)
if peakPosition[below] > 1e-8 or peakPosition[above] < 1e-8:
    errors.append(f"line contact does not start between p = {p[below]} and {p[above]}")
if errors:
    sys.exit("\n".join(errors))
print("ring collapse agrees with the elastica")
