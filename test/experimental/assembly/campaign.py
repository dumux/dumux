# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Run the benchmark repeatedly and report the minimum of each column.

A single run varies by as much as 60% on a loaded machine, which is larger than most effects worth
measuring. The minimum over repeated *processes* is the statistic that survives that: it answers
"how fast can this go", and interference can only ever make a run slower.

    python campaign.py [runs]        # from the build directory
"""

import collections
import subprocess
import sys

runs = int(sys.argv[1]) if len(sys.argv) > 1 else 9
best = collections.defaultdict(lambda: [float("inf")]*3)
meta = {}

for _ in range(runs):
    out = subprocess.run(["./bench_experimental_cvfe_analytic_jacobian"],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if not line.startswith("DATA;"):
            continue
        _, name, dofs, local, analytic, numeric, residual = line.split(";")
        meta[name] = (int(dofs), int(local))
        row = best[name]
        for i, value in enumerate((float(analytic), float(numeric), float(residual))):
            row[i] = min(row[i], value)

print(f"{'scheme':<18}{'dofs':>7}{'local':>7}{'analytic':>12}{'numeric':>11}"
      f"{'residual':>11}{'speedup':>10}")
for name, (dofs, local) in meta.items():
    analytic, numeric, residual = best[name]
    print(f"{name:<18}{dofs:>7}{local:>7}{analytic:>9.1f} ms{numeric:>8.1f} ms"
          f"{residual:>8.1f} ms{numeric/analytic:>9.2f}x")
