#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


""" Generate tables for CO2 fluid properties.

The tables are generated using the NIST (National Institute of Standards
and Technology) Standard Reference Database Number 69
(https://doi.org/10.18434/T4D303).

Copyright for NIST Standard Reference Data is governed by the Standard
Reference Data Act (https://www.nist.gov/srd/public-law).

######################################################################
In case you are using this the data generated with this script
please cite the following publications:

P.J. Linstrom and W.G. Mallard, Eds.,
NIST Chemistry WebBook, NIST Standard Reference Database Number 69,
National Institute of Standards and Technology, Gaithersburg MD, 20899,
https://doi.org/10.18434/T4D303, (retrieved [insert date]).

Span, Roland, and Wolfgang Wagner.
"A new equation of state for carbon dioxide covering
the fluid region from the triple‐point temperature
to 1100 K at pressures up to 800 MPa."
Journal of physical and chemical reference data 25.6 (1996): 1509-1596.
https://doi.org/10.1063/1.555991

######################################################################

The density and the enthalpy are calculated using the equation of Span and
Wagner (2009 "A New Equation of State for Carbon Dioxide Covering the Fluid
Region from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa").
Therefore, the maximum pressure limit is the lowest of the following values:
* 800.0000 MPa
* The pressure at which a density of 1178.5 kg/m3 is reached.

Previously, the CO2 tables in Dumux had been generated with the original Fortran
code of Span and Wagner. However, this code is not published
under an open-source license. This script is instead based on a free web service
of NIST providing access to the Span and Wagner EOS data.

The enthalpy values provided by NIST are given with respect to the IIR convention
reference state: the enthalpy is set to 200000 J/kg at 0°C for the saturated liquid.
The original Span and Wagner code uses a different reference state. Therefore, to
obtain values corresponding to the Span and Wagner code, we add a constant offset
to the enthalpy, computed as the difference between the Fortran code and the NIST
web service values at the IIR convention reference state.
Additionally, the enthalpy tables provided with Dumux contain a constant offset
(of now unknown origin) with respect to the Span and Wagner code (2.190963e+04 J/kg).
For compatibility reasons, we continue to add this offset in this script.
The total enthalpy offset amounts to -484870.2958311295 J/kg w.r.t. NIST.
"""

from datetime import date
from io import StringIO
from string import Template
import argparse
import urllib
import requests  # pylint: disable=import-error
import numpy as np  # pylint: disable=import-error

parser = argparse.ArgumentParser(
    description="This script generates tables for CO2 fluid properties \n"
    "(density and enthalpy) using the equation of Span and Wagner.\n"
)
parser.add_argument(
    "-t1", "--min_temp", required=True, type=float, help="The minimum temperature in Kelvin."
)
parser.add_argument(
    "-t2", "--max_temp", required=True, type=float, help="The maximum temperature in Kelvin."
)
parser.add_argument(
    "-nt",
    "--n_temp",
    required=True,
    type=int,
    help="The number of temperature sampling points."
    "min_temp is the first sampling point, max_temp the last.",
)
parser.add_argument(
    "-p1", "--min_press", required=True, type=float, help="The minimum pressure in Pascal."
)
parser.add_argument(
    "-p2", "--max_press", required=True, type=float, help="The maximum pressure in Pascal."
)
parser.add_argument(
    "-np",
    "--n_press",
    required=True,
    type=int,
    help="The number of pressure sampling points."
    "min_press is the first sampling point, max_press the last.",
)
cmdArgs = vars(parser.parse_args())

delta_temperature = (cmdArgs["max_temp"] - cmdArgs["min_temp"]) / (cmdArgs["n_temp"] - 1)
delta_pressure = (cmdArgs["max_press"] - cmdArgs["min_press"]) / (cmdArgs["n_press"] - 1)

DENSITY_DATA = []
ENTHALPY_DATA = []

# get the data
for i in range(cmdArgs["n_temp"]):
    temperature = cmdArgs["min_temp"] + i * delta_temperature
    query = {
        "Action": "Data",
        "Wide": "on",
        "ID": "C124389",
        "Type": "IsoTherm",
        "Digits": "12",
        "PLow": str(cmdArgs["min_press"]),
        "PHigh": str(cmdArgs["max_press"]),
        "PInc": str(delta_pressure),
        "T": str(temperature),
        "RefState": "DEF",
        "TUnit": "K",
        "PUnit": "Pa",
        "DUnit": "kg/m3",
        "HUnit": "kJ/kg",
        "WUnit": "m/s",
        "VisUnit": "uPas",
        "STUnit": "N/m",
    }
    response = requests.get(
        "https://webbook.nist.gov/cgi/fluid.cgi?" + urllib.parse.urlencode(query)
    )
    response.encoding = "utf-8"
    text = response.text
    phase = np.genfromtxt(StringIO(text), delimiter="\t", dtype=str, usecols=[-1], skip_header=1)
    values = np.genfromtxt(StringIO(text), delimiter="\t", names=True)

    # NIST provides additional samples at the transition points (if there is a
    # phase transition within the requested data range). Since the code which
    # uses the tables generated by this script can't deal with these additional
    # sample points, they are removed.
    phaseBoundaryIndices = []
    for j in range(1, len(phase) - 1):
        if phase[j] != phase[j + 1]:
            phaseBoundaryIndices += [j, j + 1]

    density = np.delete(values["Density_kgm3"], phaseBoundaryIndices)
    enthalpy = np.delete(values["Enthalpy_kJkg"], phaseBoundaryIndices)
    # transform unit (kJ/kg -> J/kg)
    enthalpy *= 1000
    # transform the reference state
    enthalpy -= 484870.2958311295
    # format the data
    DENSITY_DATA.append("    {" + ", ".join([format(x, ".12e") for x in density]) + "}")
    ENTHALPY_DATA.append("    {" + ", ".join([format(x, ".12e") for x in enthalpy]) + "}")

DENSITY_DATA = ",\n".join(DENSITY_DATA)
ENTHALPY_DATA = ",\n".join(ENTHALPY_DATA)

# write the table by filling the gaps in the template
with open("co2table.hh.template", "r") as templateFile:
    template = Template(templateFile.read())

replacements = {
    "MIN_TEMP": format(cmdArgs["min_temp"]),
    "MAX_TEMP": format(cmdArgs["max_temp"]),
    "NUM_TEMP_SAMPLES": format(cmdArgs["n_temp"]),
    "MIN_PRESS": format(cmdArgs["min_press"]),
    "MAX_PRESS": format(cmdArgs["max_press"]),
    "NUM_PRESS_SAMPLES": format(cmdArgs["n_press"]),
    "DENSITY_VALS": DENSITY_DATA,
    "ENTHALPY_VALS": ENTHALPY_DATA,
    "DATE": date.today().strftime("%B %d, %Y"),
}

with open("co2table.hh", "w") as tables:
    tables.write(template.substitute(replacements))
