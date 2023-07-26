This folder contains tools for creating pore network models from images.
SNOW and dual-SNOW from the `porespy` Python package are used and `OpenPNM` is used for some conversions.
There is a tool to generate a DUNE Grid Format file.
To try the script, there is a sample `raw` binary image in the folder. For testing run

* `python3 extract_pore_network_with_porespy.py sample.raw -r 1e-5 -v 1000 1000`
* `python3 openpnm2dgf.py sample.pnm`

which generates `sample.dgf`.

Python versions tested for the above mentioned python scripts:
- Python 3.10.6

TODO: These steps should be tested automatically in the future in a CI setup with `porespy` and `OpenPNM` installed.
