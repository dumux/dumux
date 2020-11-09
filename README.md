<img src="doc/logo/dumux_logo_hires_whitebg.png" alt="dumux logo" width="400"/>

What is DuMu<sup>x</sup>?
===============

[DuMu<sup>x</sup>][0] is a simulation toolbox mainly aimed at flow and transport
processes in porous media. DuMu<sup>x</sup> is based on the [DUNE][1]
framework and aims to provide a multitude of numerical models as well
as flexible discretization methods for complex non-linear phenomena,
such as CO2 sequestration, soil remediation, drug delivery in cancer
therapy and more. Have a look at our publications
(see below: [How to cite](#how-to-cite))
for a more detailed description of the goals and motivations behind DuMu<sup>x</sup>.


Installation
===============

Have a look at the [installation guide][3] or use the [DuMu<sup>x</sup> handbook][4],
Chapter 2.

Documentation
==============

The following resources are useful to get started with DuMu<sup>x</sup>:

* [Getting started guide](https://dumux.org/gettingstarted/) on the [DuMu<sup>x</sup> website](https://dumux.org/)
* [Handbook](https://dumux.org/handbook), a detailed DuMu<sup>x</sup> manual,
* [DuMu<sup>x</sup> course materials](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/tree/master),
* [Examples](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/tree/master/examples), with detailed description of code and results,
* [Class documentation](https://dumux.org/docs/doxygen/master/) generated from the source code,
* [Mailing list](https://listserv.uni-stuttgart.de/mailman/listinfo/dumux),
* [Changelog](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CHANGELOG.md), where all changes between different release versions are listed and explained.

Some helpful code snippets are available in the [Wiki](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/wikis/home).

License
========

DuMu<sup>x</sup> is licensed under the terms and conditions of the GNU General
Public License (GPL) version 3 or - at your option - any later
version. The GPL can be [read online][5] or in the [LICENSE.md](LICENSE.md) file
provided in the topmost directory of the DuMu<sup>x</sup> source code tree.

Please note that DuMu<sup>x</sup>' license, unlike DUNE's, does *not* feature a
template exception to the GNU General Public License. This means that
you must publish any source code which uses any of the DuMu<sup>x</sup> header
files if you want to redistribute your program to third parties. If
this is unacceptable to you, please [contact us][6] for a commercial
license.

See the file [LICENSE.md](LICENSE.md) for full copying permissions.

How to cite
============

DuMu<sup>x</sup> is research software and developed at research institutions.
If you are using DuMu<sup>x</sup> in scientific publications and in
the academic context, please cite (at least one)
of our publications:

* [Koch, T., Gläser, D., Weishaupt, K., Ackermann, S., Beck, M., Becker, B.,
  Burbulla, S., Class, H., Coltman, E., Emmert, S., Fetzer, T., Grüninger, C.,
  Heck, K., Hommel, J., Kurz, T., Lipp, M., Mohammadi, F., Scherrer, S.,
  Schneider, M., Seitz, G., Stadler, L., Utz, M., Weinhardt, F.
  & Flemisch, B. (_2020_). __DuMu<sup>x</sup> 3 – an open-source simulator for solving flow
  and transport problems in porous media with a focus on model coupling.__
  _Computers & Mathematics with Applications_,
  https://doi.org/10.1016/j.camwa.2020.02.012][7]

* [Flemisch, B., Darcis, M., Erbertseder, K., Faigle, B., Lauser, A.,
  Mosthaf, K., Müthing, S., Nuske, P., Tatomir, A., Wolff, M.,
  & Helmig, R. (_2011_). __DuMu<sup>x</sup>: DUNE for multi-{phase,component,scale,physics,…}
  flow and transport in porous media__.
  _Advances in Water Resources_, 34(9), 1102–1112.
  https://doi.org/10.1016/j.advwatres.2011.03.007][2]

You can also cite specific releases published on Zenodo:
[![zenodo badge](https://zenodo.org/badge/DOI/10.5281/zenodo.2479594.svg)](https://doi.org/10.5281/zenodo.2479594)



Automated Testing
==================

* DuMu<sup>x</sup> 3.3 / Dune 2.7: [![buildbot badge](https://git.iws.uni-stuttgart.de/buildbot/badges/dumux-rel-3-3-dune-rel-2-7-clang-full.svg)](https://git.iws.uni-stuttgart.de/buildbot/#/builders)
* DuMu<sup>x</sup> master / Dune 2.7: [![buildbot badge](https://git.iws.uni-stuttgart.de/buildbot/badges/dumux-master-dune-rel-2-7-clang-full.svg)](https://git.iws.uni-stuttgart.de/buildbot/#/builders)
* DuMu<sup>x</sup> master / Dune master: [![buildbot badge](https://git.iws.uni-stuttgart.de/buildbot/badges/dumux-master-dune-master-clang-full.svg)](https://git.iws.uni-stuttgart.de/buildbot/#/builders)

DuMu<sup>x</sup> features many tests (some unit tests and test problems) that can
be run manually. We have experimental support for automated testing with buildbot.
Click <a href="https://git.iws.uni-stuttgart.de/buildbot/#/builders" target="_blank">here (buildbot)</a>
to see the latest builds (clicking on a build
number will show a detailed overview of the build).

[![coverage report](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-coverage/badges/master/coverage.svg)](https://pages.iws.uni-stuttgart.de/dumux-repositories/dumux-coverage/)

A weekly coverage report of the test suite is created by gcovr/gcov. The report
currently doesn't include non-instantiated code, so the real coverage is likely lower. However,
only a few lines of code are never instatiated in the comprehensive test suite.


Contributing
=============

Contributions are highly welcome. Please ask questions over the [mailing list](mailto:dumux@listserv.uni-stuttgart.de).
Please review the [contribution guidelines](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CONTRIBUTING.md)
before opening issues and merge requests. For bug reports contact us
over the mailing list, or file an [issue](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/issues). For bug fixes,
feature implementations open a [merge request](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/merge_requests)
or send us formatted patches.

Backwards Compatibility
=======================

For a detailed description of the backwards compatibility policy, 
please see [contribution guidelines](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CONTRIBUTING.md).

DuMu<sup>x</sup> releases are split into major(e.g. 2.0, 3.0) and minor (e.g. 3.1, 3.2, 3.3) releases. 
Major releases are not required to maintain backwards compatibility (see below), 
but would provide a detailed guide on how to update dependent modules. 
For each minor release, maintaining backwards compatibility is strongly encouraged and recommended.

Despite the goal of maintaining backwards compatibility across minor releases,
for more complicated changes, this is decided upon on a case to case basis, due to limited developer resources. 
In the case that implementing full backwards compatibility for an update is not feasible, or would require unreasonable resources, 
the degree of backwards compatibility be decided by a vote in one of the monthly core developer meetings.

Major version update, 2.12 to 3.0
===================================

With the version update to version 3, many features have been added and a lot has been improved in DuMu<sup>x</sup>. See the
[changelog](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CHANGELOG.md) for a list of changes.
If you decide to update from version 2.12, please have a look at our small
[guide](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/wikis/Hints/Upgrade-to-Dumux-3)
on how to update an application to the new version.

[0]: https://dumux.org
[1]: https://dune-project.org/
[2]: https://dumux.org/documents/dumux_awrpaper.pdf
[3]: https://dumux.org/installation
[4]: https://dumux.org/documents/dumux-handbook-master.pdf
[5]: https://www.gnu.org/licenses/gpl-3.0.en.html
[6]: https://www.iws.uni-stuttgart.de/en/lh2/
[7]: https://doi.org/10.1016/j.camwa.2020.02.012
