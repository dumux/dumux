<img src="doc/logo/dumux_logo_hires_whitebg.png" alt="dumux logo" width="400"/>

What is DuMuX?
===============

[DuMuX][0] is a simulation toolbox mainly aimed at flow and transport
processes in porous media. DuMuX is based on the [DUNE][1]
framework and aims to provide a multitude of numerical models as well
as flexible discretization methods for complex non-linear phenomena,
such as CO2 sequestration, soil remediation, drug delivery in cancer
therapy and more. See [our publication][2] for a more detailed
description of the goals and motivations behind DuMuX.


Installation
===============

Have a look at the [installation guide][3] or use the [DuMuX handbook]
[4], chapter 2.


License
========

DuMuX is licensed under the terms and conditions of the GNU General
Public License (GPL) version 2 or - at your option - any later
version. The GPL can be [read online][5] or in the [LICENSE.md](LICENSE.md) file
provided in the topmost directory of the DuMuX source code tree.

Please note that DuMuX' license, unlike DUNE's, does *not* feature a
template exception to the GNU General Public License. This means that
you must publish any source code which uses any of the DuMuX header
files if you want to redistribute your program to third parties. If
this is unacceptable to you, please [contact us][6] for a commercial
license.

See the file [LICENSE.md](LICENSE.md) for full copying permissions.

Automated Testing
==================
[![buildbot badge](https://git.iws.uni-stuttgart.de/buildbot/badges/dumux-master-dune-release26-gcc.svg)](https://git.iws.uni-stuttgart.de/buildbot/#/builders)

DuMuX features many tests (some unit tests and test problems) that can
be run manually. We have experimental support for automated testing with buildbot.
Click <a href="https://git.iws.uni-stuttgart.de/buildbot/#/builders" target="_blank">here (buildbot)</a>
to see the latest builds (clicking on a build
number will show a detailed overview of the build).

Contributing
=============

Contributions are highly welcome. Please ask questions over the [mailing list](mailto:dumux@listserv.uni-stuttgart.de).
Please review the [contribution guidelines](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CONTRIBUTING.md)
before opening issues and merge requests. For bug reports contact us
over the mailing list, or file an [issue](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/issues). For bug fixes,
feature implementations open a [merge request](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/merge_requests)
or send us formatted patches.

Major version update, 2.12 to 3.0
===================================

With the version update to version 3, many features have been added and a lot has been improved in DuMuX. See the
[changelog](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/CHANGELOG.md) for a list of changes.
If you decide to update from version 2.12, please have a look at our small
[guide](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/wikis/Updating-programs-from-version-2.12-to-version-3.0)
on how to update an application to the new version.

[0]: http://dumux.org
[1]: http://dune-project.org
[2]: http://dumux.org/documents/dumux_awrpaper.pdf
[3]: http://www.dumux.org/installation.php
[4]: http://www.dumux.org/documents/dumux-handbook-2.12.pdf
[5]: http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
[6]: http://www.hydrosys.uni-stuttgart.de/index.en.php
