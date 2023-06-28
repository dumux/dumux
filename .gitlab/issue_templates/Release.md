<!--
SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
SPDX-License-Identifier: CC0-1.0
-->

<!--
This form is for release issue ONLY!
If you're looking for help check out the [readme](/README.md).
-->
This a release issue template with a checklist based on [Release Manager Tasks](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/wikis/Release-Manager-Tasks). Check there for details on each task.
# [5 Weeks] prior to the release:
- [ ]  Call a meeting
- [ ]  Assign developers to the major subtasks
    - Manager:
    - Lecture:
    - Doxygen:
    - Course:
    - Examples:
    - Website:
- [ ]  Create a group (dumux-repositories) milestone in GitLab
- [ ]  Fix planned [Dune] and compiler compatibility
- [ ]  Go through all existing Gitlab tasks ([issues] or [MRs])
- [ ]  Assign the open [issues] and [MRs] to the developers
- [ ]  Post a release schedule announcement on the mailing list

# [3 weeks] prior to the release:
- [ ]  Check all open MRs and Issues for severity and impact
- [ ]  Update the milestone issue
- [ ]  Announce soft feature freeze on the mailing list
- [ ]  Check in with the managers of the sub-tasks

# [2 weeks] prior to the release:
__Hard Feature Freeze!__
- [ ]  Check remaining MRs and Issues
- [ ]  Update the milestone issue
- [ ]  Generate example documentation with the script generate_example_docs.py
- [ ]  Update CHANGELOG
- [ ]  Announce hard feature freeze on the mailing list
- [ ]  Header check (run `make headercheck`)
- [ ]  Re-generate parameter list (`bin/doc/generate_parameterlist.py`) and check log file (maybe manual updates to `bin/doc/docextra/parameters.json` are needed)
- [ ]  Make sure the CMakeLists.txt in `dumux` subfolder are up-to-date (generated with `bin/create_cmakelists.py`, dumux `make install` should result in a usable installed dumux version)
- [ ]  Create release branches
- [ ]  Configure CI to test release branch
- [ ]  Local Testing (different compilers and dependency setups)

# During the week of the release:
__Final testing!__
- [ ]  Re-run tests
- [ ]  Update copyright
- [ ]  Update License
- [ ]  Test Lecture and Course
- [ ]  Update Website
- [ ]  Create a release candidate
- [ ]  Call for testing

# Releasing DuMu<sup>x</sup>
__Important:__ These steps are normally done together with Bernd. Make sure to schedule the time properly.
- [ ]  Protect release branch (in GitLab)
- [ ]  Create new tags
- [ ]  Prepare a Zenodo citation
- [ ]  Include the Zenodo citation to the website (only major releases)

# After the release
- [ ]  Bump version in dune.module to next version
- [ ]  Write a release email to dumux@listserv.uni-stuttgart.de about updates on supported features and upcoming changes
