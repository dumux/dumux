<!--
This form is for release issue ONLY!
If you're looking for help check out the [readme](/README.md).
-->
This a release issue template with a checklist based on [Release Manager Tasks](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/wikis/Release-Manager-Tasks). Check there for details on each task.
# [5 Weeks] prior to the release:
- [ ]  __Call a meeting__
- [ ]  __Assign developers to the major subtasks__
    - Manager:
    - Lecture:
    - Doxygen:
    - Handbook:
    - Course:
    - Examples:
    - Website:
- [ ]  __Create a group (dumux-repositories) milestone in GitLab__
- [ ]  __Fix planned [Dune] and compiler compatibility__
- [ ]  __Go through all existing Gitlab tasks ([issues] or [MRs])__
- [ ]  __Assign the open [issues] and [MRs] to the developers__
- [ ]  __Open a milestone issue__
- [ ]  __Post a release schedule announcement on the mailing list__

# [3 weeks] prior to the release:
- [ ]  __Check all open MRs and Issues for severity and impact__
- [ ]  __Update the milestone issue__
- [ ]  __Announce soft feature freeze on the mailing list__
- [ ]  __Check in with the managers of the sub-tasks__

# [2 weeks] prior to the release:
__Hard Feature Freeze!__
- [ ]  __Check remaining MRs and Issues__
- [ ]  __Update the milestone issue__
- [ ]  __Update CHANGELOG__
- [ ]  __Announce hard feature freeze on the mailing list__
- [ ]  __Update all install scripts and the install text in the handbook__
- [ ]  __Header check__ (run `make headercheck`)
- [ ]  __Check and update the runtime parameters (see `bin/doc/getparameterlist.py`)
- [ ]  __Make sure the CMakeLists.txt in `dumux` subfolder are up-to-date__ (generated with `bin/utils/create_cmakelists.py`, dumux `make install` should result in a useable installed dumux version)
- [ ]  __Create release branches__
- [ ]  __Configure CI to test release branch__
- [ ]  __Local Testing__ (different compilers and dependency setups)

# During the week of the release:
__Final testing!__
- [ ]  __Re-run tests__
- [ ]  __Update copyright__
- [ ]  __Update License__
- [ ]  __Test Lecture and Course__
- [ ]  __Update Website__
- [ ]  __Create a release candidate__
- [ ]  __Call for testing__

# Releasing DuMu<sup>x</sup>
__Important__ These steps are normally done together with Bernd. Make sure to schedule the time properly.
- [ ]  __Protect release branch (in GitLab)__
- [ ]  __Create new tags__
- [ ]  __Prepare a Zenodo citation__
- [ ]  __Include the Zenodo citation to the website__
- [ ]  __Write a release email to dumux@listserv.uni-stuttgart.de__

# After the release
- [ ]  __Bump version in dune.module to next version__
- [ ]  __Sent an E-Mail to the dumux mailing list about updates on supported features and upcoming changes__
