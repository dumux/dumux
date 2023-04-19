<!--
SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
SPDX-License-Identifier: CC0-1.0
-->

<!--
Thanks for considering to open a merge request!  
Before asking for a review of your MR, please read the [contributing guidelines](/CONTRIBUTING.md)
-->
**What this MR does / why does DuMux need it**:

TODO: insert text here

<!--
Is there a corresponding issue? Add "Fixes hashtag issuenumber" which will automatically close the issue when this MR is merged. Add "Related to hashtag issuenumber" if it's related but doesn't fix the issue completely.
-->

**Notes for the reviewer**

TODO: insert text here

<!--
Keep the following TODO list in the merge request description for documentation.
Bullet points marked with _(if not applicable remove)_ may be removed.
-->

Before you request a review from someone, make sure to revise the following points:

- [ ] does the new code follow the [style guide](doc/styleguide.md)?
- [ ] do the test pipelines pass? (see guide on [how to run pipelines for a merge request](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/wikis/Running-test-pipelines-for-merge-requests))
- [ ] is the code you changed and/or the new code you wrote covered in the test suite? (if not, extend the existing tests or write new ones)
- [ ] does your change affect public interfaces or behavior, or, does it introduce a new feature? If so, document the change in `CHANGELOG.md`.
- [ ] is the list of the header includes complete? ("include what you use")
- [ ] all files have to end with a `\n` character. Make sure there is no `\ No newline at end of file` comment in "Changes" of this MR.
- [ ] (if not applicable remove) are newly introduced or modified physical values/functions backed up with a scientific reference (including doi) in the docs?
- [ ] (if not applicable remove) if the examples are modified, is the documentation regenerated (using [`generate_example_docs.py`](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/generate_example_docs.py))

<!--
The following aspects might also come up during review:

* Does the change reduce the performance of the code (more CPU time or more memory) and is this justified by the benefits
* Does the change improve the performance? (if yes, add this aspect to the MR description)
* Is the code is a gross violation of programming best practices such as DRY (don't repeat yourself / code duplication, see https://de.wikipedia.org/wiki/Don%E2%80%99t_repeat_yourself, the SOLID principles (https://en.wikipedia.org/wiki/SOLID), or the C++ Core Guidelines (https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)?
* Is the code well-documented, concise, easily readable? (e.g. variables are well-named, the logic is split into small & well-named functions)
-->
