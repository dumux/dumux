# Developing DuMux

## Communicate with DuMux Developers

### Issues and Bug Tracking
The bug-tracking system via [GitLab Issues](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/issues)
offers the possibility to report bugs, discuss new development requests, or ask questions.
Feel free to register (if you don't have an account yet) and to contribute by opening an issue.

### Commits, Merges requests
To be up-to-date with the latest changes made to any git-repository, you can use RSS Feeds.
Simply navigate to [Issues](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/issues)
or [Activity](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/activity)
and hit the RSS subscribe button on the top right.

### Automatic Testing
The automatic testing using `Gitlab-CI` helps to constantly check the
DuMux problems for compiling and running correctly. An overview over
[recently run pipelines is available on GitLab](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/pipelines).
`Gitlab-CI` is configured via the `.gitlab-ci.yml` file in the project root directory.

### Mailing List
If you have questions, specific problems (which you really struggle to solve on your own),
or hints for the DuMux-developers, please contact the mailing list `dumux@iws.uni-stuttgart.de`.
You can [subscribe to the mailing list](https://listserv.uni-stuttgart.de/mailman/listinfo/dumux)
to be informed about upcoming releases, events, and see other peoples questions and answers.

### Coding Guidelines
Writing code in a readable manner is very important, especially
for future code developers (e.g. for adding features, debugging, etc.).
Therefore, we have a [contribution guide](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/CONTRIBUTING.md)
iuncluding a [style guide](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/doc/styleguide.md)
and tips on how to contribute effectively.

## Various tips and tricks
In the following, we mention a couple of tips & tricks. There is also
a [small Wiki](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/wikis/home) containing more.

### Optimized computation vs debugging
Dune and DuMux are built with the help of `dunecontrol`.
Per default, DuMux is compiled using optimization options, which leads to faster runtimes but is unsuitable
for debugging. For debug opts you can set `CMAKE_BUILD_TYPE` to `Debug` or `RelWithDebInfo`
in your options file. You can also do this in any of the `CMakeLists.txt` in DuMux by adding the line

```cmake
set(CMAKE_BUILD_TYPE Debug)
```

Afterwards rerun `cmake <path-to-build-dir>`.

### Dunecontrol for selected modules
A complete build using `dunecontrol` takes some time. In many cases not all modules need to be re-built.
Pass the flag `--only=dumux` to `dunecontrol` for configuring or building only DuMux. A more
complex example would be a case in which you have to configure and build only e.g. `dune-grid`
and DuMux. This is achieved by adding `--only=dune-grid,dumux`.

### Using Dune Debug Streams
Dune provides a helpful feature for keeping your debug-output organized.
It uses simple streams like `std::cout`, but they can be switched on and off
for the whole project. You can choose five different levels of severity:

```
5 - grave (dgrave)
4 - warning (dwarn)
3 - info (dinfo)
2 - verbose (dverb)
1 - very verbose (dvverb)
```

They are used as shown in the following example

```cpp
#include <dune/common/stdstreams.hh>
// define the minimal debug level somewhere in your code
#define DUNE_MINIMAL_DEBUG_LEVEL 4
Dune::dgrave << "message"; // will be printed
Dune::dwarn << "message"; // will be printed
Dune::dinfo << "message"; // will NOT be printed
```

### Make headercheck

When developing C++ it is important to always include what you use and not to rely on transitive includes.
To check one header file for all necessary includes to compile the contained code, use `make headercheck`.
Include the option `-DENABLE_HEADERCHECK=1` in your opts file and run `dunecontrol`.
Then go to the top level in your build-directory and type `make headercheck` to check all headers
or use bash auto-completion to find the target to build.
