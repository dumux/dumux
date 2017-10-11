# Exercise #4 (DuMuX course)

This exercise describes how to create a new DuMuX module 
and how to create a corresponding GitLab project.

This is the suggested
workflow to develop code on top of DuMuX. 

### Task 1: Create new dune module
<hr>

* Execute the following command (bash environment) in the top-folder, i.e. above the dumux folder

```bash
./dune-common/bin/duneproject
```

* Follow the introductions and specify
    * as name of the new module: `dumux-example`
    * as module dependencies: `dumux`
    * a version at your choice
    * your email address

<hr><br><br>
### Task 2: Rerun dunecontrol to configure your new project 
<hr>

The following command will configure your new module

```bash
./dune-common/bin/dunecontrol --opts=<opts file> --only=dumux-example all
```

<hr><br><br>
### Task 3: Create a new test case within your new DuMuX module
<hr>

* Create a new folder (in your module folder), e.g. `appl`

```bash
mkdir appl
```

* Copy some test case from the dumux module, e.g. test_box1p
    * Copy the problem, spatialparams, cc source file, input file

* Adjust the CMakeLists.txt file to include your new subdirectory

* Add a new CMakeLists.txt in the folder `appl` with the content

```cmake
# add a new box 1p test
dune_add_test(NAME test_box1p
              SOURCES test_box1p.cc)

# link the input file to the build folder
dune_symlink_to_source_files(FILES test_box1p.input)
 
```

* Reconfigure your module by running in the topmost directory of your new module

```bash
cmake build-cmake
```

* Build and execute the test problem

```bash
cd build-cmake
make build_tests
cd appl
./test_box1p
```

<hr><br><br>
### Task 4: Create a new GitLab project
<hr>

* Login with your username and password at https://git.iws.uni-stuttgart.de/

Note: If you don't have an account create one. We allow anyone to host repositories
on our GitLab instance as long as it is DuMuX related.

* Click the **New project** button

* Specify your project name and click the **Create project** button

* Follow the given instructions for an *existing folder* 

**Important**: Before executing the `git add .` command, you should add your cmake build folder to `.gitignore`. 
The easiest way to do so is to copy the `.gitignore` file from the dumux module into your module path. If everything
worked, executing `git status` should not show `build-cmake` anymore. Never put your executables or other build files
under version control. Only source files (`*.hh`, `*.cc`, `*.input`, `CMakeLists.txt`) should be under version control.
