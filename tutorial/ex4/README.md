# Exercise #4 (DuMuX course)

This exercise describes how to create a new DuMuX module 
and how to create a corresponding gitlab project. 

## Create new DuMuX module

### 1. Execute the following command (bash environment):
```bash
$ ./dune-common/bin/duneproject
```
Follow the introductions and specify
* Name of your module, e.g. dumux-ex
* Module dependency, which is dumux

### 2. Run dunecontrol 
here you can use `--only=Module-Name`


## Create a new test case within your new DuMuX module

### 1. Create a new folder (in your module folder), e.g. appl

### 2. Copy some test case from the dumux module, e.g. test_box1p

### 3. Incorporate this test case into your cmake files

### 4. Re-run **dunecontrol** 

### 5. Execute your test problem


## Create a new gitlab project

### 1. Login with your username and password (<https://git.iws.uni-stuttgart.de>)

### 2. Click the **New project** button

### 3. Follow the given instructions for an *existing folder* 

**Important**: Before executing the `git add .` command, you should add your cmake build folder to gitignore. 
The easiest way to do so is to just copy the *.gitignore* file from your dumux module into your module path. 
