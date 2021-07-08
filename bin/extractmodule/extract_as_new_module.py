#!/usr/bin/env python3
"""
This script extracts some specified applications into a separate Dune module.
For example make a dumux-pub repository accompanying a scientific paper.
"""

import sys
import os
import subprocess
import os.path
import argparse
import shutil
from distutils.dir_util import copy_tree
import multiprocessing as mp
import itertools
from pathlib import Path
from functools import partial
from util import getPersistentVersions
from util import versionTable
from makeinstallscript import makeInstallScript
from makeinstallscript import python_or_bash
try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))
    from getmoduleinfo import getDependencies
    from common import callFromPath, runCommand, query_yes_no
    from common import get_included_project_headers_cpp, findMatchingFiles
except Exception:
    sys.exit('Could not import common modul or getModuleInfo')


def get_dumux_pub_project_url(project_name):
    base_url = 'https://git.iws.uni-stuttgart.de/dumux-pub/'
    return base_url + '{}.git'.format(project_name.lower())


def query_empty_remote_repo_url():
    while True:
        if query_yes_no("Is your repository hosted in dumux-pub?"):
            project = input("Please provide the name of your project "
                            "(usually AuthorLastNameYearX, e.g. Luigi2020b)\n")
            remote = get_dumux_pub_project_url(project)
        else:
            remote = input("Please provide the URL of your repository:\n")

        try:
            remote_content = runCommand(
                'git ls-remote {}'.format(remote),
                suppressTraceBack=True
            )
        except Exception:
            print("Could not find your repo at {}. ".format(remote))
            print("Please double-check and reenter correct information")
            continue

        if (remote_content == ''):
            return remote
        raise Exception("Remote reposity is not empty.")


def check_input_folders(mod_dir, sub_dirs):
    if not os.path.isdir(mod_dir):
        raise Exception(
            f"Module folder {mod_dir} not found. "
            f"Make sure to run this script from one level above {mod_dir}."
        )

    for sub_dir in sub_dirs:
        path = Path(mod_dir) / Path(sub_dir)
        errMsg = 'Cannot handle the given folder {}'.format(str(path))
        if not path.exists():
            raise Exception(errMsg + ' because it does not exist.')
        if not path.is_dir():
            raise Exception(errMsg + ' because it is not a directory.')


def is_in_subtree(file, base):
    return Path(base).resolve() in Path(file).resolve().parents


def remove_redundant_subfolders(subfolders):
    return [
        sf for sf in subfolders
        if not any(is_in_subtree(sf, base) for base in subfolders)
    ]


def extract_source_files(module_dir, subfolders):
    sources = []
    for folder in subfolders:
        folderPath = os.path.abspath(os.path.join(module_dir, folder))
        curSources = findMatchingFiles(folderPath, "*.cc")
        sources += [os.path.join(folderPath, s) for s in curSources]
    return sources


def check_module(module_name):
    mod_file = f"{module_name}/dune.module"
    if not Path(mod_file).exists():
        raise Exception(
            f"Could not find module file in {mod_file}"
        )


def add_folders_to_cmake_lists(cml_file, folders):
    with open(cml_file) as cml:
        lines_bw = [line for line in reversed(cml.readlines())]

        idx = -1
        for line in lines_bw:
            if line.startswith("add_subdirectory"):
                idx = lines_bw.index(line)
                break

        new_lines = lines_bw[0:idx]
        new_lines += [f"add_subdirectory({b})\n" for b in reversed(subfolders)]
        new_lines += lines_bw[idx:-1]
        return "".join(line for line in reversed(new_lines))


def make_list_string(items, indentation="   "):
    return os.linesep.join([indentation + i for i in items])


###################################################################
# Some general information for users of the script
###################################################################
def info_explanations(module_dir, module_path, subfolders, source_files):
    sources_list = make_list_string(source_files)
    subfolder_list = make_list_string(
        [os.path.relpath(f, module_path) for f in subfolders]
    )

    return f"""
This script will extract the following subfolders of
the module '{module_dir}':
{subfolder_list}

and all headers contained in '{module_dir}'
that are required to build the exectutables from the sources:
{sources_list}

The extracted files are copied into a new DUNE module retaining the directory
structure. The files required for creating a working DUNE module (such as
CMakeLists.txt) will be created and/or updated.

In the next step, the 'duneproject' script will be executed to guide the
creation of your new DUNE module. Please answer all upcoming queries to the
best of your knowledge.

Important: the new module should NOT depend on the module '{module_dir}'
"""


###################################################################
# Main part of README.md
###################################################################
def info_readme_main(module_dir, subfolders, source_files):
    subfolders_str = ''.join([f"*   `{d}`\n" for d in subfolders])
    sources_str = ''.join([f"*   `{s}`\n" for s in source_files])
    return f"""
This file has been created automatically. Please adapt it to your needs.

## Content

The content of this DUNE module has been extracted from the module `{module_dir}`.
In particular, the following subfolders of `{module_dir}` have been extracted:
{subfolders_str}
Additionally, all headers in `{module_dir}` that are required to build the
executables from the sources
{sources_str}
have been extracted. You can configure the module just like any other DUNE
module by using `dunecontrol`. For building and running the executables,
please go to the build folders corresponding to the sources listed above.\n
"""


###################################################################
# Installation part of README.md
###################################################################
def info_readme_installation(remoteurl, install_script_name, new_module_name):
    install_script_path = os.path.join(new_module_name, install_script_name)
    return f"""

## Installation

The installation procedure is done as follows :
Create a root folder, e.g. `DUMUX`, enter the previously created folder,
clone the remote and use the install script `{install_script_name}`
provided in this repository to install all dependent modules.

```sh
mkdir DUMUX
cd DUMUX
git clone {remoteurl}
chmod u+x {install_script_path}
./{install_script_path}
```

This will clone all modules into the directory `DUMUX`,
configure your module with `dunecontrol` and build tests.

"""


###################################################################
# Infos on how to create the install script manually
###################################################################
def info_make_install(new_module_name):
    return f"""
========================================================================

The extracted module is contained in the subfolder '{new_module_name}'.
You can configure it with
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts —only={new_module_name} all
"""


if __name__ == "__main__":

    # set script parameters
    epilog = '''
    -----------------------------------------------------------
    The script has to be called one level above module_dir.
    At least one of the subfolders (FOLDER_1 [FOLDER_2 ...]) has
    to contain a source file *.cc of an executable for which
    you would like to timber a table in dumux-pub.)
    -----------------------------------------------------------

    Example usage:
    ./dumux/bin/extractmodule/extract_as_new_module.py dumux-fracture appl test

    (extracts the subfolders appl and test from the module dumux-fracture)

    '''

    parser = argparse.ArgumentParser(
        prog='extract_as_new_module.py',
        usage='./dumux/bin/extractmodule/extract_as_new_module.py'
              ' module_dir SUBFOLDER_1 [SUBFOLDER_2 ...]',
        description='This script extracts subfolders of a given DUNE module'
                    ' into a new DUNE module.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog)
    parser.add_argument('module_dir',
                        help='Module from which the subfolder is extracted')
    parser.add_argument('subfolder', nargs='+',
                        help='subfolder(s) of "module_dir" to be extracted')

    args = vars(parser.parse_args())
    module_dir = args['module_dir'].strip(os.sep)
    module_path = os.path.abspath(module_dir)
    subfolders = remove_redundant_subfolders(list(set(args['subfolder'])))

    try:
        check_input_folders(module_dir, subfolders)
    except Exception as e:
        sys.exit(f"Error when checking input folders: {e}")

    source_files = extract_source_files(module_dir, subfolders)
    if not source_files:
        sys.exit("No sources found in the provided subfolders.",
                 f" Run '{os.path.abspath(__file__)} --help' for details.")

    dune_project = shutil.which('duneproject', path="dune-common/bin")
    if dune_project is None:
        sys.exit("Could not find 'duneproject' in dune-common/bin.")

    print(info_explanations(
        module_dir, module_path, subfolders, source_files
    ))
    input("Please read the above carefully and press [Enter] to proceed...")

    try:
        subprocess.call([dune_project])
    except Exception:
        sys.exit("Failed to generate new module with duneproject!")

    try:
        new_module_name = max(
            [d for d in os.listdir() if os.path.isdir(d)],
            key=os.path.getmtime
        )
        check_module(new_module_name)
    except Exception as e:
        sys.exit(f"Validity check on new module {new_module_name} returned error: {e}")

    print(
        "\n-- duneproject done --\n"
        f"\nFound new module {new_module_name}\n"
        "Copying source files..."
    )

    new_module_path = os.path.join(
        os.path.abspath(os.path.join(module_dir, '../')),
        new_module_name
    )

    try:
        for b in subfolders:
            path_in_old_module = os.path.join(module_path, b)
            path_in_new_module = os.path.join(new_module_path, b)
            copy_tree(path_in_old_module, path_in_new_module)
    except Exception as e:
        sys.exit(f"Error copying the subfolders to new module: {e}")

    try:
        new_cmake_lists_ = add_folders_to_cmake_lists(
            os.path.join(new_module_path, "CMakeLists.txt"),
            subfolders
        )

        with open(os.path.join(new_module_path, "CMakeLists.txt"), "w") as cml:
            cml.write(new_cmake_lists_)
    except Exception as e:
        sys.exit(f"Couldn't create the new CMakeLists.txt file, error: {e}")

    try:
        with mp.Pool() as p:
            headers = itertools.chain.from_iterable(p.map(
                partial(get_included_project_headers_cpp,
                        project_base=module_path),
                source_files
            ))
        headers = list(set(headers))
    except Exception:
        sys.exit("Could not find the headers included by the sources")

    try:
        header_dirs = []
        for header in headers:
            header_dir = os.path.dirname(os.path.realpath(header))
            path_in_new_module = header_dir.replace(module_dir, new_module_name, 1)
            header_dirs.append(path_in_new_module)
            os.makedirs(path_in_new_module, exist_ok=True)
            shutil.copy(header, path_in_new_module)
    except Exception as e:
        sys.exit(f"Could not copy all required headers to the new module: {e}")

    dumux_gitignore_file = "dumux/.gitignore"
    shutil.copy(dumux_gitignore_file, new_module_path)


    try:
        new_source_files = [
            s.replace(module_dir, new_module_name, 1) for s in source_files
        ]
        for b in subfolders:
            def matching_path(path, list_of_paths, header_dirs):
                for p in list_of_paths:
                    if path in p or path in header_dirs:
                        return True
                return False

            no_source_folder = []
            for path, dirs, files in os.walk(os.path.join(new_module_path, b)):
                for d in dirs:
                    dir_path = os.path.join(path, d)
                    keep = matching_path(dir_path, new_source_files, header_dirs)
                    if not keep:
                        rel_dir_path = os.path.relpath(dir_path, new_module_path)
                        no_source_folder.append(rel_dir_path)

            if no_source_folder:
                no_source_lists = '\n'.join(dir for dir in no_source_folder)
                yes_to_all = query_yes_no(
                    "Could not automatically determine if the following directories "
                    "contain data essential for the extracted applications:\n{0}.\n"
                    "Do you want to copy all of them in the new module {1}?\n"
                    "(By choosing no you need to consider each folder seperately)"
                    .format(no_source_lists, new_module_name),
                    default="yes"
                )

                if not yes_to_all:
                    for rel_dir_path in no_source_folder:
                        answer_is_yes = query_yes_no(
                            f"Could not automatically determine if {rel_dir_path} contains data essential for the extracted applications.\n"
                            f"Do you want to copy the files in {rel_dir_path} to the extracted module {new_module_name}",
                            default="yes"
                        )
                        if not answer_is_yes:
                            # remove copy of directory
                            shutil.rmtree(os.path.join(new_module_path, rel_dir_path))
                            # remove entry from CMakeLists.txt
                            cml_path = os.path.join(module_path, "CMakeLists.txt")
                            with open(cml_path, "r+") as cml:
                                content = cml.read()
                                cml.seek(0)
                                content = content.replace(
                                    "add_subdirectory({d})\n", ""
                                )
                                cml.write(content)
                                cml.truncate()
    except Exception:
        sys.exit("Could not handle folders without sources")

    delete_folders = ['dune', 'src']
    try:
        with open(os.path.join(new_module_path, "CMakeLists.txt"), "r+") as cml:
            content = cml.read()
            cml.seek(0)
            for del_folder in delete_folders:
                if del_folder not in subfolders:
                    shutil.rmtree(os.path.join(new_module_path, del_folder))
                    content = content.replace("add_subdirectory(dune)\n", "")
            cml.write(content)
            cml.truncate()
    except Exception as e:
        sys.exit("Could not delete the folders: {}, error: {}"
                 .format(', '.join(delete_folders), e))

    os.remove(os.path.join(new_module_path, 'README'))
    readme_path = os.path.join(new_module_path, "README.md")
    with open(readme_path, "w") as readme_file:
        readme_file.write(
            info_readme_main(module_dir, subfolders, source_files)
        )

    if query_yes_no("Write detailed version information"
                    " (folder/branch/commits/dates) into README.md?\n"):
        print("Looking for the dune modules in path: "
              + os.path.abspath(".")
              + "...")

        versions = getPersistentVersions(
            [dep['folder'] for dep in getDependencies(module_path)], True
        )

        with open(readme_path, "a") as readme_file:
            readme_file.write(
                "\n## Version Information\n\n" +
                versionTable(versions)
            )

    language = python_or_bash()
    install_script_name = 'install_' + new_module_name + '.%s'%("sh" if language == "bash" else "py")
    try:
        makeInstallScript(new_module_path, ignoreUntracked=True, skipFolders=new_module_name,
                          suppressHints=True, topFolderName=None, language=language)
    except Exception:
        sys.exit(
            f"Automatically generate install script {install_script_name} failed."
            "\nTo create the script, use the 'makeinstallscript.py' script in the same folder."
            "\nRun 'python3 makeinstallscript.py --help' for more detailed information"
        )

    shutil.move(install_script_name,
                os.path.join(new_module_path, install_script_name))
    print(info_make_install(new_module_name))

    remote_url = "{$remoteurl$} (needs to be manuelly adapted later)"
    has_repo = query_yes_no(
        "Do you have an empty remote repository to push the code to (recommended)?"
    )

    if has_repo:
        try:
            remote_url = query_empty_remote_repo_url()
        except Exception as e:
            sys.exit(f"Error when trying to get the remote URL: {e}")

    with open(readme_path, "a") as readme_file:
        readme_file.write(info_readme_installation(
            remote_url, install_script_name, new_module_name
        ))

    def run_git_cmd(cmd):
        print(f"Running {cmd}")
        run_from_mod = callFromPath(new_module_path)(runCommand)
        run_from_mod(cmd)

    try:
        run_git_cmd('git init')
        run_git_cmd('git add .')
        run_git_cmd('git commit -m "Initial commit"')
        if has_repo:
            run_git_cmd('git remote add origin {}'.format(remote_url))
            run_git_cmd('git push -u origin master')
        else:
            print("\nPlease remember to replace '$remoteurl' in the installation\n"
                  "part of the README file with the correct URL of your remote git repository.")
    except Exception as e:
        sys.exit(f"Error when pushing the new module: {e}")
