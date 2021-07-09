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
from makeinstallscript import makeInstallScript, supportedLanguages
from makeinstallscript import getScriptExtension
from makeinstallscript import filterDependencies, addDependencyVersions
from makeinstallscript import addDependencyPatches

try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))
    from getmoduleinfo import getDependencies
    from common import callFromPath, runCommand, userQuery, query_yes_no
    from common import get_included_project_headers_cpp, findMatchingFiles
except Exception:
    sys.exit('Could not import common modul or getModuleInfo')


def get_readme_file_name():
    return "README.md"


def make_string_list(items, indentation="   "):
    return os.linesep.join([(indentation + '- ' + str(it)) for it in items])


def replace_file_content(file, new_content):
    with open(file, 'w') as new_file:
        new_file.write(new_content)


def append_file_content(file, content):
    with open(file, "a") as new_file:
        new_file.write(content)


def is_in_subtree(file, base):
    return Path(base).resolve() in Path(file).resolve().parents


def remove_redundant_folders(folders):
    return [
        sf for sf in folders
        if not any(is_in_subtree(sf, base) for base in folders)
    ]


def check_module_folder(mod_dir):
    if not os.path.isdir(mod_dir):
        raise Exception(
            f"Module folder {mod_dir} not found. "
            f"Make sure to run this script from one level above {mod_dir}."
        )


def check_sub_folders(mod_dir, sub_dirs):
    for sub_dir in sub_dirs:
        path = Path(mod_dir) / Path(sub_dir)
        errMsg = 'Cannot handle the given folder {}'.format(str(path))
        if not path.exists():
            raise Exception(errMsg + ' because it does not exist.')
        if not path.is_dir():
            raise Exception(errMsg + ' because it is not a directory.')


def extract_source_files(mod_dir, subfolders):
    sources = []
    for folder in subfolders:
        folderPath = os.path.abspath(os.path.join(mod_dir, folder))
        curSources = findMatchingFiles(folderPath, "*.cc")
        sources += [os.path.join(folderPath, s) for s in curSources]

    if not sources:
        raise Exception(
            "No sources found in the provided subfolders.",
            f" Run '{os.path.abspath(__file__)} --help' for details."
        )

    return sources


def run_dune_project():
    dune_project = shutil.which('duneproject', path="dune-common/bin")
    if dune_project is None:
        raise IOError("Could not find 'duneproject' in dune-common/bin.")
    subprocess.call([dune_project])


def detect_new_module():
    print("\nDetecting the newly created module")
    new_mod = max(
        [d for d in os.listdir() if os.path.isdir(d)],
        key=os.path.getmtime
    )
    print(f"Found module {new_mod}")

    new_mod_file = f"{new_mod}/dune.module"
    if not os.path.exists(new_mod_file):
        raise Exception(
            f"Could not find module file {new_mod_file}"
        )
    print(f"Successfully found the newly created module {new_mod}")
    return new_mod


def copy_sub_folders(subfolders, old_path, new_path):
    for b in subfolders:
        path_in_old_module = os.path.join(old_path, b)
        path_in_new_module = os.path.join(new_path, b)
        copy_tree(path_in_old_module, path_in_new_module)


def add_folders_to_cmake_lists(mod_path, folders):
    cml_file = os.path.join(mod_path, "CMakeLists.txt")
    with open(cml_file, "r") as cml:
        lines_bw = [line for line in reversed(cml.readlines())]

    idx = -1
    for line in lines_bw:
        if line.startswith("add_subdirectory"):
            idx = lines_bw.index(line)
            break

    new_lines = lines_bw[0:idx]
    new_lines += [f"add_subdirectory({b})\n" for b in reversed(subfolders)]
    new_lines += lines_bw[idx:-1]
    new_content = "".join(line for line in reversed(new_lines))

    replace_file_content(cml_file, new_content)


def find_headers(mod_path, source_files):
    with mp.Pool() as p:
        headers = itertools.chain.from_iterable(p.map(
            partial(get_included_project_headers_cpp,
                    project_base=mod_path),
            source_files
        ))
    return list(set(headers))


def copy_files(files, old_path, new_path):
    new_files = []
    for file in files:
        file = os.path.realpath(file)
        new_file = file.replace(old_path, new_path, 1)
        new_files.append(new_file)
        if not os.path.exists(new_file):  # exist_ok below fails (due to mode?)
            os.makedirs(new_file, exist_ok=True)
            shutil.copy(file, new_file)
    return new_files


def find_folders_without_source_files(mod_path, check_subfolders, sources):
    source_dirs = [os.path.dirname(s) for s in sources]
    source_dirs = list(set(source_dirs))

    def has_child_source_dir(dir):
        return any(is_in_subtree(s, dir) for s in source_dirs)

    def is_no_source_dir(dir):
        return dir not in source_dirs and not has_child_source_dir(dir)

    no_source_dirs = []
    for sf in check_subfolders:
        for root, dirs, _ in os.walk(os.path.join(mod_path, sf)):
            for dir in dirs:
                dir_path = os.path.join(root, dir)
                if is_no_source_dir(dir_path):
                    no_source_dirs.append(dir_path)
    no_source_dirs = remove_redundant_folders(list(set(no_source_dirs)))

    def remove_empty_parents():
        folder_map = {}
        for f in no_source_dirs:
            parent = os.path.dirname(f)
            if parent not in folder_map:
                folder_map[parent] = []
            folder_map[parent].append(f)

        for parent in folder_map:
            found = set(folder_map[parent])
            for root, dirs, files in os.walk(parent):
                dirs = [os.path.join(root, d) for d in dirs]
                if set(dirs) == found and is_no_source_dir(parent):
                    for entry in found:
                        no_source_dirs.remove(entry)
                    no_source_dirs.append(parent)
        return no_source_dirs
    return remove_empty_parents()


def remove_folder(mod_path, subfolder):

    def remove_add_subdir_cmd(cmake_lists, folder):
        with open(cmake_lists, "r") as cml:
            content = cml.read()

        key = f"add_subdirectory({folder})\n"
        if key in content:
            replace_file_content(cmake_lists, content.replace(key, ""))
            return True
        return False

    subfolderpath = os.path.abspath(os.path.join(mod_path, subfolder))
    subfoldername = os.path.basename(subfolderpath.rstrip(os.sep))

    main_cml = os.path.join(mod_path, "CMakeLists.txt")
    parent_cml = os.path.join(subfolderpath, '../CMakeLists.txt')
    if not remove_add_subdir_cmd(main_cml, subfolder):
        if os.path.exists(parent_cml):
            if not remove_add_subdir_cmd(parent_cml, subfoldername):
                raise Exception(
                    "Could not find folder {subfolderpath} in CMakeLists.txt"
                )
    shutil.rmtree(subfolderpath)


def guide_folder_deletion(mod_path, candidates):
    candidate_list = make_string_list(candidates)
    print(
        "\n"
        "Could not automatically determine if the following directories\n"
        "contain data that is essential for the extracted applications:\n"
        "\n"
        f"{candidate_list}\n"
    )

    do_remove = query_yes_no(
        "Do you want to remove some of them "
        "(by choosing 'no' they are all preserved)?\n",
        default="no"
    )

    deleted = []
    if do_remove:
        for folder in folders_without_sources:
            do_remove_folder = query_yes_no(
                f"Do you want to delete the folder {folder}?", default="yes"
            )

            if do_remove_folder:
                rel_path = os.path.relpath(folder, mod_path)
                remove_folder(new_module_path, rel_path)
                deleted.append(rel_path)
    return deleted


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
            print("Checking the repo (you may have to introcude credentials):")
            remote_content = runCommand(
                'git ls-remote {}'.format(remote),
                suppressTraceBack=True
            )
        except Exception:
            print(" - Could not find your repo at {}. ".format(remote))
            print(" - Please revisit the provided information.")
            continue

        if (remote_content == ''):
            return remote
        raise Exception("Remote reposity is not empty.")


def run_git_cmd(path, cmd):
    print(f"Running {cmd} (you may have to provide your credentials)")
    callFromPath(path)(runCommand)(cmd)


def push_repository(mod_path, remote_url):
    run_git_cmd(mod_path, 'git push -u origin master')


def guide_repository_initialization(mod_path):

    has_repo = query_yes_no(
        "Do you have an empty remote repository to push the code to?\n"
        "(If not, we recommend that you create one and answer with 'yes'.)"
    )
    remote_url = None if not has_repo else query_empty_remote_repo_url()

    run_git_cmd(mod_path, 'git init')
    run_git_cmd(mod_path, 'git add .')
    run_git_cmd(mod_path, 'git commit -m "Initial commit"')

    if has_repo:
        run_git_cmd(mod_path, 'git remote add origin {}'.format(remote_url))
        push_repository(mod_path, remote_url)

    return remote_url


def guide_versions_in_readme(mod_path, readme=None):

    write_version_info = query_yes_no(
        "Write detailed version information"
        f" (folder/branch/commits/dates) into {get_readme_file_name()}?\n"
    )
    if write_version_info:
        try:
            deps = getDependencies(mod_path)
            versions = getPersistentVersions(
                [dep['folder'] for dep in deps],
                ignoreUntracked=True
            )
        except Exception as e:
            sys.exit(f"Error when determining version info: {e}")

        if not readme:
            readme = os.path.join(mod_path, get_readme_file_name())

        append_file_content(
            readme, "\n## Version Information\n\n" + versionTable(versions)
        )


def guide_install_script_generation(mod_path, script_name_body, skip=[]):
    language = userQuery(
        'In which language would you like to generate the install script?',
        supportedLanguages()
    )
    ext = getScriptExtension(language)
    inst_script_name = script_name_body + ext

    def getModDependencies():
        try:
            return getDependencies(mod_path)
        except Exception as e:
            raise Exception(f"Error when determining dependencies: {e}")

    def processDependencies(deps):
        try:
            deps = filterDependencies(deps, skip+[mod_path])
            deps = addDependencyVersions(deps, ignoreUntracked=True)
            deps = addDependencyPatches(deps)
            return deps
        except Exception as e:
            raise Exception(f"Error processing the dependencies: {e}.")

    def makeScript(deps):
        try:
            makeInstallScript(
                modPath=mod_path,
                dependencies=deps,
                scriptName=inst_script_name,
                language=language,
                topFolderName=''
            )
        except Exception as e:
            raise Exception(f"Error during install script generation: {e}")

    deps = getModDependencies()
    deps = processDependencies(deps)
    if not deps:
        print("No dependencies found. Skipping install script generation.")
        return ''

    makeScript(deps)
    return inst_script_name


def process_install_script(script_name, mod_path, remote_url):

    new_script = os.path.join(mod_path, script_name)
    shutil.move(script_name, new_script)
    subprocess.call(['chmod', 'u+x', new_script])

    run_git_cmd(mod_path, f'git add {script_name}')
    run_git_cmd(mod_path, 'git commit -m "add install script"')

    readme = os.path.join(mod_path, get_readme_file_name())
    append_file_content(
        readme,
        info_readme_installation(remote_url, script_name, new_module_name)
    )
    run_git_cmd(mod_path, f'git commit -m "update readme" {readme}')

    if remote_url:
        push_repository(mod_path, remote_url)
    else:
        print(
            "\n"
            "Please remember to manually fix the installation instructions\n"
            "once you know the remote URL where your repository will be hosted"
        )


def no_remote_url_warning(new_mod, old_mod):
    return f"""

Warning: no remote registered for new module {new_mod}. We will therefore use
         the old module {old_mod} to determine dependencies. Please make sure
         that the dependencies listed in the dune.module files match (you
         should have been asked for the dependencies earlier when creating
         of the new module).

"""


# Some general information for users of the script
def info_initial(module_dir, module_path, subfolders, source_files):
    sources_list = make_string_list(source_files)
    subfolder_list = make_string_list(
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


# Main part of readme
def info_readme_main(module_dir, subfolders, source_files):
    def rel_path(p):
        return os.path.relpath(p, module_dir)

    subfolders_str = ''.join(
        [f"*   `{rel_path(d)}`\n" for d in subfolders]
    )
    sources_str = ''.join(
        [f"*   `{rel_path(s)}`\n" for s in source_files]
    )

    return f"""
This file has been created automatically. Please adapt it to your needs.

## Content

The content of this DUNE module was extracted from the module `{module_dir}`.
In particular, the following subfolders of `{module_dir}` have been extracted:
{subfolders_str}
Additionally, all headers in `{module_dir}` that are required to build the
executables from the sources
{sources_str}
have been extracted. You can configure the module just like any other DUNE
module by using `dunecontrol`. For building and running the executables,
please go to the build folders corresponding to the sources listed above.\n
"""


# Installation part of readme
def info_readme_installation(remote_url, install_script_name, new_module_name):
    install_script_path = os.path.join(new_module_name, install_script_name)
    remote_hints = ''
    if not remote_url:
        remote_url = "$remoteurl"
        remote_hints = "\nImportant: $remoteurl has to be corrected!\n"
    return f"""

## Installation

The installation procedure is done as follows :
Create a root folder, e.g. `DUMUX`, enter the previously created folder,
clone the remote and use the install script `{install_script_name}`
provided in this repository to install all dependent modules.
{remote_hints}
```sh
mkdir DUMUX
cd DUMUX
git clone {remote_url}
./{install_script_path}
```

This will clone all modules into the directory `DUMUX`,
configure your module with `dunecontrol` and build tests.

"""


# Some final remarks
def info_final(new_mod_name):
    return f"""
========================================================================

The module was extracted successfully!

The extracted module is contained in the subfolder '{new_mod_name}'.
You can configure it with
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts —only={new_mod_name} all
"""


###############
# Main script #
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

    # prepare input
    args = vars(parser.parse_args())
    module_dir = args['module_dir'].strip(os.sep)
    module_path = os.path.abspath(module_dir)
    base_folder = os.path.abspath(os.path.join(module_dir, '../'))

    # find executable applications
    subfolders = remove_redundant_folders(list(set(args['subfolder'])))
    check_module_folder(module_dir)
    check_sub_folders(module_dir, subfolders)
    source_files = extract_source_files(module_dir, subfolders)

    # guide user through new module creation
    print(info_initial(module_dir, module_path, subfolders, source_files))
    input("Please read the above carefully and press [Enter] to proceed...")

    run_dune_project()
    new_module_name = detect_new_module()
    new_module_path = os.path.join(base_folder, new_module_name)

    # prepare all data in new module
    copy_sub_folders(subfolders, module_path, new_module_path)
    add_folders_to_cmake_lists(new_module_path, subfolders)

    headers = find_headers(module_path, source_files)
    new_headers = copy_files(headers, module_path, new_module_path)
    new_source_files = copy_files(source_files, module_path, new_module_path)
    new_git_ignore = copy_files(["dumux/.gitignore"], "dumux", new_module_path)

    # guide user through deletion of possibly unneccessary folders
    folders_without_sources = find_folders_without_source_files(
        new_module_path, subfolders, new_headers + new_source_files
    )

    actual_subfolders = subfolders
    if folders_without_sources:
        deleted_folders = guide_folder_deletion(
            new_module_path, folders_without_sources
        )
        actual_subfolders = [s for s in subfolders if s not in deleted_folders]

    # remove stuff that is created when running duneproject
    if os.path.join(new_module_path, 'dune') not in folders_without_sources:
        remove_folder(new_module_path, 'dune')
    if os.path.join(new_module_path, 'src') not in folders_without_sources:
        remove_folder(new_module_path, 'src')

    # prepare new readme file
    dune_readme = os.path.join(new_module_path, "README")
    if os.path.exists(dune_readme):
        os.remove(dune_readme)

    new_readme = os.path.join(new_module_path, get_readme_file_name())
    replace_file_content(
        new_readme,
        info_readme_main(module_dir, actual_subfolders, source_files)
    )

    # try to initialize repo (to use its url in later steps)
    remote_url = guide_repository_initialization(new_module_path)
    if not remote_url:
        print(no_remote_url_warning(new_module_name, module_dir))

    # make install script & finalize readme
    check_path = new_module_path if remote_url else module_path
    skip_path = new_module_path if not remote_url else module_path

    guide_versions_in_readme(check_path, new_readme)
    iscript = 'install_' + new_module_name
    iscript = guide_install_script_generation(check_path, iscript, [skip_path])
    if iscript:
        process_install_script(iscript, new_module_path, remote_url)

    print(info_final(new_module_name))
