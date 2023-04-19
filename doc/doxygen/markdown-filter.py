#!/usr/bin/python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


import sys
import string
import subprocess
from os.path import join, dirname, abspath, exists


THIS_DIR = dirname(__file__)
DOC_DIR = dirname(THIS_DIR)
TOP_LEVEL_DIR = dirname(DOC_DIR)
MAIN_README = abspath(join(TOP_LEVEL_DIR, "README.md"))
assert exists(MAIN_README)


def _filter_characters(text: str) -> str:
    return "".join(filter(lambda c: c not in string.punctuation and c not in string.digits, text))


def _add_header_label(line: str) -> str:
    if not line.startswith("#"):
        return line
    line = line.rstrip("\n")
    label = _filter_characters(line)
    label = label.strip(" ").replace(" ", "-").lower()
    return f"{line} {{#{label}}}\n"


def _invoke_and_retrieve_output(cmd) -> str:
    return subprocess.check_output(cmd, text=True)


def _is_main_readme_file(path) -> bool:
    return path == MAIN_README


def _enable_doxygen_only_content(line: str) -> str:
    if line.startswith("<!--DOXYGEN_ONLY "):
        return line.replace("<!--DOXYGEN_ONLY ", "").replace("-->", "")
    return line


def _change_verbatim(line: str) -> bool:
    return line.count("```") % 2 != 0


assert len(sys.argv[1]) > 1
filePath = abspath(sys.argv[1])

# Invoke math filter to translate GitLab flavoured math to such supported by Doxygen
math_conversion_script = join(dirname(abspath(__file__)), "markdown-math-filter.pl")
result = _invoke_and_retrieve_output(["perl", "-0777", "-p", math_conversion_script, filePath])

# for the main readme, prepend "Introduction" header
if _is_main_readme_file(filePath):
    result = f"# Introduction\n\n{result}"

# Give all headers anchors (labels) s.t. doxygen cross-linking works
# correctly (may be fixed in the most recent Doxygen version)
result_lines = []
verbatim_environment = False
for line in result.split("\n"):
    if _change_verbatim(line):
        verbatim_environment = not verbatim_environment
    result_lines.append(
        line if verbatim_environment else _enable_doxygen_only_content(_add_header_label(line))
    )

# Print the final result for Doxygen to pick it up
print("\n".join(result_lines))
