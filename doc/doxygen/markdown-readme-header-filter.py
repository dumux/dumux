#!/usr/bin/python3

# Process the main readme (GitLab Markdown) to render properly with doxygen
# -> Prepend level-1 header (doxygen will use this as caption for the section)
#    and turn all level-1 headers into level-2 ones. These will then be sub-sections.
# -> Add labels to headers such that section links within the readme work.

import sys
import string


def _remove_punctuations(text: str) -> str:
    return "".join(filter(lambda c: c not in string.punctuation, text))


def _modify_header_level(line: str) -> str:
    is_level_one_header = line.startswith("# ")
    return f"#{line}" if is_level_one_header else f"{line}"


def _add_header_label(line: str) -> str:
    if not line.startswith("#"):
        return line
    line = line.rstrip("\n")
    label = _remove_punctuations(line)
    label = label.strip(" ").replace(" ", "-").lower()
    return f"{line} {{#{label}}}\n"


assert len(sys.argv) > 1
lines = ["# Introduction\n"]
with open(sys.argv[1]) as readme:
    for line in readme:
        lines.append(_add_header_label(_modify_header_level(line)))
print("".join(lines))
