#!/usr/bin/python3

import string
import argparse
import subprocess
from os.path import join, dirname, abspath


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
    return subprocess.run(cmd, capture_output=True, text=True).stdout


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process markdown files to be consumed by Doxygen")
    parser.add_argument("-f", "--file", required=True, help="The markdown file to be processed")
    parser.add_argument(
        "-p",
        "--prepend-header",
        required=False,
        help="Adds the given header at the top of the file"
    )
    args = vars(parser.parse_args())

    # Invoke math filter to translate GitLab flavoured math to such supported by Doxygen
    math_conversion_script = join(dirname(abspath(__file__)), "markdown-math-filter.pl")
    result = _invoke_and_retrieve_output(["perl", "-0777", "-p", math_conversion_script, args["file"]])

    # (maybe) prepend the given header
    if args["prepend_header"]:
        result = f"# {args['prepend_header']}\n\n{result}"

    # Give all headers anchors (labels) s.t. doxygen cross-linking works
    # correctly (may be fixed in the most recent Doxygen version)
    result_lines = []
    for line in result.split("\n"):
        result_lines.append(_add_header_label(line))

    # Print the final result for Doxygen to pick it up
    print("\n".join(result_lines))
