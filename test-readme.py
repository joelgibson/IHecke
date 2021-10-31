"""
This file extracts sections of code in the README.md file and tests them by running them through
an actual Magma session. The point of this is to keep the tutorial in the README up-to-date and
always valid.

Each test is delimited by the HTML comments
    <!-- BEGIN TEST {test_name} -->
    ...
    <!-- END TEST {test_name} -->
Any markdown code examples between these are read. Each must start with a command-line invocation
of magma, i.e '$ magma'. After that, lines starting with '> ' are considered as input, and other
lines are considered as output. All the input lines are extracted and then run through a Magma
session in one go, then the outputs are compared.

The output of this tool is a diff-like-thing (it's not the easiest thing to read, and could be
vastly improved).
"""

import bisect
import difflib
import itertools
import operator
import re
import subprocess

TEST_PATTERN = re.compile(r'<!-- BEGIN TEST ([^ ]+?) -->\n(.*)\n<!-- END TEST \1 -->', re.DOTALL)


def filter_empty(lines):
    return [line for line in lines if line]


def run_test(test_name, lines, line_number, file_name):
    # Ensure the first line is a command we know.
    commands = {
        '$ magma': ['magma', '-b']
    }
    assert lines[0] in commands

    # Get all the input lines, i.e. those starting with '> '. The output lines
    # are all the others.
    input_lines = [line[2:] for line in lines if line.startswith('> ')] + ['quit;']
    expected_output_lines = [line for line in lines[1:] if not line.startswith('> ')]

    proc = subprocess.Popen(
        commands[lines[0]],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        stdout, stderr = proc.communicate(''.join(line + '\n' for line in input_lines).encode('utf-8'), timeout=1)
    except subprocess.TimeoutExpired:
        proc.kill()
        stdout, stderr = proc.communicate()

    # Magma doesn't output to stderr at all I don't think...
    output_lines = stdout.decode('utf-8').split('\n')
    delta = difflib.Differ().compare(
        filter_empty(expected_output_lines),
        filter_empty(output_lines),
    )
    if all(line.startswith('  ') for line in delta):
        print(f"Test '{test_name}' passed")
    else:
        print(f"Test '{test_name}' failed: {file_name}:{line_number}")
        print(f"Showing diff from expected to actual output:")
        print('\n'.join(delta))


def run_tests_in_file(fname):
    with open(fname) as f:
        content = f.read()

    # Cumulative sum of character lengths per line, so that we can figure out what line a
    # character is on.
    linechars = list(itertools.accumulate((len(line) + 1 for line in content.split('\n')), operator.add))

    for match in TEST_PATTERN.finditer(content):
        name, content = match.groups()
        code_lines = [
            line[4:] for line in content.split('\n')
            if line.startswith(' '*4)
        ]
        line_number = 1 + bisect.bisect_right(linechars, match.start())
        run_test(name, code_lines, line_number, fname)


if __name__ == '__main__':
    run_tests_in_file('README.md')