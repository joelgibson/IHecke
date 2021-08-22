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
"""

import re
import subprocess
import difflib

TEST_PATTERN = re.compile(r'<!-- BEGIN TEST ([^ ]+?) -->\n(.*)\n<!-- END TEST \1 -->', re.DOTALL)


def filter_empty(lines):
    return [line for line in lines if line]


def run_test(test_name, lines):
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
        print(f"{test_name} passed")
    else:
        print(f"{test_name} failed: showing diff from expected to actual output:")
        print('\n'.join(delta))


def run_tests_in_file(fname):
    with open(fname) as f:
        for match in TEST_PATTERN.finditer(f.read()):
            name, content = match.groups()
            code_lines = [
                line[4:] for line in content.split('\n')
                if line.startswith(' '*4)
            ]
            run_test(name, code_lines)


if __name__ == '__main__':
    run_tests_in_file('README.md')