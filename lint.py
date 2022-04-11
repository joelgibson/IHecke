"""
Some code quality checks.
"""

import glob
import re


SOURCE_FILES = [
    *glob.glob("Package/*.m"),
    *glob.glob("MuCoeffs/*.m"),
    *glob.glob("PCanBases/*.m"),
    *glob.glob("Tests/*.m"),
    *glob.glob("Examples/*.m"),
]

def check_file(fname: str):
    lines = list(open(fname))
    for i, line in enumerate(lines):
        line_no = i + 1
        line = line.rstrip('\n')

        if '\t' in line:
            print(f"{fname}:{line_no} tab in line")

        if line.rstrip() != line:
            print(f"{fname}:{line_no} trailing space")

        if "New(BasisIHke)" in line:
            print(f"{fname}:{line_no} new BasisIHke created (abstract type should not be created)")

        if 'todo' in line.lower():
            print(f"{fname}:{line_no} TODO found: {line}")

        # Look for CartanName() being called, not preceded by a try statement.
        if re.search(r'(?<!_IHke)CartanName\(', line) and not (i > 0 and 'try' in lines[i-1]):
            print(f"{fname}:{line_no} calls CartanName(), it should call _IHkeCartanName()")


for fname in SOURCE_FILES:
    check_file(fname)


def search_file(pattern, fname):
    with open(fname) as f:
        match = re.search(pattern, f.read(), re.IGNORECASE)

    return match if match is None else match.group(1)


# Check that version numbers line up
code_date = search_file(r'return.*version.*(\d\d\d\d-\d\d-\d\d)', 'Package/FModIHke.m')
readme_date = search_file(r'(\d\d\d\d-\d\d-\d\d).*current', 'README.md')
if code_date is None:
    print("No version date found in the code")
elif readme_date is None:
    print("No version date found in the readme")
elif code_date != readme_date:
    print(f"Code date {code_date} and readme date {readme_date} disagree.")
