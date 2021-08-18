# p-canonical bases

These p-canonical bases were calculated offline, using an algorithm due largely to Thorge Jensen and
Geordie Williamson.

Since it's difficult to load files on a path in Magma (I really want to be able to refer to a path
relative to the current source file, but I seem to only be able to get paths relative to where the
user is executing Magma from), rather than loading these files directly, running the script

    python3 generate.py > PCanDB.m

will consolidate all of the `Type-*-P-*` files into a Magma source file, which can then be used
without needing to worry about paths. This database file just maps file names to lists of lines.