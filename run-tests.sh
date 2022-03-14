#!/usr/bin/env bash

# Run each test file
for f in Tests/Test*.m; do
    echo -n "Running $f ... "

    # For some reason passing batch:=true here (or any other argument) causes magma to segfault.
    magma -b "$f"
    if [ $? -ne 0 ]; then
        exit
    fi
done

# Run each example file at least once, to make sure they all execute properly.
declare -a examples=(
    "magma -b type:=A3 Examples/PrintCanonicalBasis.m"
    "magma -b type:=A3 tabular:=true Examples/PrintCanonicalBasis.m"
    "magma -b type:=G2 Examples/GraphLeftCells.m"
    "magma -b type:=G2 prime:=2 Examples/GraphLeftCells.m"
    "magma -b type:=G2 prime:=3 Examples/GraphLeftCells.m"
    "magma -b type:=A3 Examples/BenchmarkCanInStd.m"
    "magma -b type:=A3 notable:=true Examples/BenchmarkCanCanMultiplication.m"
    "magma -b type:=A5 Tests/LongTestMusCorrect.m"
)
for cmd in "${examples[@]}"; do
    if ! $cmd > /dev/null; then
        echo "Example failed, rerun:" $cmd
        exit
    else
        echo "Example successful:" $cmd
    fi
done

echo "All tests successful."
