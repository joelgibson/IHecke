#!/usr/bin/env bash

# Run each test file
for f in Tests/*; do
    # if magma -b batch:=true $f; then
    #     true
    # else
    #     echo $f $'failed. Run\n    magma $f\nto debug.';
    # fi;
    echo -n "Running $f ... "
    if ! magma -b batch:=true $f; then
        exit
    fi
done

# Run a select set of examples, hiding output
declare -a examples=(
    "magma -b type:=A3 Examples/PrintCanonicalBasis.m"
    "magma -b type:=A3 tabular:=true Examples/PrintCanonicalBasis.m"
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