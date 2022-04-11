# Tables of μ-coefficients

This directory contains tables of μ-coefficients (i.e. W-graphs) for some finite groups of medium rank.
These coefficients can be used to implement fast multiplications CxH -> C and CxC -> C without needing to calculate Kazhdan-Lusztig polynomials first.
This gives a great speed boost to some use-cases, like calculating products C(s1) ... C(sn), or calculating cell decompositions of bespoke bases.
The tables are calculated using a different program written in Go and fashioned after Coxeter3 by Fokku du Cloux.

Each file contains three pieces of data:

1. The Coxeter matrix used to calculate the table, used to double-check we are about to load the right table.
2. A partial multiplication table, used to communicate a bijection `N: W -> [1, |W|]` labelling group elements by integers. The partial multiplication table is a sparse integer matrix `M` of dimensions `|W| x |S|`, where the entry `M[N(w), s] = N(sw)` if `w < sw` and `s` is the first left descent of `sw`, and `M[N(w), s] = 0` otherwise. The bijection will always satisfy `N(id) = 1`, but is otherwise arbitrary.
3. A sparse integer matrix of dimensions `|W| x |W|`, where the row with index `N(w)` contains the μ-coefficients `μ(x, w)` for all `x < w` with `l(x) <= l(w) - 3`.

Each file can be loaded into MAGMA quickly by a simple `eval`.

The label `N(w)` can be calculated on any element `w` recursively, by finding the first left descent `s`, calculating `N(sw)`, and then applying the partial multiplication `M[N(sw), s]`.
The base case for the recursion is `N(id) = 1`.
In practice (at least in MAGMA) this recursive approach is too slow: we just build the bijection as the file is loaded.

Whenever `x` is covered by `w` in the Bruhat order, meaning that `x <= w` and `l(x) = l(w) - 1`, we have `μ(x, w) = 1`.
Since the set of elements covered by `w` is easy to compute, we have chosen to not store these coefficients in the table, instead the row for `w` stores only those `x` such that `l(x) <= l(w) - 3`.
(Recall that if the length difference between `x` and `w` is even, then `μ(x, w) = 0`).
The reader of the table is responsible for completing these missing mu-coefficients.


## Why make tables of μ-coefficients?

The only algorithm I know for computing μ-coefficients is the straightforward one: calculate all Kazhdan-Lusztig polynomials, and take the coefficients of `v` in each.
Calculating all Kazhdan-Lusztig polynomials fast requires using a language other than MAGMA, which runs too slowly for this purpose (it suffers from the same speed drawbacks as other interpreted languages like Python).
Until I find a better algorithm (perhaps inducing W-graphs?) which can be run in MAGMA fast enough, loading tables is fine.


## Other technical details

The tables are loaded in as sparse integer matrices, partly because MAGMA has a special fast mode for reading in integer sequences.
We only provide tables for the cases where IHecke takes over a second to compute mu-coefficients for the whole group.
By running the benchmarks

    magma -b type:=B4 notable:=true Examples/BenchmarkCanInStd.m

we can find the frontier.

- A5 takes 0.45 seconds but A6 takes 58 seconds. We include tables for A5, A6, and A7 (A5 is included so that we have something to run `Tests/LongTestMusCorrect` with, that doesn't take too long).
- B4 takes 0.25 seconds but B5 takes 50 seconds. We include tables for B5 and B6.
- D4 takes 0.05 seconds but D5 takes 6 seconds. We include tables for D5 and D6.
- E6 is too large to time using Magma. We include a table for E6.
- F4 takes 3 seconds, so we include a table for F4.
- H3 takes 0.03 seconds but H4 takes
