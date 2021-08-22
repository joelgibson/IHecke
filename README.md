# IHecke

`IHecke` is a package for [Magma] implementing Iwahori-Hecke algebras. It implements the standard
and canonical bases of the Hecke algebra (Using the Soergel convention, so a reader of the book
[Introduction to Soergel Bimodules][SBim] should feel at home). In addition, it comes bundled with
some p-canonical bases for small rank groups (all groups of rank at most 4, excluding F4), and can
calculate left, right, and two-sided cells and p-cells.

The term "Hecke algebra" is already used within Magma for a different kind of algebra, hence all of
the names in this package are prefixed with the letter I (standing for "Iwahori"), for example the
intrinsic`IHeckeAlgebra()` creates an element of type `AlgIHke`. Whenever we refer to the "Hecke
algebra", we really mean the Iwahori-Hecke algebra.

Two main goals this package has are to (A) present a neat, easy-to-use interface to the user, and
(B) have well-organised and easily-extendable internals for a package implementer. The user-facing
design should be clear enough that it could be used by someone with minimal Magma knowledge, and
for example could be used as an aid while reading some chapters of *Introduction to Soergel
Bimodules*. The internal design should make it straightforward to add new bases, without having to
change any existing code (see the implementation of the p-canonical basis in `AlhIHkePCan.m` for an
example of this). Having outstanding performance is *not* a goal at this time, so don't expect to be
able to use this package to generate Kazhdan-Lusztig polynomials for giant groups. The internal
design is discussed more at [Internals.md](Internals.md).

[Magma]: http://magma.maths.usyd.edu.au/magma/
[SBim]: https://www.springer.com/gp/book/9783030488253


- [Tutorial](#tutorial)
  * [Installation](#installation)
  * [Working with Coxeter groups in Magma](#working-with-coxeter-groups-in-magma)
  * [The standard basis](#the-standard-basis)
  * [The canonical basis](#the-canonical-basis)
  * [Operations on bases and elements](#operations-on-bases-and-elements)
  * [Cells](#cells)
  * [The p-canonical basis](#the-p-canonical-basis)
- [Examples](#examples)
  * [Printing basis elements](#printing-basis-elements)
  * [Displaying cells and p-cells](#displaying-cells-and-p-cells)
- [Changelog](#changelog)
- [TODO](#todo)


# Tutorial

## Installation

Download and unpack the `IHecke` directory somewhere. (You can get a zip of the current code by
going to the Code button in the top-right and clicking "Download zip", otherwise you can check the
code out using Git). If you're using it alongside your existing Magma files, your directory
structure should look something like this:

    ExistingMagmaFile.m
    IHecke/
        Examples/
        Package/
        PCanBases/
        Tests/
        IHecke.spec
        run-tests.sh
        ...

Inside of `ExistingMagmaFile.m`, the library can be loaded using `AttachSpec("IHecke/IHecke.spec")`.

The tests and examples are designed to be run from within the `IHecke` directory. Change into the
`IHecke` directory, and run

    $ bash run-tests.sh

You should see some output ending in:

    All tests successful.


## Working with Coxeter groups in Magma

<!-- BEGIN TEST CoxeterGroups -->

Before starting with `IHecke`, some familiarity with how to wrangle Coxeter groups inside of Magma
is helpful. Start up Magma and define a Coxeter group of your favourite type:

    $ magma
    > W := CoxeterGroup(GrpFPCox, "B3");

The name `GrpFPCox` is instructing Magma to build a Coxeter group backed by the finitely-presented
(`FP`) group machinery, which we need for some operations on words. The numbering convention being
used by Magma can be checked by calling `CoxeterDiagram` on the group:

    > CoxeterDiagram(W);

    B3    1 - 2 === 3

The identity of the group is `W.0`, and the three coxeter generators are `W.1`, `W.2`, and `W.3`
respectively. We can check some of the Coxeter relations:

    > W.1 * W.2;
    $.1 * $.2
    > (W.1 * W.2)^2;
    $.2 * $.1
    > (W.1 * W.2)^3;
    Id($)
    > (W.2 * W.3)^4 eq W.0;
    true

Another way of constructing elements of `W` is by coercing from a list of integers, which has the
same effect as multiplying generators:

    > W.1 * W.2 * W.3 eq W ! [1, 2, 3];
    true

In order to turn a group element back into a list of integers, call `Eltseq` on the element:

    > LongestElement(W);
    $.1 * $.2 * $.1 * $.3 * $.2 * $.1 * $.3 * $.2 * $.3
    > Eltseq(LongestElement(W));
    [ 1, 2, 1, 3, 2, 1, 3, 2, 3 ]

In small examples, it can be helpful to assign names like `s, t, u` to the generators. This can be
done by using `AssignNames`:

    > AssignNames(~W, ["s", "t", "u"]);
    > LongestElement(W);
    s * t * s * u * t * s * u * t * u

It can also be done at group construction time. The statement

    > W<s,t,u> := CoxeterGroup(GrpFPCox, "B3");

has the effect of setting `W` to the new Coxeter group, assigning the names `s`, `t`, and `u` to the
generators, and also binding `s`, `t`, and `u` to variables representing the generators, so we can
now do

    > s*t*s eq t*s*t;
    true

Binding `s`, `t`, and `u` could have also been done by calling `s := W.1`, `t := W.2`, and
`u := W.3`.

The length of an element is given by the `Length` or `#` function:

    > #(s*t*s);
    3
    > #(s*t*s*t);
    2

The ordering operators `lt` and `le` (and `gt` and `ge`) compare words lexicographically. To compare
in the Bruhat partial ordering instead, use `BruhatLessOrEqual(x, y)`.

More information can be found in the [Magma handbook](http://magma.maths.usyd.edu.au/magma/handbook/),
split between the "Finitely-Presented Groups" and "Lie Theory -> Coxeter Groups" sections.


<!-- END TEST CoxeterGroups -->


## The standard basis

<!-- BEGIN TEST Bases -->

To get going, attach the spec file. Depending on where `magma` is being run from (inside or outside
of the `IHecke` directory), the path to the spec file may change. Verify things are correctly loaded
by checking the version:

    $ magma
    > AttachSpec("IHecke/IHecke.spec");
    > IHeckeVersion();
    IHecke version 2021-08-18

Now, create a `GrpFPCox` of your favourite type:

    > W := CoxeterGroup(GrpFPCox, "B4");

Calling `IHeckeAlgebra(W)` creates a Hecke algebra object. This object is fairly useless on its own,
and just acts as a kind of "basis factory". Call `IHeckeAlgebraStd` on the Hecke algebra to create
the standard basis.

    > HAlg := IHeckeAlgebra(W);
    > HAlg;
    Iwahori-Hecke algebra of type B4
    > H := IHeckeAlgebraStd(HAlg);
    > H;
    Standard basis of Hecke algebra for Coxeter group of type B4, symbol H

A basis object like `H` models the Hecke algebra *together with a choice of basis*. Elements of the
basis can be created in two ways. The first is by using the `.` operator, which manufactures basis
elements. Passing a Coxeter group element to `H.` creates that basis element:

    > H.(W.1 * W.2 * W.1);
    (1)H(121)

There are a few shortcuts: when `i` is an integer, `H.i` has the same effect as `H.(W.i)`, so for
example `H.0` is the basis element corresponding to the identity, and `H.1`, `H.2`, and so on the
basis elements corresponding to the Coxeter generators. Passing a list of integers like `H.[1,2,1]`
first creates the Coxeter group element `W.1 * W.2 * W.1`, then creates the basis element.

    > H.0;
    (1)H(id)
    > H.[1,2,1] eq H.(W.1 * W.2 * W.1);
    true
    > H.LongestElement(W) + H.1 * H.1;
    (1)H(1213214321432434) + (v^-1 - v)H(1) + (1)H(id)

Another way of creating elements is via *coercion*, which means the `!` operator in Magma. With an
integer or Laurent polynomial on the right, `H!` acts as the unit map of the algebra, taking a
scalar and returning the identity times that scalar.

    > H ! 2;
    (2)H(id)
    > LPoly<v> := BaseRing(H);
    > H ! (v^10 - v^-1);
    (-v^-1 + v^10)H(id)
    > H.0 * (v^10 - v^-1);
    (-v^-1 + v^10)H(id)

The main use of the coercion operator `!` will be changing between bases, which we will see below.



## The canonical basis

The canonical basis is created using `IHeckeAlgebraCan` on the parent Hecke algebra object:

    > C := IHeckeAlgebraCan(HAlg);
    > C;
    Canonical basis of Hecke algebra for Coxeter group of type B4, symbol C

Elements of the canonical basis can be created using the same `.` and `!` operators as in the
standard basis:

    > C.0 - C.1 + C.LongestElement(W) + C.1;
    (1)C(1213214321432434) + (1)C(id)

In the calculation above, nothing interesting has happened: internally every basis is treated just
as a free module over `W` (so addition, subtraction, and scalar multiplication internal to `H` is
programmed in exactly the same way as for `C`, and any other basis). The interesting stuff starts
happening when we use the coercion `!` operator to change bases:

    > H ! C.1;
    (1)H(1) + (v)H(id)
    > H ! C.[1, 2, 1];
    (1)H(121) + (v)H(21) + (v)H(12) + (v^2)H(2) + (v^2)H(1) + (v^3)H(id)

The coefficients appearing in these expansions are examples of Kazhdan-Lusztig polynomials (or
rather the Soergel normalisation of them). They can be extracted using the `Coefficient` function:

    > C_121 := C.[1, 2, 1];
    > Coefficient(C_121, W.1);
    0
    > Coefficient(H ! C_121, W.1);
    v^2

Notice that the coefficient of the Coxeter group element `W.1` in the canonical basis element
`C_121`was `0` when in the canonical basis, but `v^2` when expressed in the standard basis.

In mixed-basis arithmetic, whichever basis appears on the left is used as the output basis. In
mixed-basis comparisons, the right element is (internally) converted to the left basis, then
elements are compared as formal linear combinations.

    > H.1 + C.1;
    (2)H(1) + (v)H(id)
    > C.1 + H.1;
    (2)C(1) + (-v)C(id)
    > H.1 + C.1 eq C.1 + H.1;
    true



## Operations on bases and elements

Before moving on, we'll list the various operations that can be done on basis objects and linear
combinations.

The basis of an element can be accessed using `Parent`, and the Hecke algebra of a basis can be
accessed using `Parent` again:

    > Parent(H.1) eq H;
    true
    > Parent(H) eq HAlg;
    true

The functions `CoxeterGroup` and `BaseRing` can be called on a Hecke algebra or basis. Furthermore,
`BasisSymbol` and `BasisName` can be called on a basis.

    > CoxeterGroup(H) eq W;
    true
    > BaseRing(H);
    Laurent series field in v over Integer Ring
    > BasisName(H);
    Standard basis
    > BasisSymbol(H);
    H

The support and coefficients of a linear combination can be accessed using `Support`. The first
return value is an indexed set (set with a defined ordering, denoted by `{@ ... @}` in Magma) of
Coxeter group elements, and the second return value is a list giving the coefficients of the support
in the same order as the set.

    > H.[2,1] * H.1;
    (v^-1 - v)H(21) + (1)H(2)
    > supp, coeffs := Support(H.[2,1] * H.1);
    > supp;
    {@ $.2, $.2 * $.1 @}
    > coeffs;
    [
        1,
        v^-1 - v
    ]

In Magma, if a function returning multiple values is used in a single-value context, then the
function is treated as returning only the first argument. Hence we can check the size of the support
by using the cardinality `#` operator:

    > #Support(H.[2,1] * H.1);
    2

In order to extract a particular coefficient, use the `Coefficient` function.

    > Coefficient(H.[2,1] * H.1, W![2,1]);
    v^-1 - v
    > Coefficient(H.[2,1] * H.1, W.1);
    0

The Kazhdan-Lusztig bar involution can be calculated using `Bar`.

    > Bar(H.1);
    (1)H(1) + (-v^-1 + v)H(id)
    > Bar(v * C.[1,2,1]);
    (v^-1)C(121)


<!-- END TEST Bases -->


## Cells

<!-- BEGIN TEST Cells -->

A based algebra defines a partitioning of its indexing set into equivalence classes called *cells*,
and furthermore a partial order on those partitions called the *cell order*. There are three
different flavours of cells: *left* cells, *right* cells, and *two-sided* cells. All three of these
can be computed at once by running the function `Cells` on a basis object. We demonstrate this in
an S3 example:

    $ magma
    > SetColumns(0);
    > AttachSpec("IHecke/IHecke.spec");
    > W := CoxeterGroup(GrpFPCox, "A2");
    > HAlg := IHeckeAlgebra(W);
    > C := IHeckeAlgebraCan(HAlg);
    > left, right, twoSided := Cells(C);
    > left, right, twoSided;
    4 Left cells of Canonical basis of Hecke algebra for Coxeter group of type A2, symbol C
    4 Right cells of Canonical basis of Hecke algebra for Coxeter group of type A2, symbol C
    3 Two-sided cells of Canonical basis of Hecke algebra for Coxeter group of type A2, symbol C

Each of `left`, `right`, and `twoSided` is a bundle of information about the computed cells. The
partition can be extracted by using `Components`:

    > Components(left);
    {@
    { $.1 * $.2 * $.1 },
    { $.1, $.2 * $.1 },
    { $.2, $.1 * $.2 },
    { Id($) }
    @}


The particular cell that an element falls into can be extracted using `Component`:

    > Component(left, W.1);
    { $.1, $.2 * $.1 }

Two Coxeter group elements can be compared in the cell ordering by using `CellLe`:

    > CellLe(left, W.0, W.1);
    true
    > CellLe(left, W.1, W.0);
    false
    > CellLe(left, W.1, W.1);
    true

<!-- END TEST Cells -->

## The p-canonical basis

`IHecke` comes bundled with some pre-calculated p-canonical bases. The algorithm which actually
computes these basis elements is quite complicated and is not a part of `IHecke`, so there is no
general method of constructing a p-canonical basis element: if it's not pre-calculated you're out
of luck. Running `IHeckeAlgebraPCan` with a Cartan type and a prime will do one of three things:

1. Print `p-canonical basis ... was loaded from the database`, and return a p-canonical basis
    object, loaded from one of the pre-calculated datasets.
2. Print `p-canonical basis ... coincides with canonical basis`, and return a canonical basis object.
3. Print `p-canonical basis ... unavailable`, and throw an error.

The rules for telling when the p-canonical basis coincides with the canonical basis are just as
ad-hoc as the datasets, and more or less boil down to a collection of theoretical results from some
papers, as well as empirical results (they were calculated directly, and checked to be equal).

A p-canonical basis can be constructed in much the same way as the standard or canonical basis, but
it requires two extra arguments: a cartan type and a prime. (The cartan type cannot be inferred from
the Hecke algebra in general, since the p-canonical basis depends on the root system rather than the
Coxeter group. In short, B_n and C_n have different canonical bases).

<!-- BEGIN TEST p-canonical -->

    $ magma
    > AttachSpec("IHecke/IHecke.spec");
    > SetColumns(0);
    > W := CoxeterGroup(GrpFPCox, "B2");
    > HAlg := IHeckeAlgebra(W);
    > H := IHeckeAlgebraStd(HAlg);
    > C := IHeckeAlgebraCan(HAlg);
    > pC := IHeckeAlgebraPCan(HAlg, "B2", 2);
    The 2-canonical basis for type B2 was loaded from the database.

As with the other bases, all of the interest is in converting between bases:

    > pC.[2,1,2];
    (1)p2C(212)
    > C ! pC.[2,1,2];
    (1)C(212) + (1)C(2)

In fact, 212 is the only Coxeter element for which the 2-canonical basis of B2 differs from the
canonical basis. We can confirm this:

    > [C ! pC.w : w in EnumerateCoxeterGroup(W) | C.w ne pC.w];
    [
    (1)C(212) + (1)C(2)
    ]

<!-- END TEST p-canonical -->


# Examples

There are some short programs in the `Examples/` directory, showing how the library can be used. The
examples are designed to be run from inside the same directory as `IHecke.spec`.


## Printing basis elements

The `Examples/PrintCanonicalBasis.m` program shows how to use `IHecke` as a calculator for canonical
basis elements in terms of the standard basis, by printing out either a list or a table. To print a
list:

    $ magma -b type:=A2 Examples/PrintCanonicalBasis.m
    Numbering convention for type A2:

    A2    1 - 2

    Canonical basis:
    (1)C(id) = (1)H(id)
    (1)C(1) = (1)H(1) + (v)H(id)
    (1)C(2) = (1)H(2) + (v)H(id)
    (1)C(12) = (1)H(12) + (v)H(2) + (v)H(1) + (v^2)H(id)
    (1)C(21) = (1)H(21) + (v)H(2) + (v)H(1) + (v^2)H(id)
    (1)C(121) = (1)H(121) + (v)H(21) + (v)H(12) + (v^2)H(2) + (v^2)H(1) + (v^3)H(id)

To print a table:

    $ magma -b type:=A2 tabular:=true Examples/PrintCanonicalBasis.m
        id  1   2   12  21  121
    id  1
    1   v   1
    2   v       1
    12  v^2 v   v   1
    21  v^2 v   v       1
    121 v^3 v^2 v^2 v   v   1


## Displaying cells and p-cells

The `Examples/GraphLeftCells.m` program shows how to use `IHecke` in conjunction with
[Graphviz](https://graphviz.org/) to graphically display the left cells or p-cells for a group. (The
left cells are the left cells of the canonical basis, the left p-cells are the left cells of the
p-canonical basis). You will need to have Graphviz installed (on Linux or Mac, it should be
available from your package manager). You can check it is installed by running

    $ dot -V
    dot - graphviz version 2.43.0 (0)

The program `Examples/GraphLeftCells.m` produces a text file describing a directed graph, which the
Graphviz program `dot` then recieves and turns into an image file. To see what is happening, try
running the example without running it through `dot`:

    $ magma -b type:=G2 prime:=3 Examples/GraphLeftCells.m

You should see a fairly understandable description of a directed graph. Note that we used the `-b`
flag above: this tells Magma not to print its usual preamble and postamble, which would confuse
GraphViz at the next step. To create the image, run the output through the `dot` layout program: we
also first run it through `tred` which computes the transitive reduction of the graph, whittling the
edges down to a minimal subset which generate the cell order:

    $ magma -b type:=G2 prime:=3 Examples/GraphLeftCells.m | tred | dot -Tpng -o g2p3.png

After opening `g2p3.png`, you should see a graph much like one of the ones below:

<p align="center">
<img alt="G2, characteristic zero cells" src="Examples/GraphLeftCells-G2p0.gif?raw=true">
&nbsp;
<img alt="G2, characteristic 2 cells" src="Examples/GraphLeftCells-G2p2.gif?raw=true">
&nbsp;
<img alt="G2, characteristic 3 cells" src="Examples/GraphLeftCells-G2p3.gif?raw=true">
</p>

The cell indices (`LCell #2` for example) are not meaningful.


# Changelog

- Development version
  - Added example of how to display left cells using Graphviz.
- Version 2021-08-18 (Current)
  - Initial release, supporting bases `Std`, `Can`, and `PCan` of the Hecke algebra, as well as cell
      computations.


# TODO

- (High priority) Better handling of the cell graph
  - Perform transitive reduction, so that it is in a canonical form.
  - Make it available to the user (interface design), and an example involving GraphVis.
- (High priority) Spherical and Antispherical left/right modules.
- (Cleanup) Document the use of `eq` on `AlgIHkeBase`.
- (Low priority) Allow more custom formatting of the output, such as choosing basis names, and
    choosing formatting of Coxeter group elements.
- (Low priority) Performance:
  - Write another implementation of cells using only the mu-coefficients from the canonical basis,
      which should be faser than the naive version.
  - Try swapping out the Laurent series ring for a rational function field, and see if operations
    go faster (especially checking self-duality using the bar involution, and building cells).
    - When I tried this (the only operation which is unclear is taking the coefficient of v, but
        this can be done by coercing into the Laurent series ring and doing it there), things
        performed about on-par, a little slower if anything.
    - With the current Laurent series implementation, taking a profile of some code which double
        checks self-duality in B4, almost all the time (11 seconds) was spent in `_AddScaled`, and
        of that 11 seconds about 3 were spent in `IsDefined`, `*`, and `+`. This means that most of
        the time is simply spent interpreting the function `_AddScaled`.