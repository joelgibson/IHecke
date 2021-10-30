# IHecke

`IHecke` is a package for [Magma] implementing Iwahori-Hecke algebras. It implements the standard
and canonical bases of the Hecke algebra (Using the Soergel convention, so a reader of the book
[Introduction to Soergel Bimodules][SBim] should feel at home). In addition, it comes bundled with
some p-canonical bases for small rank groups (all groups of rank at most 4, excluding F4), and can
calculate left, right, and two-sided cells and p-cells. It also contains implementations of the
antispherical and spherical modules (and their canonical bases), for any parabolic subgroup.

Two main goals this package has are to (A) present a neat, easy-to-use interface to the user, and
(B) have well-organised and easily-extendable internals for a package implementer. The user-facing
design should be clear enough that it could be used by someone with minimal Magma knowledge, and
for example could be used as an aid while reading some chapters of *Introduction to Soergel
Bimodules*. The internal design should make it straightforward to add new bases, without having to
change any existing code (see [Defining a new basis](#defining-a-new-basis) below for an example).
Having outstanding performance is *not* a goal at this time, so don't expect to be
able to use this package to generate Kazhdan-Lusztig polynomials for giant groups. The internal
design is discussed more at [Internals.md](Internals.md).

The term "Hecke algebra" is already used within Magma for a different kind of algebra, hence all of
the names in this package are prefixed with the letter I (standing for "Iwahori"), for example the
intrinsic`IHeckeAlgebra()` creates an element of type `IHkeAlg`. Whenever we refer to the "Hecke
algebra", we really mean the Iwahori-Hecke algebra.

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
  * [Antispherical and spherical modules](#antispherical-and-spherical-modules)
- [Examples](#examples)
  * [Printing basis elements](#printing-basis-elements)
  * [Displaying cells and p-cells](#displaying-cells-and-p-cells)
  * [Defining a new basis](#defining-a-new-basis)
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
The rest of this tutorial, along with all the examples and tests, are designed to be run from within
the `IHecke/` directory. After changing into the `IHecke/` directory, run the tests using

    $ bash run-tests.sh

You should see some output ending in

    All tests successful.

These tests double-check that the library is working correctly with your version of Magma. The
library does use some newer features of Magma, such as dual iteration (`for k -> v in assocs`) which
was introduced in version 2.25 (January 2020).


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
    > AttachSpec("IHecke.spec");
    > IHeckeVersion();
    IHecke version 2021-09-10

Now, create a `GrpFPCox` of your favourite type:

    > W := CoxeterGroup(GrpFPCox, "B4");

Calling `IHeckeAlgebra(W)` creates a Hecke algebra object. This object is fairly useless on its own,
and just acts as a kind of "basis factory". Call `StandardBasis` on the Hecke algebra to create
the standard basis.

    > HAlg := IHeckeAlgebra(W);
    > HAlg;
    Iwahori-Hecke algebra of type B4
    > H := StandardBasis(HAlg);
    > H;
    Standard basis of Iwahori-Hecke algebra of type B4, symbol H

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

Note that `H.[1, 1]` is the same as the identity, since `H.[1, 1] = H.(W.1 * W.1) = H.(W.0)`. In
order to turn a list like `[1, 1]` into a product of *algebra* generators `H.1 * H.1`, use the
[reduction operator](http://magma.maths.usyd.edu.au/magma/handbook/text/117#998) `&*` together with
a list comprehension:

    > H.[1, 1];
    (1)H(id)
    > &*[H.i : i in [1, 1]];
    (v^-1 - v)H(1) + (1)H(id)

Empty lists can be used with the reduction operator `&*`, so long as Magma is informed of what
universe the elements of the list are from:

    > &*[H | ]; // Empty list with universe H, elements could be added with [H | H.1, H.2] etc.
    (1)H(id)

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

The canonical basis is created using `CanonicalBasis` on the parent Hecke algebra object:

    > C := CanonicalBasis(HAlg);
    > C;
    Canonical basis of Iwahori-Hecke algebra of type B4, symbol C

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



## Common operations

Before moving on, we'll list the various operations that can be done on basis objects and linear
combinations.

Firstly, there are some "shortcut" methods for creating the Hecke algebra, group, and bases all at
once. `ShortcutIHeckeAlgebra` can be passed a `GrpFPCox`, and will return the Hecke algebra,
standard, and canonical bases.

    > W := CoxeterGroup(GrpFPCox, "B4");
    > HAlg, H, C := ShortcutIHeckeAlgebra(W);

The shortcut methods `ShortcutIHeckeAntiSpherical(HAlg)` and `ShortcutIHeckeSpherical(HAlg)` are
similar, returning a triple of (module, standard basis, canonical basis) for the antispherical
and spherical modules.

The basis of an element can be accessed using `Parent`, and the Hecke algebra of a basis can be
accessed using the `FreeModule` intrinsic on either the basis or an element:

    > Parent(H.1) eq H;
    true
    > FreeModule(H) eq HAlg;
    true
    > FreeModule(H.1) eq HAlg;
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
    > AttachSpec("IHecke.spec");
    > W := CoxeterGroup(GrpFPCox, "A2");
    > HAlg := IHeckeAlgebra(W);
    > C := CanonicalBasis(HAlg);
    > left, right, twoSided := Cells(C);
    > left, right, twoSided;
    4 Left cells of Canonical basis of Iwahori-Hecke algebra of type A2, symbol C
    4 Right cells of Canonical basis of Iwahori-Hecke algebra of type A2, symbol C
    3 Two-sided cells of Canonical basis of Iwahori-Hecke algebra of type A2, symbol C

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

A generating set of the edges between components can be listed using `GeneratingEdges`, which can
then be used to draw a Hasse diagram of the cells (see the
[example below](#displaying-cells-and-p-cells)). Each edge is represented as a sequence containing
two components, with the cell component on the left being lower than the one on the right.

    > GeneratingEdges(left);
    [
    [
    { $.1, $.2 * $.1 },
    { $.1 * $.2 * $.1 }
    ],
    [
    { $.2, $.1 * $.2 },
    { $.1 * $.2 * $.1 }
    ],
    [
    { Id($) },
    { $.1, $.2 * $.1 }
    ],
    [
    { Id($) },
    { $.2, $.1 * $.2 }
    ]
    ]


The `GeneratingEdges` function returns a generating set but not necessarily a minimal one. In a
future version of `IHecke`, it will only return minimal generating sets for the cell order.


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
    > AttachSpec("IHecke.spec");
    > SetColumns(0);
    > W := CoxeterGroup(GrpFPCox, "B2");
    > HAlg := IHeckeAlgebra(W);
    > H := StandardBasis(HAlg);
    > C := CanonicalBasis(HAlg);
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


## Antispherical and spherical modules

The antispherical and spherical modules are two modules which are naturally associated to the Hecke
algebra and a standard parabolic subgroup. Currently in `IHecke` they are implemented as right
modules only.

In order to create the antispherical or spherical module, set up a Hecke algebra and bases as usual:

<!-- BEGIN TEST Antispherical -->

    $ magma
    > AttachSpec("IHecke.spec"); SetColumns(0);
    > W := CoxeterGroup(GrpFPCox, "B4");
    > HAlg := IHeckeAlgebra(W);
    > H := StandardBasis(HAlg);
    > C := CanonicalBasis(HAlg);

In the numbering convention of Magma, the type `A3` standard parabolic sitting inside of `B4` is
generated by the first three generators `I = [1, 2, 3]`. To form the right antispherical module
for this standard parabolic, use `IHeckeAntiSpherical` on the Coxeter group object:

    > ASMod := IHeckeASMod(W, [1, 2, 3]);
    > ASMod;
    Antispherical module of type B4, parabolic [ 1, 2, 3 ]

(*Note:* it is important that the exact same Coxeter group object is used to create the Hecke
algebra and the antispherical/spherical modules, since Magma will not automatically coerce between
distinct Coxeter groups with identical presentations).

The generating set can be retrieved again using `Parabolic`:

    > Parabolic(ASMod);
    [ 1, 2, 3 ]

The standard basis is created by calling `StandardBasis` on the `ASMod` module object:

    > aH := StandardBasis(ASMod);
    > aH;
    Standard basis of Antispherical module of type B4, parabolic [ 1, 2, 3 ], symbol aH

Bases of the antispherical and spherical modules are naturally indexed by the right cosets of the
standard parabolic subgroup `W_I` inside the Coxeter group `W`. Rather than using right cosets, we
use instead the unique minimal-length representative inside each coset, i.e. that unique element
whose left descent set does not intersect `I = [1, 2, 3]`. This can be checked using `IsMinimal`:

    > IsMinimal([1, 2, 3], W.1);
    false
    > IsMinimal([1, 2, 3], W.4 * W.3);
    true

Creating basis elements works in the same way for the antispherical and spherical modules as it does
for the Hecke algebra, except that a few operations which were previously legal are no longer legal.
The standard and canonical bases will only accept Coxeter group elements which are `I`-minimal, and
throw an error otherwise:

    > aH.1;
    Runtime error: $.1 is not minimal with respect to [ 1, 2, 3 ]
    > aH.[4, 3];
    (1)aH(43)

Coercion from scalars only works for the zero scalar.

    > aH ! 0;
    (0)aH(id)
    > aH ! 1;

    >> aH ! 1;
          ^
    Runtime error in '!': Standard basis of Antispherical module of type B4, parabolic [ 1, 2, 3 ], symbol aH is not an algebra (no unit map), so cannot insert scalar 1
    LHS: ASModIHkeStd
    RHS: RngIntElt

Other than this, basis elements can be added, subtracted, and multiplied by scalars as expected.

The action of the Hecke algebra on the module can be performed using the `*` operator, keeping in
mind that these are *right* modules:

    > aH.0 * H.[1,2,1];
    (-v^3)aH(id)
    > aH.4 * H.[3,4];
    (1)aH(434)
    > aH.0 * C.[1,2,1];
    (0)aH(id)

Attempting to multiply elements of the module with itself will fail:

    > aH.0 * aH.0;
    Runtime error: Multiplication not defined between Antispherical module of type B4, parabolic [ 1, 2, 3 ]
    and Antispherical module of type B4, parabolic [ 1, 2, 3 ]

The canonical basis of the Antispherical module can be created using `CanonicalBasis`:

    > aC := CanonicalBasis(ASMod);
    > aC;
    Canonical basis of Antispherical module of type B4, parabolic [ 1, 2, 3 ], symbol aC

Alternatively, both bases (along with the module itself) can be created using a shortcut method:

    > ASMod, aH, aC := ShortcutIHeckeASMod(W, [1, 2, 3]);

Basis conversions between the canonical and standard bases work as expected.

    > aH ! aC.[4,3,4];
    (1)aH(434) + (v)aH(43)

Furthermore, the Hecke algebra can act on the canonical basis directly (internally, this converts to
the standard basis, performs the action, and converts back).

    > aC.[4,3] * C.4;
    (1)aC(434) + (1)aC(4)

The spherical module works analagously to the antispherical module. It can be created using
`IHeckeSpherical(W, [1, 2, 3])`, with its standard and canonical bases created using
`StandardBasis` and `CanonicalBasis` respectively. There is also a shortcut function:

    > SMod, sH, sC := ShortcutIHeckeSMod(W, [1, 2, 3]);
    > SMod;
    Spherical module of type B4, parabolic [ 1, 2, 3 ]
    > sH;
    Standard basis of Spherical module of type B4, parabolic [ 1, 2, 3 ], symbol sH
    > sC;
    Canonical basis of Spherical module of type B4, parabolic [ 1, 2, 3 ], symbol sC

<!-- END TEST Antispherical -->


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


## Defining a new basis

The design of `IHecke` is "open", meaning that new bases, modules, etc are able to be defined by
users without having to modify the package itself. This example shows how to add a new (simple)
basis of the Hecke algebra in about 20 lines of code.

The "new basis" we will define is the standard basis `T(w)` as it was defined in Kazhdan-Lusztig's
seminal work [KL79]. Its quadratic relation is `(T(s) + 1)(T(s) - v^-2) = 0` - note that while
[KL79] works over the indeterminate `q`, we will stick with the indeterminate `v`, using `v^-2 = q`.

Defining a new basis requires defining four things: a new basis 'type', a function creating a basis
object of that type, and then two basis conversion functions: one from the new basis to the standard
basis, and one from the standard basis to the new basis. Have a read of the well-commented file
[`Examples/KLStandard.m`](Examples/KLStandard.m) to see how this is done.

We will demonstrate the new basis in action. Since bases are defined using intrinsics, which can
only be loaded using `.spec` files, we will load
[`Examples/KLStandard.spec`](Examples/KLStandard.spec) right after loading `IHecke`.

<!-- BEGIN TEST KLStandard -->

    $ magma
    > AttachSpec("IHecke.spec"); SetColumns(0);
    > AttachSpec("Examples/KLStandard.spec");
    > W := CoxeterGroup(GrpFPCox, "B4");
    > HAlg, H, C := ShortcutIHeckeAlgebra(W);

The new basis is creating using the function defined in `Examples/KLStandard.m`:

    > T := StandardBasisKL(HAlg);
    > T;
    KL-standard basis of Iwahori-Hecke algebra of type B4, symbol T

Multiplication is automatically defined (behind-the-scenes, the algebra elements are changed into
the standard basis, multiplied there, and then changed back into the `T` basis). We can check the
quadratic relation:

    > T.1 * T.1;
    (v^-2 + -1)T(1) + (v^-2)T(id)

We can examine canonical basis elements (our `C` corresponds to `C'` in [KL79]) in the `T` basis:

    > T ! C.[2,1,3,2];
    (v^4)T(2132) + (v^4)T(232) + (v^4)T(213) + (v^4)T(132) + (v^4)T(121) + (v^4)T(32) + (v^4)T(23) + (v^4)T(21) + (v^4)T(13) + (v^4)T(12) + (v^4)T(3) + (v^2 + v^4)T(2) + (v^4)T(1) + (v^2 + v^4)T(id)

According to [KL79], `C'(w)` should be scaled by `v^(-#w)` in order to have the coefficients coming
out as the original KL polynomials `P(y, w)`:

    > v := BaseRing(T).1;
    > T ! (v^-4 * C.[2,1,3,2]);
    (1)T(2132) + (1)T(232) + (1)T(213) + (1)T(132) + (1)T(121) + (1)T(32) + (1)T(23) + (1)T(21) + (1)T(13) + (1)T(12) + (1)T(3) + (v^-2 + 1)T(2) + (1)T(1) + (v^-2 + 1)T(id)

For example, we have the nontrivial Kazhdan-Lusztig polynomial `P(2, 2132) = v^-2 + 1 = q + 1`.

<!-- END TEST KLStandard -->

[KL79]: https://eudml.org/doc/142660


# Changelog

Breaking changes are marked with a (!).

- Version 2021-09-10 (Current)
  - Added exponentiation operator `^`, for nonnegative powers.
- Version 2021-08-24
  - Many internals relabelled.
  - Added example of how to display left cells using Graphviz.
  - Added example of defining a new basis.
  - Added antispherical and spherical modules.
  - (!) Renamed `IHeckeAlgebraStd` to `StandardBasis`, and `IHeckeAlgebraCan` to `CanonicalBasis`.
  - (!) Removed the `Parent` property from bases, and replaced it with `FreeModule`. This property
      can be accessed on a basis or element by the `FreeModule(...)` intrinsic.
- Version 2021-08-18
  - Initial release, supporting bases `Std`, `Can`, and `PCan` of the Hecke algebra, as well as cell
      computations.


# TODO

- (High priority) Allow the multiplication `C(w) * C.s` to be driven by a table of mu-coefficients
    (in other words, a W-graph) so that this special case is extremely fast.
- (Low priority) Perform transitive reduction on the cell graph.
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