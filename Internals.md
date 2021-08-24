# IHecke Internals

The internals of `IHecke` are designed to be extensible in an "open" way, so that new bases and
definitions can be added without having to change existing code. However, be aware that the
internals of this library are in some sense a medium-sized experiment with pushing the limits of
inheritance and multiple dispatch. A vague hope is that some of this patterns can be ported to
another, faster language with multiple dispatch, such as Julia.


# Types

There are three "levels" of objects defined in `IHecke`.

1. At the top level we have "free module" objects, for example the Hecke algebra itself, or the
    antispherical and spherical modules. These objects themselves do not implement any arithmetic,
    but serve as a container for some data (the Coxeter group, a parabolic subset, etc) as well as
    a cache object for generating bases. Each free module type must extend `FModIHke`.
2. At the middle level we have "basis" objects such as `IHkeAlgStd`, `IHkeAlgCan`, and so on. These
    basis objects model a free module *together with a choice of basis*. Each basis type must extend
    `BasisIHke`.
3. At the bottom level we have element objects of type `EltIHke`, which is the data of a basis
    object (its `Parent`), and a map of Coxeter group elements to scalars (the linear combination).

The types `FModIHke` and `BasisIHke` are "abstract" types, in the sense that there should be types
inheriting from these types, but to ever actually create an element of this type (eg by calling
`New(FModIHke)` or `New(BasisIHke)`) is an error. On the flip side, `EltIHke` is a concrete type
which should never need to be extended at all.

An example of how the types fit together:

    > W := CoxeterGroup(GrpFPCox, "B4");
    > HAlg := IHeckeAlgebra(W);             // HAlg has type IHkeAlg, extending FModIHke
    > H := StandardBasis(HAlg);             // H has type IHkeAlgStd, extending BasisIHke
    > C := CanonicalBasis(HAlg);            // C has type IHkeAlgCan, extending BasisIHke
    > elt := H.1;                           // elt has type EltIHke, with Parent(elt) eq H.


# Background: overloading and multiple dispatch

The extensibility of `IHecke` is implemented by overloading intrinsics in Magma. Dispatch of
intrinsics in Magma is similar in spirit to a language like Julia, or the CLOS (Common Lisp Object
System). It is helpful to know how this works, i.e. how Magma selects the implementation of an
intrinsic to run, given the types of its arguments.

Suppose we have the following type definitions, defining types A and B, with B inheriting from A:

    declare type A;
    declare type B: A;

Now we have an intrinsic which accepts an object of type `A` and an integer, and returns an integer:

    intrinsic Foo(a::A, n::RngIntElt) -> RngIntElt
    {}
        return n;
    end intrinsic;

Supposing that `a` and `b` are objects of types `A` and `B` respectively, both function calls
`Foo(a, 3)` and `Foo(b, 3)` will return `3`. The reason that `Foo(b, 3)` works is that `B` inherits
from `A`, which can be checked using `ISA(B, A)` in Magma.

We can add a more specific overload of `Foo` by defining `Foo` again with different argument types:

    intrinsic Foo(b::B, n::RngIntElt) -> RngIntElt
    {}
        return n^2;
    end intrinsic;

This time, running `Foo(a, 3)` will return `3`, but running `Foo(b, 3)` will return `9`. There was
nothing special here about the parameter being in the first argument either: Magma will dispatch on
the type of all arguments of an intrinsic (intrinsics can be defined to take up to 6 arguments).

`IHecke` takes advantage of multiple dispatch by defining internal functions which not only take an
`EltIHke` as an argument, but the relevant basis object(s) too, allowing the programmer to easily
specify operations on certain bases. The fact that these same-named intrinsics automatically fit
together and overload eachother based on specificity enables the "open" aspect, where all that needs
to be done to add a new basis is to define a new type and overload a new intrinsic (rather than, for
instance, edit a giant case-by-case statement inside of `IHecke` itself).


# Operations on free modules

When defining a free module, declare a new type inheriting from `FModIHke`, and call the intrinsic
`_FModIHkeInit` during construction. In addition, the free module should define an implementation
of `eq`, and an implementation of `StandardBasis`, which should return a standard / default basis
of the module.


# Operations on bases

When defining a basis, declare a new type inheriting from `BasisIHke`, and call `_BasisIHkeInit`
during construction. A new basis *must* define a basis conversion to and from the standard basis
of its free module, after which all operations implemented on the standard basis will be
automatically implemented for the new basis. More specific versions of these operations can be
implemented *à la carte*, which one might want to do for efficiency (some operation is much faster
when implemented internally, like the bar operation on the canonical basis), or correctness checking
(comparing two different basis-specific implementations of the same abstract map).


## Changing bases

`IHecke` provides the user-facing function `ToBasis` for changing bases:

    intrinsic ToBasis(A::BasisIHke, B::BasisIHke, eltB::EltIHke) -> EltIHke
    {Express eltB, which must be in the B basis, in the A basis.}

This function should not be overloaded directly, but it calls out to the intrinsic `_ToBasis` to do
the work, which should be overloaded for each pair of bases the programmer wants to supply a
conversion for. When adding a new basis, there *must* be at least a conversion to and from the
standard basis, since the user-facing function `ToBasis` will fall back to converting via the
standard basis if it does not find a direct conversion from B to A.

In order to hook into the `ToBasis` function, the author should define exactly one of the intrinsics

    intrinsic _ToBasis(A::BasisIHke, B::BasisIHke, w::GrpFPCoxElt) -> EltIHke
    intrinsic _ToBasis(A::BasisIHke, B::BasisIHke, eltB::EltIHke) -> EltIHke

replacing the two abstract types `BasisIHke` mentioned with specific bases. For example, to define
a change of basis from the canonical basis into the standard basis, the an intrinsic with the
following signature should be defined:

    intrinsic _ToBasis(H::IHkeAlgStd, C::IHkeAlgCan, w::GrpFPCoxElt) -> EltIHke

After defining `_ToBasis`, basis conversions during coercions like `H ! C.[1,2,1]`, as well as basis
conversions used to implement operations inside `IHecke`, will work automatically.


## Multiplication and action maps

`IHecke` overloads `'*'(::EltIHke, ::EltIHke)`, which calls out to the intrinsic

    intrinsic _Multiply(A::BasisIHke, eltA::EltIHke, B::BasisIHke, eltB::EltIHke) -> EltIHke

When `_Multiply` is called, we will have `A eq Parent(eltA)` and `B eq Parent(eltB)`. The result
should be in the `A` basis.

In order to implement multiplication in a new free module, or the action of one free module on
another, the programmer should overload `_Multiply` for the standard bases of that module. After
doing this, `_Multiply` may be overridden for other pairs of bases as required.


## Unit map

`IHecke` will automatically coerce `0` into any free module, but will only coerce a nonzero scalar
if a "unit" map (as in the unit map of an algebra) is provided. One must be provided for the
standard basis, and other bases may provide their own if desired.

    intrinsic _Unit(A::BasisIHke) -> EltIHke

Once this map is provided, nonzero scalars can be coerced into the basis, and empty products of
basis elements will yield the unit object.


## Bar involution

`IHecke` provides the function `Bar` to perform the bar-involution on modules:

    intrinsic Bar(elt::EltIHke) -> EltIHke

Internally, this intrinsic calls

    intrinsic _Bar(A::BasisIHke, eltA::EltIHke) -> EltIHke

which should be overloaded first for the standard basis (if `Bar` should be implemented on a new
free module), and then may be overloaded basis-by-basis as desired.


# Internal assumptions

## Testing

It is important for testing sanity purposes that arithmetic, comparisons, etc do not change basis if
both elements involved are already in the same basis. For example, `C.1 + C.2 eq C.2 + C.1` should
never invoke a basis change.

## Unitriangularity

We assume that every basis is unitriangular to the standard basis. (Wait, do we actually? I'm not
sure that much of the code relies on this).

- The assumption that `B.0` is the identity for every basis is used in the coercion-from-scalar code.
- It would be best if `B!1` always manufactured the identity element, for consistency.


# Notes

- I tried a prototype that did not have a top level free module object, but making some of the
    user-facing operations nice was hard. For example, to multiply in an arbitrary basis, the
    fallback function will convert to the standard basis, multiply there, and convert back. For this
    conversion to happen, a standard basis object needs to be "available" to the expression `x * y`,
    so we call `StandardBasis(FreeModule(x))`, which returns the standard basis object which
    probably already exists.