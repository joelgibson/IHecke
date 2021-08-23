# IHecke Internals

A lot of thought has gone into designing the internals of this library in an "open" way, so that new
bases and definitions can be added without having to change any existing code. I explain below how
this is done (and for seeing an example of success, check the implementation of the p-canonical
basis in `Package/IHkeAlgPCan.m`). My hope is that I can take some of these patterns over to a
faster language (I have Julia in mind) and implement them there.


## Types

There are three "levels" of objects defined in `IHecke`.

1. At the top level, an `IHkeAlg` object bundles the data of a Coxeter group together with a cache
    for constructing bases, so that calling `IHeckeAlgebraCan(HAlg)` on the same `HAlg` object twice
    will yield the same basis object (which is relevant, since basis objects uses caches to speed up
    their operation).
2. At the middle level we have basis objects such as `IHkeAlgStd`, `IHkeAlgCan`, and so on. These
    types all extend a common base type `BasisIHke`, which allows us to neatly define common
    operations and fallback methods. There is never an actual object created of type `BasisIHke`.
3. At the bottom level we have element objects of type `EltIHke`, which is the data of a basis
    object (its parent), and a map of Coxeter group elements to scalars (the linear combination).


Defining the types in this way (especially that middle level) allows us to take full advantage of
Magma's multiple dispatch intrinsics, so that we can define a "default" or "fallback" implementation
of an operation using the base type `BasisIHke`, but specific bases can define their own versions
of these functions, which get automatically used in preference to the default or fallback. I'm
calling the extensibility built out of this pattern a *protocol*.


Some notes:

- I tried a prototype that did not have a top level `IHkeAlg` object, but making some of the
    user-facing operations nice was hard. For example, to multiply in an arbitrary basis, the
    fallback function will convert to the standard basis, multiply there, and convert back. For this
    conversion to happen, a standard basis object needs to be "available" to the expression `x * y`,
    so we call `IHeckeAlgebraStd(Parent(Parent(x))`, which returns the standard basis object which
    probably already exists.


## Protocols

The pattern I am calling "protocols" is relying on the multiple dispatch aspect of Magma intrinsics
in order to have a clean and efficient way of extending `IHecke` with new bases. More or less this
boils down to calling some specially named functions and including the basis as an explicit argument
so that Magma will look up first if any function is implemented for that specific basis before
falling back to some kind of default implementation.

Two exceptions are the standard and canonical bases: these participate in the protocols just like
any other basis, but some of the top-level functions will fall back to them for implementations (eg
multiplication falls back to multiplication in the standard basis). Therefore the library considers
these bases to be "special".


### Example (Bar involution)

Consider the case of the bar involution, a unary operation on the elements of the Hecke algebra. Its
type signature should be:

    intrinsic Bar(elt::EltIHke) -> EltIHke

Furthermore we require that it should perform the bar involution on `elt`, and return the result in
the same basis as `elt`. However, that type signature does not mention the basis, so Mamga cannot
choose a different implementation of `Bar` depending on the basis of `elt`. To get around this, we
define a second intrinsic:

    intrinsic _IHkeProtBar(B::BasisIHke, elt::EltIHke) -> EltIHke
    {Fallback implementation of the bar involution.}
        return false;
    end intrinsic;

The underscore at the start of `_IHkeProtBar` is a convention which indicates an internal function
which users should not call directly, and the `IHkeProt` part indicates that the function is part
of a protocol. We can now define the body of the user-facing `Bar` function:

    intrinsic Bar(elt::EltIHke) -> EltIHke
    {Perform the bar involution on elt, returning the result in the same basis as elt.}
        result := _IHkeProtBar(Parent(elt), elt);

        if Type(result) eq EltIHke then
            return result;
        end if;

        // In the actual library, we convert to the canonical basis, perform bar, and convert back,
        // rather than throwing an error.
        error "No Bar involution implemented for", Parent(elt);
    end intrinsic;

With only these two intrinsics defined, any call to `Bar` reaches the fallback implementation of
`_IHkeProtBar`, and so an error is thrown.

When defining a new basis, the library author only needs to write a version of the `_IHkeProtBar`
for their basis:

    intrinsic _IHkeProtBar(H::IHkeAlgStd, elt::EltIHke) -> EltIHke
    {Perform the bar involution on an element of the standard basis.}

and now any time `Bar` is called on an element of the standard basis, the more specific intrinsic
takes precedence over the fallback intrinsic, and we get an answer.

I take the "Bar protocol" to mean the collection of intrinsics `Bar` and `_IHkeProtBar`, together
with the contract about how these things will be invoked: for example, whenever `_IHkeProtBar` is
invoked on `(A, elt)`, we have `A eq Parent(elt)`.


### Basis conversion

The non-coercive function for changing bases is

    intrinsic ChangeBasis(A::BasisIHke, B::BasisIHke, elt::EltIHke) -> EltIHke
    {Express elt, which must be in the B basis, in the A basis.}

The contract for calling this function is that we must have `B eq Parent(elt)` (an error will be
thrown otherwise). In order to hook into the `ChangeBasis` function, the author should define the
intrinsic

    intrinsic _IHkeProtToBasis(A::BasisIHke, B::BasisIHke, w::GrpFPCoxElt) -> EltIHke

replacing the two abstract types `BasisIHke` mentioned with specific bases. For example, to define
a change of basis from the canonical basis into the standard basis, the an intrinsic with the
following signature should be defined:

    intrinsic _IHkeProtToBasis(H::IHkeAlgStd, C::IHkeAlgCan, w::GrpFPCoxElt) -> EltIHke

When defining a new basis, the author *must* define at least a conversion from the new basis into
either the standard or canonical basis, and a conversion from either of the standard or canonical
basis into the new basis. The `ChangeBasis` function will automatically fall back to going via these
bases if a direct conversion from `B` into `A` is not defined by `_IHkeProtToBasis`.


### Multiplication

The function for multiplying two elements is `'*'`, which calls the protocol function

    intrinsic _IHkeProtMult(A::BasisIHke, eltA::EltIHke, B::BasisIHke, eltB::EltIHke) -> EltIHke

When `_IHkeProtMult` is called, we will have `A eq Parent(eltA)` and `B eq Parent(eltB)`. The result
should be in the `A` basis.

The implementation of `'*'` falls back to a conversion into the standard basis, multiplying there,
and converting back into the basis of `eltA`.


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