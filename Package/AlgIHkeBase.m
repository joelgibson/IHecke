import "EltIHke.m": _LaurentPolyRing, _EltIHkeConstruct, _AddScaled, _RemoveZeros, _AddScaledTerm;


declare type AlgIHkeBase;
declare attributes AlgIHkeBase:
    // AlgIHke
    Parent,

    // A short string naming the basis, used for printing eg H(id) + H(121).
    BasisSymbol,

    // A human-readable name describing the basis.
    BasisName;


////////////////////
// Initialisation

// An AlgIHkeBase should never be created, but rather types extending it should be created.
// This procedure initialises the common attributes.
procedure _AlgIHkeBaseInit(~basis, alg, basisSymbol, basisName)
    assert ISA(Type(alg), AlgIHke);
    assert ISA(Type(basis), AlgIHkeBase);
    assert Type(basisSymbol) eq MonStgElt;
    assert Type(basisName) eq MonStgElt;

    basis`Parent := alg;
    basis`BasisSymbol := basisSymbol;
    basis`BasisName := basisName;
end procedure;


/////////////////////
// Accessor functions

intrinsic Print(A::AlgIHkeBase)
{}
    printf "%o of Hecke algebra for Coxeter group of type %o, symbol %o",
        BasisName(A),
        CartanName(CoxeterGroup(A)),
        BasisSymbol(A);
end intrinsic;

intrinsic Parent(A::AlgIHkeBase) -> AlgIHke
{The parent Hecke algebra.}
    return A`Parent;
end intrinsic;

intrinsic CoxeterGroup(A::AlgIHkeBase) -> GrpFPCox
{The underlying Coxeter group of the Hecke algebra.}
    return A`Parent`CoxeterGroup;
end intrinsic;

intrinsic BasisSymbol(A::AlgIHkeBase) -> MonStgElt
{A short string naming the basis.}
    return A`BasisSymbol;
end intrinsic;

intrinsic BasisName(A::AlgIHkeBase) -> MonStgElt
{A human-readable name describing the basis.}
    return A`BasisName;
end intrinsic;

intrinsic BaseRing(A::AlgIHkeBase) -> Rng
{Return the base ring (a Laurent series ring).}
    return _LaurentPolyRing;
end intrinsic;


//////////////////////////////
// Creation of basis elements

intrinsic '.'(A::AlgIHkeBase, w::GrpFPCoxElt) -> EltIHke
{The basis element indexed by the Coxeter group element w.}
    require Parent(w) eq CoxeterGroup(A):
        "Group element", w, "is not a member of", CoxeterGroup(A);

    terms := AssociativeArray(CoxeterGroup(A));
    terms[w] := _LaurentPolyRing ! 1;
    return _EltIHkeConstruct(A, terms);
end intrinsic;

intrinsic '.'(A::AlgIHkeBase, s::RngIntElt) -> EltIHke
{The basis element indexed by the Coxeter group element W.s, equivalent to A.(W.s).}
    w := CoxeterGroup(A) . s;
    terms := AssociativeArray(CoxeterGroup(A));
    terms[w] := _LaurentPolyRing ! 1;
    return _EltIHkeConstruct(A, terms);
end intrinsic;

intrinsic '.'(A::AlgIHkeBase, word::SeqEnum[RngIntElt]) -> EltIHke
{The basis element indexed by the Coxeter group word, equivalent to A.(W!word).}
    w := CoxeterGroup(A) ! word;
    terms := AssociativeArray(CoxeterGroup(A));
    terms[w] := _LaurentPolyRing ! 1;
    return _EltIHkeConstruct(A, terms);
end intrinsic;


//////////////////////
// Coercion of scalars

intrinsic IsCoercible(A::AlgIHkeBase, r::RngElt) -> BoolElt, EltIHke
{The unit map of the algebra: coerce a scalar to the identity multiplied by that scalar.}
    ok, result := IsCoercible(_LaurentPolyRing, r);
    if ok then
        return true, A.0 * r;
    end if;
end intrinsic;

intrinsic IsCoercible(A::AlgIHkeBase, x::EltIHke) -> BoolElt, EltIHke
{Change the basis of x into A.}
    // Hecke algebras must be the same.
    if Parent(A) ne Parent(Parent(x)) then
        return false;
    end if;

    // Easy case: if x is already in the A basis, return x.
    if A eq Parent(x) then
        return true, x;
    end if;

    // Otherwise, perform a basis conversion.
    return true, ChangeBasis(A, Parent(x), x);
end intrinsic;


/////////////////////
// Predicates

intrinsic 'eq'(A::AlgIHkeBase, B::AlgIHkeBase) -> BoolElt
{Two Hecke algebra bases are 'equal' if they refer to the same named basis. They may or may not be
 actually equal as a set of W-indexed elements in the Hecke algebra. This function should be
 overridden.}
    return false;
end intrinsic;


//////////////////////////
// Basis change (protocol)
//
// Basis change is done by coercion, eg to put an element of the canonical basis into the standard
// basis, use something like (H ! C.w).
//
// Each basis type can specialise the intrinsic _IHkeProtToBasis for any bases it cares about, but it
// must implement at least a conversion to the standard or canonical basis, and a conversion from
// the standard or canonical basis.
//
// A basis type does not need to specialise the intrinsic for conversions from itself to itself,
// this is handled below.

intrinsic ChangeBasis(A::AlgIHkeBase, B::AlgIHkeBase, elt::EltIHke) -> EltIHke
{Express elt, which must be in the B basis, in the A basis.}
    assert Parent(A) eq Parent(B);
    assert B eq Parent(elt);
    W := CoxeterGroup(A);

    // Test to see if there is a direct basis conversion A <- B.
    if Type(_IHkeProtToBasis(A, B, W.0)) eq EltIHke then
        terms := AssociativeArray(W);
        for w -> coeff in elt`Terms do
            _AddScaled(~terms, _IHkeProtToBasis(A, B, w)`Terms, coeff);
        end for;
        _RemoveZeros(~terms);
        return _EltIHkeConstruct(A, terms);
    end if;

    // Otherwise, see if there is an indirect conversion A <- X <- B, where X is either the
    // standard basis or the canonical basis.
    H := IHeckeAlgebraStd(Parent(A));
    if Type(_IHkeProtToBasis(A, H, W.0)) eq EltIHke and Type(_IHkeProtToBasis(H, B, W.0)) eq EltIHke then
        return ChangeBasis(A, H, ChangeBasis(H, B, elt));
    end if;

    C := IHeckeAlgebraCan(Parent(A));
    if Type(_IHkeProtToBasis(A, C, W.0)) eq EltIHke and Type(_IHkeProtToBasis(C, B, W.0)) eq EltIHke then
        return ChangeBasis(A, C, ChangeBasis(C, B, elt));
    end if;

    error Sprintf("No basis conversion found to (%o) from (%o)", A, B);
end intrinsic;


intrinsic _IHkeProtToBasis(A::AlgIHkeBase, B::AlgIHkeBase, w::GrpFPCoxElt) -> EltIHke
{Fallback function for expressing B(w) in the A basis.}
    if A eq B then
        return B.w;
    end if;

    return false;
end intrinsic;
