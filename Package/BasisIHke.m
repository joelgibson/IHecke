import "EltIHke.m":
    _AddScaled,
    _RemoveZeros,
    _AddScaledTerm;


declare type BasisIHke;
declare attributes BasisIHke:
    FreeModule,     // A free module type, i.e. an object of a type inheriting from FModIHke.
    BasisSymbol,    // A short string naming the basis, used for printing eg H(id) + H(121).
    BasisName;      // A human-readable name describing the basis.


////////////////////
// Initialisation

intrinsic _BasisIHkeInit(~basis::BasisIHke, fmod::FModIHke, symbol::MonStgElt, name::MonStgElt)
{Bases (types extending BasisIHke) should call _BasisIHkeInit to initialise common fields.}
    require Type(basis) ne BasisIHke: "BasisIHke is an abstract type, and should not be created.";
    basis`FreeModule := fmod;
    basis`BasisSymbol := symbol;
    basis`BasisName := name;
end intrinsic;


/////////////////////
// Accessor functions

intrinsic Print(A::BasisIHke)
{}
    printf "%o of %o, symbol %o",
        BasisName(A),
        Name(FreeModule(A)),
        BasisSymbol(A);
end intrinsic;

intrinsic FreeModule(A::BasisIHke) -> IHkeAlg
{The free module for which this is a basis (the Hecke algebra, antispherical module, etc).}
    return A`FreeModule;
end intrinsic;

intrinsic CoxeterGroup(A::BasisIHke) -> GrpFPCox
{The underlying Coxeter group of the Hecke algebra.}
    return CoxeterGroup(A`FreeModule);
end intrinsic;

intrinsic BasisSymbol(A::BasisIHke) -> MonStgElt
{A short string naming the basis.}
    return A`BasisSymbol;
end intrinsic;

intrinsic BasisName(A::BasisIHke) -> MonStgElt
{A human-readable name describing the basis.}
    return A`BasisName;
end intrinsic;

intrinsic BaseRing(A::BasisIHke) -> Rng
{Return the base ring (a Laurent series ring).}
    return BaseRing(A`FreeModule);
end intrinsic;


intrinsic 'eq'(A::BasisIHke, B::BasisIHke) -> BoolElt
{By default, two bases compare equal if they have the same type, and their free modules
 compare equal. Bases parameterised on some additional data (eg p-canonical bases) should override
 this intrinsic.}
    return Type(A) eq Type(B) and FreeModule(A) cmpeq FreeModule(B);
end intrinsic;

intrinsic _IHkeProtToBasis(A::BasisIHke, B::BasisIHke, w::GrpFPCoxElt) -> EltIHke
{Fallback function for expressing B(w) in the A basis. Bases must override this at least twice, to
 convert to and from the default basis.}
    return false;
end intrinsic;

// TODO: Re-evaluate whether this was a good choice.
intrinsic _IHkeProtToBasisElt(A::BasisIHke, B::BasisIHke, elt::EltIHke) -> EltIHke
{Default implementation, lifting the map w -> elt to a map elt -> elt.}
    W := CoxeterGroup(A);
    if Type(_IHkeProtToBasis(A, B, W.0)) ne EltIHke then
        return false;
    end if;


    terms := AssociativeArray(W);
    for w -> coeff in elt`Terms do
        _AddScaled(~terms, _IHkeProtToBasis(A, B, w)`Terms, coeff);
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(A, terms);
end intrinsic;


intrinsic _IHkeProtUnit(A::BasisIHke) -> EltIHke
{If A is the basis of an algebra, this should return the unit element.}
    return false;
end intrinsic;



//////////////////////////////
// Creation of basis elements

intrinsic '.'(A::BasisIHke, w::GrpFPCoxElt) -> EltIHke
{The basis element indexed by the Coxeter group element w.}
    require Parent(w) eq CoxeterGroup(A):
        "Group element", w, "is not a member of", CoxeterGroup(A);

    terms := AssociativeArray(CoxeterGroup(A));
    terms[w] := BaseRing(A) ! 1;
    return EltIHkeConstruct(A, terms);
end intrinsic;

intrinsic '.'(A::BasisIHke, s::RngIntElt) -> EltIHke
{The basis element indexed by the Coxeter group element W.s, equivalent to A.(W.s).}
    w := CoxeterGroup(A) . s;
    terms := AssociativeArray(CoxeterGroup(A));
    terms[w] := BaseRing(A) ! 1;
    return EltIHkeConstruct(A, terms);
end intrinsic;

intrinsic '.'(A::BasisIHke, word::SeqEnum[RngIntElt]) -> EltIHke
{The basis element indexed by the Coxeter group word, equivalent to A.(W!word).}
    w := CoxeterGroup(A) ! word;
    terms := AssociativeArray(CoxeterGroup(A));
    terms[w] := BaseRing(A) ! 1;
    return EltIHkeConstruct(A, terms);
end intrinsic;


///////////
// Coercion

intrinsic IsCoercible(A::BasisIHke, elt::RngElt) -> BoolElt, EltIHke
{The unit map of the algebra: coerce a scalar to the identity multiplied by that scalar.}
    // Check that r coerces into the base ring.
    ok, r := IsCoercible(BaseRing(A), elt);
    if not ok then
        return false, Sprintf("%o does not coerce into the base ring %o", elt, BaseRing(A));
    end if;

    // If r is zero, return the zero element (this is valid in any free module, while other scalars
    // are only valid if the module also provides a unit map).
    if r eq 0 then
        return true, EltIHkeConstruct(A, AssociativeArray(CoxeterGroup(A)));
    end if;

    // Check if this specific basis defines a unit.
    unit := _IHkeProtUnit(A);
    if Type(unit) eq EltIHke then
        return true, unit * r;
    end if;

    // Check if the default basis defines a unit.
    def := DefaultBasis(FreeModule(A));
    unit := _IHkeProtUnit(def);
    if Type(unit) eq EltIHke then
        return true, BasisChange(A, def, unit * r);
    end if;

    // We assume that there is no algebra structure.
    return false, Sprintf("No unit map provided for", A, "to use to insert the scalar", r);
end intrinsic;

intrinsic IsCoercible(A::BasisIHke, x::EltIHke) -> BoolElt, EltIHke
{Change the basis of x into A.}
    if FreeModule(A) ne FreeModule(x) then
        return false, Sprintf("Cannot coerce from %o into %o (different free modules)", FreeModule(x), FreeModule(A));
    end if;

    // Easy case: if x is already in the A basis, return x.
    if A eq Parent(x) then
        return true, x;
    end if;

    // Otherwise, perform a basis conversion.
    return true, ChangeBasis(A, Parent(x), x);
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

intrinsic ChangeBasis(A::BasisIHke, B::BasisIHke, elt::EltIHke) -> EltIHke
{Express elt, which must be in the B basis, in the A basis.}
    assert FreeModule(A) eq FreeModule(B);
    assert B eq Parent(elt);
    W := CoxeterGroup(A);

    // Simple case: same basis.
    if A eq B then
        return elt;
    end if;

    // Test to see if there is a direct basis conversion A <- B.
    result := _IHkeProtToBasisElt(A, B, elt);
    if Type(result) eq EltIHke then
        return result;
    end if;

    // Otherwise, convert via the default basis.
    D := DefaultBasis(FreeModule(A));
    if Type(_IHkeProtToBasis(A, D, W.0)) eq EltIHke and Type(_IHkeProtToBasis(D, B, W.0)) eq EltIHke then
        return ChangeBasis(A, D, ChangeBasis(D, B, elt));
    end if;

    error Sprintf("No basis conversion found to (%o) from (%o), i.e. %o from %o", A, B, Type(A), Type(B));
end intrinsic;
