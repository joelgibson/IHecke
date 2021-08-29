import "EltIHke.m": _AddScaled, _AddScaledTerm, _RemoveZeros;

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

intrinsic Parent(A::BasisIHke) -> IHkeAlg
{}
    require false: "Call FreeModule() on a basis instead of Parent()";
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

intrinsic _ToBasis(A::BasisIHke, B::BasisIHke, w::GrpFPCoxElt) -> EltIHke
{Fallback implementation. Bases should override either this, or the EltIHke version.}
    return false;
end intrinsic;

intrinsic _ToBasis(A::BasisIHke, B::BasisIHke, eltB::EltIHke) -> EltIHke
{Fallback implementation, which uses _ToBasis(::BasisIHke, ::BasisIHKe, ::GrpFPCoxElt) if available.
 Bases should override either this, or the GrpFPCoxElt version.}
    W := CoxeterGroup(A);
    terms := AssociativeArray(W);
    for w -> coeff in eltB`Terms do
        summand := _ToBasis(A, B, w);
        if Type(summand) ne EltIHke then
            return false;
        end if;
        _AddScaled(~terms, summand`Terms, coeff);
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(A, terms);
end intrinsic;

intrinsic ToBasis(A::BasisIHke, B::BasisIHke, eltB::EltIHke) -> EltIHke
{Change eltB, which must be in the B basis, into the A basis.}
    require FreeModule(A) eq FreeModule(B): "Free modules", FreeModule(A), "and", FreeModule(B), "incompatible.";
    require B eq Parent(eltB): "The basis", B, "disagrees with the element basis", Parent(eltB);

    // Simple case: if the bases are equal, nothing needs to be done.
    if A eq B then return eltB; end if;

    // Is there a direct basis conversion A <- B defined?
    eltA := _ToBasis(A, B, eltB);
    if Type(eltA) eq EltIHke then
        error if A ne Parent(eltA), "Basis conversion output", Parent(eltA), "instead of", A;
        return eltA;
    end if;

    // Convert via the standard basis of the free module.
    StdBasis := StandardBasis(FreeModule(A));
    eltStd := _ToBasis(StdBasis, B, eltB);
    error if Type(eltStd) ne EltIHke, "No basis conversion defined for", B, "into", StdBasis;
    error if Parent(eltStd) ne StdBasis, "Basis conversion output", Parent(eltStd), "instead of", StdBasis;

    eltA := _ToBasis(A, StdBasis, eltStd);
    error if Type(eltA) ne EltIHke, "No basis conversion defined for", StdBasis, "into", A;
    error if Parent(eltA) ne A, "Basis conversion output", Parent(eltA), "instead of", A;

    return eltA;
end intrinsic;

intrinsic _Unit(A::BasisIHke) -> EltIHke
{If A is the basis of an algebra, this should return the unit element.}
    return false;
end intrinsic;


//////////////////////////////
// Creation of basis elements

intrinsic '.'(A::BasisIHke, w::GrpFPCoxElt) -> EltIHke
{The basis element indexed by the Coxeter group element w.}
    require Parent(w) eq CoxeterGroup(A): "Group element", w, "is not a member of", CoxeterGroup(A);

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
    unit := _Unit(A);
    if Type(unit) eq EltIHke then
        return true, unit * r;
    end if;

    // Check if the standard basis defines a unit.
    def := StandardBasis(FreeModule(A));
    unit := _Unit(def);
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
    return true, ToBasis(A, Parent(x), x);
end intrinsic;