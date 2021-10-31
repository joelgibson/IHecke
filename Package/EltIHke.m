// An EltIHke a formal linear combination of Coxeter group elements, along with a Parent object
// which determines a basis. For example, (H, 1*w) would be w in the standard basis, while (C, 1*w)
// would be w in the canonical basis.
declare type EltIHke;
declare attributes EltIHke:
    Parent, // A basis, i.e. a type inheriting from BasisIHke.
    Terms;  // An AssociativeArray mapping Coxeter group elements to ring elements. This array is
            // normalised so that there are no elements mapping to zero.


////////////////////////////
// Construction

intrinsic _EltIHkeValidate(B::BasisIHke, elt::EltIHke)
{Bases should override this intrinsic and throw an error if elt contains any illegal terms. For
 example, bases for the antispherical and spherical modules should reject non-minimal elements.}
end intrinsic;

intrinsic EltIHkeConstruct(B::BasisIHke, terms::Assoc) -> EltIHke
{Constructor for EltIHke (only this function should ever be used).}
    // Check that zeros are normalised out.
    assert forall{w : w -> coeff in terms | coeff ne 0};
    assert Universe(terms) eq CoxeterGroup(B);

    elt := New(EltIHke);
    elt`Parent := B;
    elt`Terms := terms;
    _EltIHkeValidate(B, elt);
    return elt;
end intrinsic;


//////////////////////////
// Accessors

function FmtElt(w)
    return #w eq 0
        select "id"
        else   &cat[IntegerToString(i) : i in Eltseq(w)];
end function;

intrinsic Print(elt::EltIHke)
{Print a Hecke algebra or module element in a regular form. Suppose the basis is "H". Then,
a scalar multiple of H(w) is printed as (scalar)H(w), with the identity element of the Coxeter
group element being rendered as 'id'. Zero is printed as (0)H(id).}
    symbol := BasisSymbol(elt`Parent);
    keys := Sort(Setseq(Keys(elt`Terms)));
    if IsZero(elt) then
        printf "(%o)%o(%o)", 0, symbol, FmtElt(CoxeterGroup(elt`Parent).0);
    else
        printf Join([
                Sprintf("(%o)%o(%o)", elt`Terms[key], symbol, FmtElt(key))
                : key in Reverse(Sort(Setseq(Keys(elt`Terms))))
            ], " + ");
    end if;
end intrinsic;

intrinsic Parent(elt::EltIHke) -> BasisIHke
{The parent structure (a type inheriting from BasisIHke or ModIHke).}
    return elt`Parent;
end intrinsic;

intrinsic FreeModule(elt::EltIHke) -> FModIHke
{The free module to which this element belongs.}
    return FreeModule(Parent(elt));
end intrinsic;

intrinsic Coefficient(elt::EltIHke, w::GrpFPCoxElt) -> RngElt
{Return the coefficient of w in the linear combination.}
    if IsDefined(elt`Terms, w) then
        return elt`Terms[w];
    else
        return BaseRing(Parent(elt)) ! 0;
    end if;
end intrinsic;

intrinsic Support(elt::EltIHke) -> SetIndx[GrpFPCoxElt], SeqEnum[RngElt]
{Return two parallel sequences giving the group elements and coefficients in the linear combination.
 The support (the group elements) are returned as an indexed set, which can be converted to a
 sequence by using Setseq(...), the ordering is the inbuilt ordering on GrpFPCox, a lexicographic
 ordering refining the Bruhat ordering.}
    keys := {@ x : x in Sort(Setseq(Keys(elt`Terms))) @};
    return keys, [elt`Terms[k] : k in keys];
end intrinsic;


/////////////////////
// Predicates

intrinsic IsZero(elt::EltIHke) -> BoolElt
{Returns true if the element is zero.}
    return #elt`Terms eq 0;
end intrinsic;

intrinsic 'eq'(elt1::EltIHke, elt2::EltIHke) -> BoolElt
{Compares whether two elements are equal. This has the same semantics as IsZero(elt1 - elt2).}
    // Change into the left basis if necessary.
    if Parent(elt1) ne Parent(elt2) then
        elt2 := ToBasis(Parent(elt1), Parent(elt2), elt2);
    end if;

    // Associative arrays have no 'eq' function. Instead we check manually: two functions defined on
    // finite domains are equal if:
    //   1: The domains have equal cardinality,
    //   2: Every element of the first domain is in the second domain, and
    //   3: The values agree at all of these points.
    a := elt1`Terms;
    b := elt2`Terms;
    if #a ne #b then
        return false;
    end if;

    for w -> coeff in a do
        if not IsDefined(b, w) or b[w] ne coeff then
            return false;
        end if;
    end for;

    return true;
end intrinsic;


////////////////////////////////////
// AssociativeArray helpers

// Remove all k -> v pairs where v is zero.
procedure _RemoveZeros(~terms)
    keys := [k : k -> v in terms | IsZero(v)];
    for k in keys do
        Remove(~terms, k);
    end for;
end procedure;

// termsA <- termsA + scalar * termsB
procedure _AddScaled(~termsA, termsB, scalar)
    for w -> v in termsB do
        if IsDefined(termsA, w) then
            termsA[w] +:= v * scalar;
        else
            termsA[w] := v * scalar;
        end if;
    end for;
end procedure;

// termsA <- termsA + Basis(w) * scalar
procedure _AddScaledTerm(~terms, w, scalar)
    // Early exit is an optimisation, unnecessary for correctness. We don't want to be creating a
    // bunch of (w -> 0) entries in the associative array only to delete them later.
    if scalar eq 0 then
        return;
    end if;

    if IsDefined(terms, w) then
        terms[w] +:= scalar;
    else
        terms[w] := scalar;
    end if;
end procedure;

// Check whether the terms of the associative array are "unitriangular to w", i.e. w appears with
// coefficient 1, and all other nonzero terms are lower in the Bruhat order.
function _IsUniTriangular(terms, w)
    return IsDefined(terms, w)
        and terms[w] eq 1
        and forall{u : u -> _ in terms | BruhatLessOrEqual(u, w)};
end function;


////////////////////////
// R-module operations

intrinsic '-'(elt::EltIHke) -> EltIHke
{The negation of an element.}
    terms := AssociativeArray(CoxeterGroup(Parent(elt)));
    for w -> coeff in elt`Terms do
        terms[w] := -coeff;
    end for;
    return EltIHkeConstruct(Parent(elt), terms);
end intrinsic;

intrinsic '+'(elt1::EltIHke, elt2::EltIHke) -> EltIHke
{The sum of two elements. Will throw an error if adding over different free modules.}
    error if FreeModule(elt1) ne FreeModule(elt2),
        "Cannot add in different free modules", FreeModule(elt1), "and", FreeModule(elt2);

    // Change into the left basis if necessary.
    if Parent(elt1) ne Parent(elt2) then
        elt2 := ToBasis(Parent(elt1), Parent(elt2), elt2);
    end if;

    W := CoxeterGroup(elt1`Parent);
    assocs := AssociativeArray(W);
    _AddScaled(~assocs, elt1`Terms, 1);
    _AddScaled(~assocs, elt2`Terms, 1);
    _RemoveZeros(~assocs);
    return EltIHkeConstruct(Parent(elt1), assocs);
end intrinsic;

intrinsic '-'(elt1::EltIHke, elt2::EltIHke) -> EltIHke
{The difference of two elements. Will throw an if subtracting over different free modules.}
    error if FreeModule(elt1) ne FreeModule(elt2),
        "Cannot subtract in different free modules", FreeModule(elt1), "and", FreeModule(elt2);

    // Change into the left basis if necessary.
    if Parent(elt1) ne Parent(elt2) then
        elt2 := ToBasis(Parent(elt1), Parent(elt2), elt2);
    end if;

    W := CoxeterGroup(elt1`Parent);
    assocs := AssociativeArray(W);
    _AddScaled(~assocs, elt1`Terms, 1);
    _AddScaled(~assocs, elt2`Terms, -1);
    _RemoveZeros(~assocs);
    return EltIHkeConstruct(Parent(elt1), assocs);
end intrinsic;

intrinsic '*'(elt::EltIHke, scalar::RngElt) -> EltIHke
{Scale an element by a scalar (Laurent polynomial or integer).}
    ok, r := IsCoercible(BaseRing(Parent(elt)), scalar);
    error if not ok,
        "Cannot scalar multiply by", scalar, "since the base ring of", Parent(elt), "is", BaseRing(Parent(elt));

    W := CoxeterGroup(elt`Parent);
    assocs := AssociativeArray(W);
    _AddScaled(~assocs, elt`Terms, r);
    _RemoveZeros(~assocs);
    return EltIHkeConstruct(Parent(elt), assocs);
end intrinsic;

intrinsic '*'(scalar::RngElt, elt::EltIHke) -> EltIHke
{Scale an element by a scalar (Laurent polynomial or integer).}
    return elt * scalar; // Defer to other implementation
end intrinsic;


////////////////////////////
// Multiplication (protocol)

intrinsic '*'(eltA::EltIHke, eltB::EltIHke) -> EltIHke
{}
    // In the case we are multiplying by zero, do nothing (might trigger a basis change otherwise).
    if IsZero(eltA) or IsZero(eltB) then
        return Parent(eltA) ! 0;
    end if;

    result := _Multiply(Parent(eltA), eltA, Parent(eltB), eltB);
    if Type(result) eq EltIHke then
        return result;
    end if;

    // Convert to standard basis and multiply.
    stdA := StandardBasis(FreeModule(eltA));
    stdB := StandardBasis(FreeModule(eltB));
    product := _Multiply(stdA, ToBasis(stdA, Parent(eltA), eltA), stdB, ToBasis(stdB, Parent(eltB), eltB));

    if Type(product) ne EltIHke then
        error "Multiplication not defined between", FreeModule(eltA), "and", FreeModule(eltB);
    end if;

    // Change back into either the A or B basis (we don't know which a priori, eg if this is a module action
    // rather than a multiplication). Prefer the left.
    if FreeModule(product) cmpeq FreeModule(eltA) then
        return ToBasis(Parent(eltA), Parent(product), product);
    elif FreeModule(product) cmpeq FreeModule(eltB) then
        return ToBasis(Parent(eltB), Parent(product), product);
    end if;

    error
        "Multiplication result in module", FreeModule(result),
        "was incompatible with either of", FreeModule(eltA),
        "or", FreeModule(eltB);
end intrinsic;

intrinsic _Multiply(A::BasisIHke, eltA::EltIHke, B::BasisIHke, eltB::EltIHke) -> EltIHke
{Fallback implementation for multiplication.}
    return false;
end intrinsic;


////////////////////////////
// Exponentiation

intrinsic '^'(elt::EltIHke, n::RngIntElt) -> EltIHke
{The nth power of elt.}
    require n ge 0: "The exponent", n, "must be nonnegative.";
    return &*[Parent(elt) | elt : _ in [1..n]];
end intrinsic;


////////////////////////////
// Bar involution (protocol)

intrinsic Bar(elt::EltIHke) -> EltIHke
{Perform the bar involution on elt, returning the result in the same basis as elt.}
    // First check if there is a bar involution defined on this basis.
    result := _Bar(Parent(elt), elt);
    if Type(result) eq EltIHke then return result; end if;

    // Convert to canonical basis, perform Bar, convert back.
    C := CanonicalBasis(FreeModule(elt));
    result := ToBasis(Parent(elt), C, Bar(ToBasis(C, Parent(elt), elt)));
    error if Type(result) ne EltIHke, "Bar involution not defined on", C;
    return result;
end intrinsic;

intrinsic _Bar(A::BasisIHke, eltA::EltIHke) -> EltIHke
{Fall-back implementation of Bar involution.}
    return false;
end intrinsic;
