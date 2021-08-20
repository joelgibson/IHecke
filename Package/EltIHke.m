import "Base.m":
    _LaurentPolyRing,
    _v;

// An EltIHke is a formal Laurent-polynomial-linear combination of Coxeter group elements. The
// interpretation of this linear combination depends on what the Parent structure is: it is either
// a type inheriting from AlgIHkeBase (a basis in the Hecke algebra), or ModIHke (a basis in a spherical
// or antispherical module).
declare type EltIHke;
declare attributes EltIHke:
    // An instance of a type inheriting from AlgIHkeBase or ModIHke.
    Parent,

    // An associative array of Coxeter group elements to Laurent polynomials.
    // This is always normalised so that there are no elements mapping to zero.
    Terms;


////////////////////////////
// Construction

// An EltIHke is always constructed via this function. Other files defining intrinsics should
// import this function.
function _EltIHkeConstruct(parent, terms)
    assert forall{w : w -> coeff in terms | coeff ne 0};
    assert Universe(terms) eq CoxeterGroup(parent);

    elt := New(EltIHke);
    elt`Parent := parent;
    elt`Terms := terms;
    return elt;
end function;


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
        printf "(%o)%o(%o)", _LaurentPolyRing!0, symbol, FmtElt(CoxeterGroup(elt`Parent).0);
    else
        printf Join([
                Sprintf("(%o)%o(%o)", elt`Terms[key], symbol, FmtElt(key))
                : key in Reverse(Sort(Setseq(Keys(elt`Terms))))
            ], " + ");
    end if;
end intrinsic;

intrinsic Parent(elt::EltIHke) -> .
{The parent structure (a type inheriting from AlgIHkeBase or ModIHke).}
    return elt`Parent;
end intrinsic;

intrinsic Coefficient(elt::EltIHke, w::GrpFPCoxElt) -> RngElt
{Return the coefficient of w in the linear combination.}
    if IsDefined(elt`Terms, w) then
        return elt`Terms[w];
    else
        return _LaurentPolyRing ! 0;
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
        elt2 := ChangeBasis(Parent(elt1), Parent(elt2), elt2);
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

intrinsic '+'(elt1::EltIHke, elt2::EltIHke) -> EltIHke
{The sum of two elements. Will throw an error if adding over different free modules.}
    error if Parent(Parent(elt1)) ne Parent(Parent(elt2)),
        "Cannot add in different free modules", Parent(Parent(elt1)), "and", Parent(Parent(elt2));

    // Change into the left basis if necessary.
    if Parent(elt1) ne Parent(elt2) then
        elt2 := ChangeBasis(Parent(elt1), Parent(elt2), elt2);
    end if;

    W := CoxeterGroup(elt1`Parent);
    assocs := AssociativeArray(W);
    _AddScaled(~assocs, elt1`Terms, 1);
    _AddScaled(~assocs, elt2`Terms, 1);
    _RemoveZeros(~assocs);
    return _EltIHkeConstruct(Parent(elt1), assocs);
end intrinsic;

intrinsic '-'(elt1::EltIHke, elt2::EltIHke) -> EltIHke
{The difference of two elements. Will throw an if subtracting over different free modules.}
    error if Parent(Parent(elt1)) ne Parent(Parent(elt2)),
        "Cannot subtract in different free modules", Parent(Parent(elt1)), "and", Parent(Parent(elt2));

    // Change into the left basis if necessary.
    if Parent(elt1) ne Parent(elt2) then
        elt2 := ChangeBasis(Parent(elt1), Parent(elt2), elt2);
    end if;

    W := CoxeterGroup(elt1`Parent);
    assocs := AssociativeArray(W);
    _AddScaled(~assocs, elt1`Terms, 1);
    _AddScaled(~assocs, elt2`Terms, -1);
    _RemoveZeros(~assocs);
    return _EltIHkeConstruct(Parent(elt1), assocs);
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
    return _EltIHkeConstruct(Parent(elt), assocs);
end intrinsic;

intrinsic '*'(scalar::RngElt, elt::EltIHke) -> EltIHke
{Scale an element by a scalar (Laurent polynomial or integer).}
    return elt * scalar; // Defer to other implementation
end intrinsic;


////////////////////////////
// Multiplication (protocol)

intrinsic '*'(eltA::EltIHke, eltB::EltIHke) -> EltIHke
{}
    result := _IHkeProtMult(Parent(eltA), eltA, Parent(eltB), eltB);
    if Type(result) eq EltIHke then
        return result;
    end if;

    // Convert to standard basis and multiply.
    H := IHeckeAlgebraStd(Parent(Parent(eltA)));
    product := _IHkeProtMult(H, ChangeBasis(H, Parent(eltA), eltA), H, ChangeBasis(H, Parent(eltB), eltB));
    return ChangeBasis(Parent(eltA), H, product);
end intrinsic;

intrinsic _IHkeProtMult(A::AlgIHkeBase, eltA::EltIHke, B::AlgIHkeBase, eltB::EltIHke) -> EltIHke
{Fallback implementation for multiplication.}
    return false;
end intrinsic;


////////////////////////////
// Bar involution (protocol)

intrinsic Bar(elt::EltIHke) -> EltIHke
{Perform the bar involution on elt, returning the result in the same basis as elt.}
    return _IHkeProtBar(Parent(elt), elt);
end intrinsic;

intrinsic _IHkeProtBar(A::AlgIHkeBase, elt::EltIHke) -> EltIHke
{Fall-back implementation of the bar involution: convert to canonical, apply involution, convert back.}
    C := IHeckeAlgebraCan(Parent(A));
    return ChangeBasis(A, C, Bar(ChangeBasis(C, A, elt)));
end intrinsic;


/////////////////////
// Input

// TODO: Move to p-can stuff? Re-evaluate at least.
intrinsic ReadEltIHke(parent::., eltStr::MonStgElt) -> EltIHke
{Read a printed Hecke element. The symbol used for printing the basis ("H", "C", etc) are ignored,
 and instead the given parent is used as the basis in which to interpret the string.}
    W := CoxeterGroup(parent);
    v := _v; // Needed for evalling the coefficient.
    terms := AssociativeArray(W);
    line := eltStr;
    while true do
        ok, _, matches := Regexp("^\\(([^)]+)\\)[^(]*\\(([^)]+)\\)( [+] (.*))?$", eltStr);
        if not ok then
            error if not Regexp("^[ ]*$", eltStr), "Could not parse the element", line;
        end if;

        w := matches[2] eq "id"
            select W.0
            else W ! [StringToInteger(matches[2][i]) : i in [1..#matches[2]]];

        _AddScaledTerm(~terms, w, eval matches[1]);
        if #matches eq 2 then
            break;
        end if;
        eltStr := matches[4];
    end while;
    _RemoveZeros(~terms);
    return _EltIHkeConstruct(parent, terms);
end intrinsic;

