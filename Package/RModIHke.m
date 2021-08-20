import "Base.m":
    _LaurentPolyRing;
import "EltIHke.m":
    _EltIHkeConstruct,
    _AddScaled,
    _AddScaledTerm,
    _RemoveZeros;

declare type RModIHke;
declare attributes RModIHke:
    // An AlgIHke
    Parent,

    // A sequence of integers giving the parabolic quotient.
    Para;


declare type RModIHkeBase;
declare attributes RModIHkeBase:
    // AlgIHke
    Parent,

    // A short string naming the basis, used for printing.
    BasisSymbol,

    // A human-readable name describing the basis.
    BasisName,

    // The sequence of integers giving the parabolic quotient. This determines the set of Coxeter
    // group elements which is "allowed".
    Para;


procedure _RModIHeckeBaseInit(~basis, alg, basisSymbol, basisName, para)
    assert ISA(Type(alg), AlgIHke);
    assert ISA(Type(basis), RModIHkeBase);
    assert Type(basisSymbol) eq MonStgElt;
    assert Type(basisName) eq MonStgElt;
    assert para subset [1..Rank(CoxeterGroup(alg))];

    basis`Parent := alg;
    basis`BasisSymbol := basisSymbol;
    basis`BasisName := basisName;
    basis`Para := Sort(Setseq({s : s in para}));
end procedure;


/////////////////////
// Accessor functions

intrinsic Print(M::RModIHkeBase)
{}
    printf "%o of right module over Hecke algebra for Coxeter group of type %o, symbol %o, parabolic quotient %o",
        BasisName(M),
        CartanName(CoxeterGroup(M)),
        BasisSymbol(M),
        Parabolic(M);
end intrinsic;

intrinsic Parent(M::RModIHkeBase) -> AlgIHke
{The parent Hecke algebra.}
    return M`Parent;
end intrinsic;

intrinsic CoxeterGroup(M::RModIHkeBase) -> GrpFPCox
{The underlying Coxeter group of the Hecke algebra.}
    return M`Parent`CoxeterGroup;
end intrinsic;

intrinsic BasisSymbol(M::RModIHkeBase) -> MonStgElt
{A short string naming the basis.}
    return M`BasisSymbol;
end intrinsic;

intrinsic BasisName(M::RModIHkeBase) -> MonStgElt
{A human-readable name describing the basis.}
    return M`BasisName;
end intrinsic;

intrinsic BaseRing(M::RModIHkeBase) -> Rng
{Return the base ring (a Laurent series ring).}
    return _LaurentPolyRing;
end intrinsic;

// TODO: Is this the right name?
intrinsic Parabolic(M::RModIHkeBase) -> SeqEnum[RngIntElt]
{}
    return M`Para;
end intrinsic;


/////////////////////////////
// Creation of basis elements
//
// The difference from AlgIHkeBase is that we only allow elements which are minimal in their right
// W_I coset, where I is the generating set of the parabolic subgroup.

intrinsic '.'(M::RModIHkeBase, w::GrpFPCoxElt) -> EltIHke
{The basis element indexed by the Coxeter group element w.}
    require Parent(w) eq CoxeterGroup(M):
        "Group element", w, "is not a member of", CoxeterGroup(M);
    require IsMinimal(Parabolic(M), w):
        "Group element", w, "is not minimal in W_I*x with I =", Parabolic(M);

    terms := AssociativeArray(CoxeterGroup(M));
    terms[w] := _LaurentPolyRing ! 1;
    return _EltIHkeConstruct(M, terms);
end intrinsic;

intrinsic '.'(M::RModIHkeBase, s::RngIntElt) -> EltIHke
{The basis element indexed by the Coxeter group element W.s, equivalent to A.(W.s).}
    require s notin Parabolic(M):
        "The generator", s, "is in the quotient group W_I with I =", Parabolic(M);

    w := CoxeterGroup(M) . s;
    terms := AssociativeArray(CoxeterGroup(M));
    terms[w] := _LaurentPolyRing ! 1;
    return _EltIHkeConstruct(M, terms);
end intrinsic;

intrinsic '.'(M::RModIHkeBase, word::SeqEnum[RngIntElt]) -> EltIHke
{The basis element indexed by the Coxeter group word, equivalent to A.(W!word).}
    w := CoxeterGroup(M) ! word;
    require IsMinimal(Parabolic(M), w):
        "Group element", w, "is not minimal in W_I*x with I =", Parabolic(M);

    terms := AssociativeArray(CoxeterGroup(M));
    terms[w] := _LaurentPolyRing ! 1;
    return _EltIHkeConstruct(M, terms);
end intrinsic;


//////////////////////
// Equality (protocol)

intrinsic 'eq'(M::RModIHkeBase, N::RModIHkeBase) -> BoolElt
{}
    return false;
end intrinsic;

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
