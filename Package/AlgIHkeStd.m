import "Base.m":
    _LaurentPolyRing,
    _v;

import "EltIHke.m":  _EltIHkeConstruct, _AddScaled, _RemoveZeros, _AddScaledTerm;
import "AlgIHkeBase.m": _AlgIHkeBaseInit;


declare type AlgIHkeStd[EltIHke]: AlgIHkeBase;
declare attributes AlgIHkeStd:
    // An associative array, used as a cache for the bar involution on basis elements.
    BarCache;


//////////////////////////////////
// Creation

intrinsic IHeckeAlgebraStd(alg::AlgIHke) -> AlgIHkeStd
{The standard basis of the Hecke algebra.}
    if not IsDefined(alg`BasisCache, AlgIHkeStd) then
        basis := New(AlgIHkeStd);
        _AlgIHkeBaseInit(~basis, alg, "H", "Standard basis");
        basis`BarCache := AssociativeArray(CoxeterGroup(alg));
        alg`BasisCache[AlgIHkeStd] := basis;
    end if;
    return alg`BasisCache[AlgIHkeStd];
end intrinsic;


/////////////////////////////////
// Overrides

intrinsic 'eq'(H1::AlgIHkeStd, H2::AlgIHkeStd) -> BoolElt
{}
    // If we have two bases of type AlgIHkeStd of the same Hecke algebra, they must be equal.
    return Parent(H1) eq Parent(H2);
end intrinsic;

intrinsic _IHkeProtUnit(H::AlgIHkeStd) -> EltIHke
{}
    return H.0;
end intrinsic;



//////////////////////////////////
// Multiplication

// Right-multiply an element in the standard basis by H(s).
function _RightMultStdGen(elt, s)
    assert ISA(Type(Parent(elt)), AlgIHkeStd);

    W := CoxeterGroup(Parent(elt));
    terms := AssociativeArray();
    for w -> coeff in elt`Terms do
        ws := w * (W.s);
        if #w lt #ws then
            _AddScaledTerm(~terms, ws, coeff);
        else
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, (_v^-1 - _v) * coeff);
        end if;
    end for;
    _RemoveZeros(~terms);
    return _EltIHkeConstruct(Parent(elt), terms);
end function;

intrinsic _IHkeProtMult(A::AlgIHkeStd, eltA::EltIHke, B::AlgIHkeStd, eltB::EltIHke) -> EltIHke
{Multiplication in the standard basis.}
    require A eq B: "Parents must be equal";

    acc := AssociativeArray();
    for w -> coeff in eltB`Terms do
        piece := eltA;
        for s in Eltseq(w) do
            piece := _RightMultStdGen(piece, s);
        end for;
        _AddScaled(~acc, piece`Terms, coeff);
    end for;
    _RemoveZeros(~acc);
    return _EltIHkeConstruct(A, acc);
end intrinsic;


////////////////////////////
// Bar involution (protocol)

// Right-multiply an element in the standard basis by H(s)^-1.
function _RightMultStdGenInv(elt, s)
    assert ISA(Type(Parent(elt)), AlgIHkeStd);

    W := CoxeterGroup(Parent(elt));
    terms := AssociativeArray();
    for w -> coeff in elt`Terms do
        ws := w * (W.s);
        if #w lt #ws then
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, (_v - _v^-1) * coeff);
        else
            _AddScaledTerm(~terms, ws, coeff);
        end if;
    end for;
    _RemoveZeros(~terms);
    return _EltIHkeConstruct(Parent(elt), terms);
end function;

// Return the bar-involution of the basis element H(w), i.e. the element H(w^-1)^-1.
// Do this by picking a right descent s, and using that H(w^-1)^-1 = H(ws^-1)^-1 * H(s)^-1.
function _BarInvolutionStd(H, w)
    assert ISA(Type(H), AlgIHkeStd);

    if IsDefined(H`BarCache, w) then
        return H`BarCache[w];
    end if;

    if #w eq 0 then
        return H.0;
    end if;

    W := CoxeterGroup(H);
    s := Rep(RightDescentSet(W, w));
    H`BarCache[w] := _RightMultStdGenInv(_BarInvolutionStd(H, w * W.s), s);
    return H`BarCache[w];
end function;

intrinsic _EltIHkeBar(H::AlgIHkeStd, elt::EltIHke) -> EltIHke
{The bar involution on elt, mapping p(v)H(w) to p(v^-1)H(w^-1)^-1.}
    assert H eq Parent(elt);

    twist := hom<_LaurentPolyRing -> _LaurentPolyRing | _v^-1>;
    terms := AssociativeArray(CoxeterGroup(H));
    for w -> coeff in elt`Terms do
        _AddScaled(~terms, _BarInvolutionStd(H, w)`Terms, twist(coeff));
    end for;
    _RemoveZeros(~terms);
    return _EltIHkeConstruct(H, terms);
end intrinsic;