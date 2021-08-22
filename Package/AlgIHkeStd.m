import "Base.m":
    _LaurentPolyRing,
    _v;

import "EltIHke.m": _AddScaled, _RemoveZeros, _AddScaledTerm;
import "AlgIHkeBase.m": _AlgIHkeBaseInit;


//////////////////////////////////////
// Standard basis of the Hecke algebra

declare type AlgIHkeStd[EltIHke]: AlgIHkeBase;
declare attributes AlgIHkeStd:
    // An associative array, used as a cache for the bar involution on basis elements.
    BarCache;

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

intrinsic _IHkeProtUnit(H::AlgIHkeStd) -> EltIHke
{}
    return H.0;
end intrinsic;


//////////////////////////////////
// Multiplication

// Right-multiply an element in the standard basis by H(s).
// Whenever this is called for the Hecke algebra, we will have I=[] and eig is irrelevant.
// In the antispherical or spherical modules, I is the parabolic quotient, and eig is (-v) or v^-1.
function _RightMultStdGen(elt, s, I, eig)
    W := CoxeterGroup(Parent(elt));
    terms := AssociativeArray();
    for w -> coeff in elt`Terms do
        ws := w * (W.s);
        if #I ne 0 and not IsMinimal(I, ws) then
            _AddScaledTerm(~terms, w, eig * coeff);
        elif #w lt #ws then
            _AddScaledTerm(~terms, ws, coeff);
        else
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, (_v^-1 - _v) * coeff);
        end if;
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(Parent(elt), terms);
end function;

intrinsic _IHkeProtMult(A::AlgIHkeStd, eltA::EltIHke, B::AlgIHkeStd, eltB::EltIHke) -> EltIHke
{Multiplication in the standard basis.}
    require A eq B: "Parents must be equal";

    acc := AssociativeArray();
    for w -> coeff in eltB`Terms do
        piece := eltA;
        for s in Eltseq(w) do
            piece := _RightMultStdGen(piece, s, [], 0);
        end for;
        _AddScaled(~acc, piece`Terms, coeff);
    end for;
    _RemoveZeros(~acc);
    return EltIHkeConstruct(A, acc);
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
    return EltIHkeConstruct(Parent(elt), terms);
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
    return EltIHkeConstruct(H, terms);
end intrinsic;


/////////////////////////////////////////////
// Standard basis of the antispherical module

declare type ASModIHkeStd[EltIHke]: AlgIHkeBase;
declare attributes ASModIHkeStd:
    // An associative array, used as a cache for the bar involution on basis elements.
    BarCache;

intrinsic IHeckeAntiSphericalStd(asmod::ASphIHke) -> ASModIHkeStd
{The standard basis of the right antispherical module.}
    if not IsDefined(asmod`BasisCache, ASModIHkeStd) then
        basis := New(ASModIHkeStd);
        _AlgIHkeBaseInit(~basis, asmod, "aH", "Standard basis");
        basis`BarCache := AssociativeArray(CoxeterGroup(asmod));
        asmod`BasisCache[ASModIHkeStd] := basis;
    end if;
    return asmod`BasisCache[ASModIHkeStd];
end intrinsic;

intrinsic _IHkeProtValidateElt(aH::ASModIHkeStd, elt::EltIHke)
{Only allow I-minimal elements.}
    I := Parabolic(Parent(aH));
    error if not forall(w){w : w -> _ in elt`Terms | IsMinimal(I, w)},
        w, "is not minimal with respect to", I;
end intrinsic;

intrinsic _IHkeProtMult(aH::ASModIHkeStd, asElt::EltIHke, H::AlgIHkeStd, hElt::EltIHke) -> EltIHke
{Action of the standard basis of the Hecke algebra on the standard basis of the antispherical module.}
    require Parent(aH)`CoxeterMatrix eq Parent(H)`CoxeterMatrix: "incompatible module";

    I := Parabolic(Parent(aH));
    acc := AssociativeArray();
    for w -> coeff in hElt`Terms do
        piece := asElt;
        for s in Eltseq(w) do
            piece := _RightMultStdGen(piece, s, I, -_v);
        end for;
        _AddScaled(~acc, piece`Terms, coeff);
    end for;
    _RemoveZeros(~acc);
    return EltIHkeConstruct(aH, acc);
end intrinsic;


/////////////////////////////////////////
// Standard basis of the spherical module

declare type SModIHkeStd[EltIHke]: AlgIHkeBase;
declare attributes SModIHkeStd:
    // An associative array, used as a cache for the bar involution on basis elements.
    BarCache;

intrinsic IHeckeSphericalStd(smod::SphIHke) -> SModIHkeStd
{The standard basis of the right spherical module.}
    if not IsDefined(smod`BasisCache, SModIHkeStd) then
        basis := New(SModIHkeStd);
        _AlgIHkeBaseInit(~basis, smod, "sH", "Standard basis");
        basis`BarCache := AssociativeArray(CoxeterGroup(smod));
        smod`BasisCache[SModIHkeStd] := basis;
    end if;
    return smod`BasisCache[SModIHkeStd];
end intrinsic;

intrinsic _IHkeProtValidateElt(sH::SModIHkeStd, elt::EltIHke)
{Only allow I-minimal elements.}
    I := Parabolic(Parent(sH));
    error if not forall(w){w : w -> _ in elt`Terms | IsMinimal(I, w)},
        w, "is not minimal with respect to", I;
end intrinsic;

intrinsic _IHkeProtMult(sH::SModIHkeStd, sElt::EltIHke, H::AlgIHkeStd, hElt::EltIHke) -> EltIHke
{Action of the standard basis of the Hecke algebra on the standard basis of the spherical module.}
    require Parent(sH)`CoxeterMatrix eq Parent(H)`CoxeterMatrix: "incompatible module";

    I := Parabolic(Parent(sH));
    acc := AssociativeArray();
    for w -> coeff in hElt`Terms do
        piece := sElt;
        for s in Eltseq(w) do
            piece := _RightMultStdGen(piece, s, I, _v^-1);
        end for;
        _AddScaled(~acc, piece`Terms, coeff);
    end for;
    _RemoveZeros(~acc);
    return EltIHkeConstruct(sH, acc);
end intrinsic;