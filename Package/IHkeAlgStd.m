// This file implements the standard bases of the Hecke algebra, antispherical module, and spherical
// module. We'll call any one of these three, considered as a right module, a "parabolic module".
// This treats multiplication and the bar involution somewhat uniformly across the three.

import "EltIHke.m": _AddScaled, _AddScaledTerm, _RemoveZeros, _IsUniTriangular;

// Abstract type, used for the standard basis of all three modules.
declare type IHkeAlgBaseStd: BasisIHke;
declare attributes IHkeAlgBaseStd:
    BarCache,   // An associative array, used as a cache for the bar involution on basis elements
    Para,       // Parabolic subset (equals [] for the full Hecke algebra)
    Eig;        // Irrelevant for the full Hecke algebra, (-v) or (v^-1) for the anti/spherical mods

// Factory function for basis types extending IHkeAlgBaseStd.
function _GetOrCreateBasis(fmod, basisType, symbol, name, para, eig)
    assert ISA(basisType, IHkeAlgBaseStd);
    if not IsDefined(fmod`BasisCache, basisType) then
        basis := New(basisType);
        _BasisIHkeInit(~basis, fmod, symbol, name);
        basis`BarCache := AssociativeArray();
        basis`Para := para;
        basis`Eig := eig;
        fmod`BasisCache[basisType] := basis;
    end if;

    return fmod`BasisCache[basisType];
end function;

// (M, eltM) is an element inside a right module over the Hecke algebra (either the algebra itself,
// or a parabolic module). Return the right action of H(s) on eltM.
function _RightMultStdGen(M, eltM, s)
    W := CoxeterGroup(M);
    v := BaseRing(M).1;
    terms := AssociativeArray(W);
    for w -> coeff in eltM`Terms do
        ws := w * (W.s);
        if not IsMinimal(M`Para, ws) then
            _AddScaledTerm(~terms, w, M`Eig * coeff);
        elif #w lt #ws then
            _AddScaledTerm(~terms, ws, coeff);
        else
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, (v^-1 - v) * coeff);
        end if;
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(M, terms);
end function;

// (M, eltM) is an element inside a right module over the Hecke algebra (either the algebra itself,
// or a parabolic module). This function returns eltM * H(terms), where terms is interpreted as a
// linear combination of standard basis elements of the Hecke algebra.
function _RightAction(M, eltM, terms)
    acc := AssociativeArray(CoxeterGroup(M));
    for w -> coeff in terms do
        piece := eltM;
        for s in Eltseq(w) do
            piece := _RightMultStdGen(M, piece, s);
        end for;
        _AddScaled(~acc, piece`Terms, coeff);
    end for;
    _RemoveZeros(~acc);
    return EltIHkeConstruct(M, acc);
end function;

// Return the bar-involution of the basis element M(w), i.e. the image of the element H(w^-1)^-1
// inside the module M. (Pick a right descent s, and use H(w^-1)^-1 = H(ws^-1)^-1 * H(s)^-1).
function _BarInvolutionStd(M, w)
    if IsDefined(M`BarCache, w) then
        return M`BarCache[w];
    end if;

    if #w eq 0 then
        return M.0;
    end if;

    W := CoxeterGroup(M);
    v := BaseRing(M).1;
    s := Rep(RightDescentSet(W, w));
    bar_ws := _BarInvolutionStd(M, w * W.s);
    bar_w := _RightMultStdGen(M, bar_ws, s) + (v - v^-1)*bar_ws;
    assert _IsUniTriangular(bar_w`Terms, w);
    M`BarCache[w] := bar_w;
    return bar_w;
end function;

intrinsic _IHkeProtBar(H::IHkeAlgBaseStd, elt::EltIHke) -> EltIHke
{The bar involution on elt, mapping p(v)H(w) to p(v^-1)H(w^-1)^-1.}
    assert H eq Parent(elt);

    R := BaseRing(H);
    twist := hom<R -> R | (R.1)^-1>;
    terms := AssociativeArray(CoxeterGroup(H));
    for w -> coeff in elt`Terms do
        _AddScaled(~terms, _BarInvolutionStd(H, w)`Terms, twist(coeff));
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(H, terms);
end intrinsic;


////////////////////////////////////////////
// Standard basis of the Hecke algebra

declare type IHkeAlgStd[EltIHke]: IHkeAlgBaseStd;

intrinsic IHeckeAlgebraStd(alg::IHkeAlg) -> IHkeAlgStd
{The standard basis of the Hecke algebra.}
    return _GetOrCreateBasis(alg, IHkeAlgStd, "H", "Standard basis", [], 0);
end intrinsic;

intrinsic _IHkeProtUnit(H::IHkeAlgStd) -> EltIHke
{Unit element in the standard basis.}
    return H.0;
end intrinsic;

intrinsic _IHkeProtMult(H1::IHkeAlgStd, elt1::EltIHke, H2::IHkeAlgStd, elt2::EltIHke) -> EltIHke
{Multiplication inside the standard basis of the Hecke algebra.}
    return _RightAction(H1, elt1, elt2`Terms);
end intrinsic;


/////////////////////////////////////////////
// Standard basis of the antispherical module

declare type ASModIHkeStd[EltIHke]: IHkeAlgBaseStd;

intrinsic IHeckeAntiSphericalStd(ASMod::IHkeASMod) -> ASModIHkeStd
{The standard basis of the right spherical module.}
    v := BaseRing(ASMod).1;
    return _GetOrCreateBasis(ASMod, ASModIHkeStd, "aH", "Standard basis", Parabolic(ASMod), -v);
end intrinsic;

intrinsic _EltIHkeValidate(aH::ASModIHkeStd, elt::EltIHke)
{Only allow I-minimal elements.}
    I := Parabolic(Parent(aH));
    error if not forall(w){w : w -> _ in elt`Terms | IsMinimal(I, w)},
        w, "is not minimal with respect to", I;
end intrinsic;

intrinsic _IHkeProtMult(M::ASModIHkeStd, eltM::EltIHke, H::IHkeAlgStd, eltH::EltIHke) -> EltIHke
{Action of the standard basis of the Hecke algebra on the standard basis of the antispherical module.}
    return _RightAction(M, eltM, eltH`Terms);
end intrinsic;


/////////////////////////////////////////
// Standard basis of the spherical module

declare type SModIHkeStd[EltIHke]: IHkeAlgBaseStd;

intrinsic IHeckeSphericalStd(SMod::IHkeSMod) -> SModIHkeStd
{The standard basis of the right spherical module.}
    v := BaseRing(SMod).1;
    return _GetOrCreateBasis(SMod, SModIHkeStd, "sH", "Standard basis", Parabolic(SMod), v^-1);
end intrinsic;

intrinsic _EltIHkeValidate(sH::SModIHkeStd, elt::EltIHke)
{Only allow I-minimal elements.}
    I := Parabolic(Parent(sH));
    error if not forall(w){w : w -> _ in elt`Terms | IsMinimal(I, w)},
        w, "is not minimal with respect to", I;
end intrinsic;

intrinsic _IHkeProtMult(M::SModIHkeStd, eltM::EltIHke, H::IHkeAlgStd, eltH::EltIHke) -> EltIHke
{Action of the standard basis of the Hecke algebra on the standard basis of the spherical module.}
    return _RightAction(M, eltM, eltH`Terms);
end intrinsic;