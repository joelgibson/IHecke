import "EltIHke.m": _AddScaled, _AddScaledTerm, _RemoveZeros, _IsUniTriangular;

// Abstract type, used for the canonical basis of all three modules.
declare type IHkeAlgBaseCan: BasisIHke;
declare attributes IHkeAlgBaseCan:
    CanInStdCache,
    StdInCanCache,
    MuCache,
    Para,
    Eig;

// Factory function for basis types extending IHkeAlgBaseCan.
function _GetOrCreateBasis(fmod, basisType, symbol, name, para, eig)
    assert ISA(basisType, IHkeAlgBaseCan);
    if not IsDefined(fmod`BasisCache, basisType) then
        basis := New(basisType);
        _BasisIHkeInit(~basis, fmod, symbol, name);
        basis`CanInStdCache := AssociativeArray();
        basis`StdInCanCache := AssociativeArray();
        basis`MuCache := AssociativeArray();
        basis`Para := para;
        basis`Eig := eig;
        fmod`BasisCache[basisType] := basis;
    end if;

    return fmod`BasisCache[basisType];
end function;

// Assume that elt is an element written in terms of the standard basis, and right-multiply
// by the element B(s). I is the parabolic quotient ([] if we are in the Hecke algebra), and eig
// is (-v) for antispherical, v^-1 for spherical, and irrelevant if in the Hecke algebra.
function _RightMultCanGen(elt, s, I, eig)
    W := CoxeterGroup(Parent(elt));
    v := BaseRing(Parent(elt)).1;
    terms := AssociativeArray();
    for w -> coeff in elt`Terms do
        ws := w * (W.s);
        if #I ne 0 and not IsMinimal(I, ws) then
            _AddScaledTerm(~terms, w, (v + eig) * coeff);
        elif #w lt #ws then
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, v * coeff);
        else
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, v^-1 * coeff);
        end if;
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(Parent(elt), terms);
end function;

// An associative array mapping group elements u to the coefficient mu(u, w).
function _MuCoeffs(H, C, w)
    if IsDefined(C`MuCache, w) then
        return C`MuCache[w];
    end if;

    W := CoxeterGroup(C);
    if #w eq 0 then
        return AssociativeArray(W);
    end if;

    Cw := _IHkeProtToBasis(H, C, w);
    mu := AssociativeArray(W);
    for u -> coeff in Cw`Terms do
        v_coeff := Coefficient(coeff, 1);
        if v_coeff ne 0 then
            mu[u] := v_coeff;
        end if;
    end for;

    C`MuCache[w] := mu;
    return mu;
end function;

// Express C(w) in the standard basis.
function _CanToStd(H, C, w, I, eig)
    if IsDefined(C`CanInStdCache, w) then
        return C`CanInStdCache[w];
    end if;

    if #w eq 0 then
        return H.0;
    end if;

    // Take an arbitrary right descent of w.
    W := CoxeterGroup(H);
    s := Rep(RightDescentSet(W, w));
    ws := w * W.s;

    // C(ws) = C(w)C(s) - sum[us < u]mu(u, ws) C(u)
    Cw := _IHkeProtToBasis(H, C, ws);
    CwCs := _RightMultCanGen(Cw, s, I, eig);
    for u -> coeff in _MuCoeffs(H, C, ws) do
        if #(u*W.s) lt #u then
            CwCs -:= coeff * _IHkeProtToBasis(H, C, u);
        end if;
    end for;

    // Double-check unitriangularity
    assert2 _IsUniTriangular(CwCs`Terms, w);

    C`CanInStdCache[w] := CwCs;
    return CwCs;
end function;

// Express H(w) in the canonical basis, using unitriangularity of C in H.
function _StdToCan(H, C, w)
    if IsDefined(C`StdInCanCache, w) then
        return C`StdInCanCache[w];
    end if;

    if #w eq 0 then
        return C.0;
    end if;

    W := CoxeterGroup(H);
    terms := AssociativeArray(W);
    terms[w] := BaseRing(C) ! 1;
    Cw := _IHkeProtToBasis(H, C, w);
    for u -> coeff in Cw`Terms do
        if u eq w then
            continue;
        end if;

        Cu := _IHkeProtToBasis(C, H, u);
        _AddScaled(~terms, Cu`Terms, -coeff);
    end for;
    _RemoveZeros(~terms);
    result := EltIHkeConstruct(C, terms);

    C`StdInCanCache[w] := result;
    return result;
end function;

intrinsic _IHkeProtBar(C::IHkeAlgBaseCan, elt::EltIHke) -> EltIHke
{The bar involution of elt, fixing each basis element and twisting scalars by by v -> v^-1.}
    assert C eq Parent(elt);

    R := BaseRing(C);
    twist := hom<R -> R | (R.1)^-1>;
    terms := AssociativeArray(CoxeterGroup(C));
    for w -> coeff in elt`Terms do
        terms[w] := twist(coeff);
    end for;

    return EltIHkeConstruct(C, terms);
end intrinsic;


///////////////////////////////////////
// Canonical basis of the Hecke algebra

declare type IHkeAlgCan[EltIHke]: IHkeAlgBaseCan;

intrinsic IHeckeAlgebraCan(HAlg::IHkeAlg) -> IHkeAlgCan
{The canonical basis of the Hecke algebra.}
    return _GetOrCreateBasis(HAlg, IHkeAlgCan, "C", "Canonical basis", [], 0);
end intrinsic;

intrinsic _IHkeProtToBasis(H::IHkeAlgStd, C::IHkeAlgCan, w::GrpFPCoxElt) -> EltIHke
{Express C(w) in the standard basis.}
    return _CanToStd(H, C, w, [], 0);
end intrinsic;

intrinsic _IHkeProtToBasis(C::IHkeAlgCan, H::IHkeAlgStd, w::GrpFPCoxElt) -> EltIHke
{Express H(w) in the canonical basis.}
    return _StdToCan(H, C, w);
end intrinsic;


//////////////////////////////////////////////
// Canonical basis of the antispherical module

declare type ASModIHkeCan[EltIHke]: IHkeAlgBaseCan;

intrinsic IHeckeAntiSphericalCan(ASMod::IHkeASMod) -> ASModIHkeCan
{The canonical basis of the right antispherical module.}
    v := BaseRing(ASMod).1;
    return _GetOrCreateBasis(ASMod, ASModIHkeCan, "aC", "Canonical basis", Parabolic(ASMod), -v);
end intrinsic;

intrinsic _EltIHkeValidate(aC::ASModIHkeCan, elt::EltIHke)
{Only allow I-minimal elements.}
    I := Parabolic(Parent(aC));
    error if not forall(w){w : w -> _ in elt`Terms | IsMinimal(I, w)},
        w, "is not minimal with respect to", I;
end intrinsic;

intrinsic _IHkeProtToBasis(aH::ASModIHkeStd, aC::ASModIHkeCan, w::GrpFPCoxElt) -> EltIHke
{Express aC(w) in the standard basis.}
    return _CanToStd(aH, aC, w, aH`Para, aH`Eig);
end intrinsic;

intrinsic _IHkeProtToBasis(aC::ASModIHkeCan, aH::ASModIHkeStd, w::GrpFPCoxElt) -> EltIHke
{Express aH(w) in the canonical basis.}
    return _StdToCan(aH, aC, w);
end intrinsic;


//////////////////////////////////////////
// Canonical basis of the spherical module

declare type SModIHkeCan[EltIHke]: IHkeAlgBaseCan;

intrinsic IHeckeSphericalCan(SMod::IHkeSMod) -> SModIHkeCan
{The canonical basis of the right antispherical module.}
    v := BaseRing(SMod).1;
    return _GetOrCreateBasis(SMod, SModIHkeCan, "sC", "Canonical basis", Parabolic(SMod), v^-1);
end intrinsic;

intrinsic _EltIHkeValidate(sC::SModIHkeCan, elt::EltIHke)
{Only allow I-minimal elements.}
    I := Parabolic(Parent(sC));
    error if not forall(w){w : w -> _ in elt`Terms | IsMinimal(I, w)},
        w, "is not minimal with respect to", I;
end intrinsic;

intrinsic _IHkeProtToBasis(sH::SModIHkeStd, sC::SModIHkeCan, w::GrpFPCoxElt) -> EltIHke
{Express sC(w) in the standard basis.}
    return _CanToStd(sH, sC, w, Parabolic(Parent(sC)), (BaseRing(sC).1)^-1);
end intrinsic;

intrinsic _IHkeProtToBasis(sC::SModIHkeCan, sH::SModIHkeStd, w::GrpFPCoxElt) -> EltIHke
{Express sH(w) in the canonical basis.}
    return _StdToCan(sH, sC, w);
end intrinsic;