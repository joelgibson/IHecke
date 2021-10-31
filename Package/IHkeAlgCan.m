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

    Cw := _ToBasis(H, C, w);
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
    Cw := _ToBasis(H, C, ws);
    CwCs := _RightMultCanGen(Cw, s, I, eig);
    for u -> coeff in _MuCoeffs(H, C, ws) do
        if #(u*W.s) lt #u then
            CwCs -:= coeff * _ToBasis(H, C, u);
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
    Cw := _ToBasis(H, C, w);
    for u -> coeff in Cw`Terms do
        if u eq w then
            continue;
        end if;

        Cu := _ToBasis(C, H, u);
        _AddScaled(~terms, Cu`Terms, -coeff);
    end for;
    _RemoveZeros(~terms);
    result := EltIHkeConstruct(C, terms);

    C`StdInCanCache[w] := result;
    return result;
end function;

intrinsic _Bar(C::IHkeAlgBaseCan, elt::EltIHke) -> EltIHke
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

intrinsic CanonicalBasis(HAlg::IHkeAlg) -> IHkeAlgCan
{The canonical basis of the Hecke algebra.}
    return _GetOrCreateBasis(HAlg, IHkeAlgCan, "C", "Canonical basis", [], 0);
end intrinsic;

intrinsic _ToBasis(H::IHkeAlgStd, C::IHkeAlgCan, w::GrpFPCoxElt) -> EltIHke
{Express C(w) in the standard basis.}
    return _CanToStd(H, C, w, [], 0);
end intrinsic;

intrinsic _ToBasis(C::IHkeAlgCan, H::IHkeAlgStd, w::GrpFPCoxElt) -> EltIHke
{Express H(w) in the canonical basis.}
    return _StdToCan(H, C, w);
end intrinsic;

intrinsic _Multiply(C::IHkeAlgCan, eltA::EltIHke, B::IHkeAlgCan, eltB::EltIHke) -> EltIHke
{Intercept left or right multiplication by the identity or C(s), and implement it in terms
 of the mu-coefficients. Otherwise return false (defers to standard multiplication).}
    assert C eq B;
    H := StandardBasis(FreeModule(C));
    W := CoxeterGroup(C);
    L := BaseRing(C);

    // Two special cases: multiplication by the identity, or by C(s) for a simple generator s.
    // For multiplication by C(s) we have the formulas
    //   C(s) C(w) = (v + v^-1) C(w) if sw < w.
    //   C(s) C(w) = C(sw) + sum[z | sz < z]mu(z, w)C(z).

    if #eltA`Terms eq 1 then
        s := Rep(Keys(eltA`Terms));
        scale := eltA`Terms[s];

        if #s eq 0 then
            return scale * eltB;
        end if;
        if #s eq 1 then
            terms := AssociativeArray(W);
            for w -> coeff in eltB`Terms do
                if #(s*w) lt #w then
                    _AddScaledTerm(~terms, w, scale * coeff * (L.1 + (L.1)^-1));
                else
                    _AddScaledTerm(~terms, s*w, scale * coeff);
                    mu := _MuCoeffs(H, C, w);
                    for z -> mucoeff in mu do
                        if #(s*z) lt #z then
                            _AddScaledTerm(~terms, z, mucoeff * scale * coeff);
                        end if;
                    end for;
                end if;
            end for;
            _RemoveZeros(~terms);
            return EltIHkeConstruct(C, terms);
        end if;
    end if;

    if #eltB`Terms eq 1 then
        s := Rep(Keys(eltB`Terms));
        scale := eltB`Terms[s];

        if #s eq 0 then
            return scale * eltA;
        end if;
        if #s eq 1 then
            terms := AssociativeArray(W);
            for w -> coeff in eltA`Terms do
                if #(w*s) lt #w then
                    _AddScaledTerm(~terms, w, scale * coeff * (L.1 + (L.1)^-1));
                else
                    _AddScaledTerm(~terms, w*s, scale * coeff);
                    mu := _MuCoeffs(H, C, w);
                    for z -> mucoeff in mu do
                        if #(z*s) lt #z then
                            _AddScaledTerm(~terms, z, mucoeff * scale * coeff);
                        end if;
                    end for;
                end if;
            end for;
            _RemoveZeros(~terms);
            return EltIHkeConstruct(C, terms);
        end if;
    end if;

    return false;
end intrinsic;


//////////////////////////////////////////////
// Canonical basis of the antispherical module

declare type ASModIHkeCan[EltIHke]: IHkeAlgBaseCan;

intrinsic CanonicalBasis(ASMod::IHkeASMod) -> ASModIHkeCan
{The canonical basis of the right antispherical module.}
    v := BaseRing(ASMod).1;
    return _GetOrCreateBasis(ASMod, ASModIHkeCan, "aC", "Canonical basis", Parabolic(ASMod), -v);
end intrinsic;

intrinsic _EltIHkeValidate(aC::ASModIHkeCan, elt::EltIHke)
{Only allow I-minimal elements.}
    I := Parabolic(FreeModule(aC));
    error if not forall(w){w : w -> _ in elt`Terms | IsMinimal(I, w)},
        w, "is not minimal with respect to", I;
end intrinsic;

intrinsic _ToBasis(aH::ASModIHkeStd, aC::ASModIHkeCan, w::GrpFPCoxElt) -> EltIHke
{Express aC(w) in the standard basis.}
    return _CanToStd(aH, aC, w, aH`Para, aH`Eig);
end intrinsic;

intrinsic _ToBasis(aC::ASModIHkeCan, aH::ASModIHkeStd, w::GrpFPCoxElt) -> EltIHke
{Express aH(w) in the canonical basis.}
    return _StdToCan(aH, aC, w);
end intrinsic;


//////////////////////////////////////////
// Canonical basis of the spherical module

declare type SModIHkeCan[EltIHke]: IHkeAlgBaseCan;

intrinsic CanonicalBasis(SMod::IHkeSMod) -> SModIHkeCan
{The canonical basis of the right antispherical module.}
    v := BaseRing(SMod).1;
    return _GetOrCreateBasis(SMod, SModIHkeCan, "sC", "Canonical basis", Parabolic(SMod), v^-1);
end intrinsic;

intrinsic _EltIHkeValidate(sC::SModIHkeCan, elt::EltIHke)
{Only allow I-minimal elements.}
    I := Parabolic(FreeModule(sC));
    error if not forall(w){w : w -> _ in elt`Terms | IsMinimal(I, w)},
        w, "is not minimal with respect to", I;
end intrinsic;

intrinsic _ToBasis(sH::SModIHkeStd, sC::SModIHkeCan, w::GrpFPCoxElt) -> EltIHke
{Express sC(w) in the standard basis.}
    return _CanToStd(sH, sC, w, sH`Para, sH`Eig);
end intrinsic;

intrinsic _ToBasis(sC::SModIHkeCan, sH::SModIHkeStd, w::GrpFPCoxElt) -> EltIHke
{Express sH(w) in the canonical basis.}
    return _StdToCan(sH, sC, w);
end intrinsic;
