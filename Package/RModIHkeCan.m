import "Base.m":
    _LaurentPolyRing,
    _v;
import "EltIHke.m":
    _EltIHkeConstruct,
    _AddScaled,
    _RemoveZeros,
    _AddScaledTerm,
    _IsUniTriangular;
import "RModIHke.m":
    _RModIHeckeBaseInit;

declare type RModIHkeCan[EltIHke]: RModIHkeBase;
declare attributes RModIHkeCan:
    CanInStdCache,
    StdInCanCache,
    MuCache;

intrinsic IHeckeRModuleCan(alg::AlgIHke, para::SeqEnum[RngIntElt]) -> RModIHkeStd
{Standard basis of the antispherical right module.}
    // TODO: Cache the basis on the algebra object.
    basis := New(RModIHkeCan);
    _RModIHeckeBaseInit(~basis, alg, "aC", "Canonical basis", para);
    basis`CanInStdCache := AssociativeArray(CoxeterGroup(alg));
    basis`StdInCanCache := AssociativeArray(CoxeterGroup(alg));
    basis`MuCache := AssociativeArray(CoxeterGroup(alg));
    return basis;
end intrinsic;

intrinsic 'eq'(H1::RModIHkeCan, H2::RModIHkeCan) -> BoolElt
{}
    return Parent(H1) eq Parent(H2) and Parabolic(H1) eq Parabolic(H2);
end intrinsic;


//////////////////////////
// Basis conversion

// Assume that elt is an element written in terms of the standard basis, and right-multiply
// by the element B(s).
function _RightMultCanGen(elt, s)
    assert ISA(Type(Parent(elt)), AlgIHkeStd);

    W := CoxeterGroup(Parent(elt));
    I := Parabolic(Parent(elt));
    terms := AssociativeArray();
    for w -> coeff in elt`Terms do
        ws := w * (W.s);
        if not IsMinimal(I, ws) then
            // Do nothing, since aH(w) * B(s) is zero.
        elif #w lt #ws then
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, _v * coeff);
        else
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, _v^-1 * coeff);
        end if;
    end for;
    _RemoveZeros(~terms);
    return _EltIHkeConstruct(Parent(elt), terms);
end function;


// Return the mu-coefficient vector for w.
intrinsic MuCoeff(aH::RModIHkeStd, aC::RModIHkeCan, w::GrpFPCoxElt) -> Assoc
{}
    assert Parent(w) eq CoxeterGroup(aC);
    assert IsMinimal(Parabolic(aH), w);

    if IsDefined(aC`MuCache, w) then
        return aC`MuCache[w];
    end if;

    W := CoxeterGroup(aC);
    if #w eq 0 then
        return AssociativeArray(W);
    end if;

    aCw := _IHkeProtToBasis(aH, aC, w);
    mu := AssociativeArray(W);
    for u -> coeff in aCw`Terms do
        v_coeff := Coefficient(coeff, 1);
        if v_coeff ne 0 then
            mu[u] := v_coeff;
        end if;
    end for;

    aC`MuCache[w] := mu;
    return mu;
end intrinsic;


intrinsic _IHkeProtToBasis(aH::RModIHkeStd, aC::RModIHkeCan, w::GrpFPCoxElt) -> EltIHke
{Express aC(w) in the standard basis aH.}
    if IsDefined(aC`CanInStdCache, w) then
        return aC`CanInStdCache[w];
    end if;

    if #w eq 0 then
        return aH.0;
    end if;

    // Take an arbitrary right descent of w.
    W := CoxeterGroup(aH);
    s := Rep(RightDescentSet(W, w));
    ws := w * W.s;

    // C(ws) = C(w)C(s) - sum[us < u]mu(u, ws) C(u)
    aCw := _IHkeProtToBasis(aH, aC, ws);
    aCwCs := _RightMultCanGen(aCw, s);
    for u -> coeff in MuCoeff(aH, aC, ws) do
        if #(u*W.s) lt #u then
            aCwCs -:= coeff * _IHkeProtToBasis(aH, aC, u);
        end if;
    end for;

    // Double-check unitriangularity
    assert2 _IsUniTriangular(aCwCs`Terms, w);

    aC`CanInStdCache[w] := aCwCs;
    return aCwCs;
end intrinsic;


// intrinsic _IHkeProtToBasis(C::AlgIHkeCan, H::AlgIHkeStd, w::GrpFPCoxElt) -> EltIHke
// {Express H(w) in the canonical basis.}
//     if IsDefined(C`StdInCanCache, w) then
//         return C`StdInCanCache[w];
//     end if;

//     if #w eq 0 then
//         return C.0;
//     end if;

//     W := CoxeterGroup(H);
//     terms := AssociativeArray(W);
//     terms[w] := _LaurentPolyRing ! 1;
//     Cw := _IHkeProtToBasis(H, C, w);
//     for u -> coeff in Cw`Terms do
//         if u eq w then
//             continue;
//         end if;

//         Cu := _IHkeProtToBasis(C, H, u);
//         _AddScaled(~terms, Cu`Terms, -coeff);
//     end for;
//     _RemoveZeros(~terms);
//     result := _EltIHkeConstruct(C, terms);

//     C`StdInCanCache[w] := result;
//     return result;
// end intrinsic;
