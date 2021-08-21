import "Base.m":
    _LaurentPolyRing,
    _v;
import "EltIHke.m":
    _AddScaled,
    _RemoveZeros,
    _AddScaledTerm,
    _IsUniTriangular;
import "AlgIHkeBase.m":
    _AlgIHkeBaseInit;


declare type AlgIHkeCan[EltIHke]: AlgIHkeBase;
declare attributes AlgIHkeCan:
    CanInStdCache,
    StdInCanCache,
    MuCache;


//////////////////////////////////
// Creation

intrinsic IHeckeAlgebraCan(alg::AlgIHke) -> AlgIHkeCan
{The canonical basis of the Hecke algebra.}
    if not IsDefined(alg`BasisCache, AlgIHkeCan) then
        basis := New(AlgIHkeCan);
        _AlgIHkeBaseInit(~basis, alg, "C", "Canonical basis");
        basis`CanInStdCache := AssociativeArray(CoxeterGroup(alg));
        basis`StdInCanCache := AssociativeArray(CoxeterGroup(alg));
        basis`MuCache := AssociativeArray(CoxeterGroup(alg));
        alg`BasisCache[AlgIHkeCan] := basis;
    end if;
    return alg`BasisCache[AlgIHkeCan];
end intrinsic;


////////////
// Overrides

intrinsic 'eq'(C1::AlgIHkeCan, C2::AlgIHkeCan) -> BoolElt
{}
    // If we have two bases of type AlgIHkeCan of the same Hecke algebra, they must be equal.
    return Parent(C1) eq Parent(C2);
end intrinsic;


//////////////////////////
// Basis conversion

// Assume that elt is an element written in terms of the standard basis, and right-multiply
// by the element B(s).
function _RightMultCanGen(elt, s)
    assert ISA(Type(Parent(elt)), AlgIHkeStd);

    W := CoxeterGroup(Parent(elt));
    terms := AssociativeArray();
    for w -> coeff in elt`Terms do
        ws := w * (W.s);
        if #w lt #ws then
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, _v * coeff);
        else
            _AddScaledTerm(~terms, ws, coeff);
            _AddScaledTerm(~terms, w, _v^-1 * coeff);
        end if;
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(Parent(elt), terms);
end function;


// Return the mu-coefficient vector for w.
intrinsic MuCoeff(H::AlgIHkeStd, C::AlgIHkeCan, w::GrpFPCoxElt) -> Assoc
{}
    assert Parent(w) eq CoxeterGroup(C);

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
end intrinsic;


intrinsic _IHkeProtToBasis(H::AlgIHkeStd, C::AlgIHkeCan, w::GrpFPCoxElt) -> EltIHke
{Express C(w) in the standard basis.}
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
    CwCs := _RightMultCanGen(Cw, s);
    for u -> coeff in MuCoeff(H, C, ws) do
        if #(u*W.s) lt #u then
            CwCs -:= coeff * _IHkeProtToBasis(H, C, u);
        end if;
    end for;

    // Double-check unitriangularity
    assert2 _IsUniTriangular(CwCs`Terms, w);

    C`CanInStdCache[w] := CwCs;
    return CwCs;
end intrinsic;


intrinsic _IHkeProtToBasis(C::AlgIHkeCan, H::AlgIHkeStd, w::GrpFPCoxElt) -> EltIHke
{Express H(w) in the canonical basis.}
    if IsDefined(C`StdInCanCache, w) then
        return C`StdInCanCache[w];
    end if;

    if #w eq 0 then
        return C.0;
    end if;

    W := CoxeterGroup(H);
    terms := AssociativeArray(W);
    terms[w] := _LaurentPolyRing ! 1;
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
end intrinsic;


////////////////////////////
// Bar involution (protocol)

intrinsic _IHkeProtBar(C::AlgIHkeCan, elt::EltIHke) -> EltIHke
{The bar involution of elt, fixing each basis element and twisting scalars by by v -> v^-1.}
    assert C eq Parent(elt);

    twist := hom<_LaurentPolyRing -> _LaurentPolyRing | _v^-1>;
    terms := AssociativeArray(CoxeterGroup(C));
    for w -> coeff in elt`Terms do
        terms[w] := twist(coeff);
    end for;

    return EltIHkeConstruct(C, terms);
end intrinsic;