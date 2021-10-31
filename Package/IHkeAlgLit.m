import "EltIHke.m": _AddScaled, _AddScaledTerm, _RemoveZeros, _IsUniTriangular;

declare type IHkeAlgLit: BasisIHke;
declare attributes IHkeAlgLit:
    BasisInCan,
    CanInBasis;

intrinsic LiteralBasisInCan(HAlg::IHkeAlg, symbol::MonStgElt, name::MonStgElt) -> IHkeAlgLit
{A basis of the Hecke algebra which is declared in terms of another basis. This basis may
 also be partial (not fully defined), but we require that if w is defined, then all
 elements less than w are defined. A basis element may never be re-defined. This basis
 stores its elements in terms of the canonical basis, and must be unitriangular to the
 canonical basis: L.w = C.w + (linear combination of C.x for x < w).}
    basis := New(IHkeAlgLit);
    _BasisIHkeInit(~basis, HAlg, symbol, name);
    basis`BasisInCan := AssociativeArray(CoxeterGroup(HAlg));
    basis`CanInBasis := AssociativeArray(CoxeterGroup(HAlg));
    return basis;
end intrinsic;

intrinsic SetBasisElement(L::IHkeAlgLit, w::GrpFPCoxElt, Lw::EltIHke)
{Set L.w = Lw. Throws an error if L.w is already defined.}
    W := CoxeterGroup(L);
    C := CanonicalBasis(FreeModule(L));
    require Parent(w) eq W: "Coxeter group has incorrect parent";
    require not IsDefined(L`BasisInCan, w): "Basis is already defined at", w;
    // Should we check that it is defined on all descents or something?

    // Set the basis element, and the inverse map.
    L`BasisInCan[w] := C ! Lw;

    terms := AssociativeArray(W);
    terms[w] := BaseRing(L) ! 1;
    for u -> coeff in L`BasisInCan[w]`Terms do
        if u eq w then
            continue;
        end if;
        error if not IsDefined(L`CanInBasis, u),
            "Could not invert basis transformation (missing lower terms)";

        _AddScaled(~terms, L`CanInBasis[u]`Terms, -coeff);
    end for;
    _RemoveZeros(~terms);
    L`CanInBasis[w] := EltIHkeConstruct(L, terms);
end intrinsic;

intrinsic _ToBasis(C::IHkeAlgCan, L::IHkeAlgLit, w::GrpFPCoxElt) -> EltIHke
{Convert from the literal basis to the canonical basis.}
    error if not IsDefined(L`BasisInCan, w),
        "The basis", L`name, "is not yet defined at", w;
    return L`BasisInCan[w];
end intrinsic;

intrinsic _ToBasis(L::IHkeAlgLit, C::IHkeAlgCan, w::GrpFPCoxElt) -> EltIHke
{Convert from the canonical basis to the literal basis.}
    error if not IsDefined(L`CanInBasis, w),
        "The basis", L`name, "is not yet defined at", w;
    return L`CanInBasis[w];
end intrinsic;


// We also need to define a conversion into the standard basis.
intrinsic _ToBasis(H::IHkeAlgStd, L::IHkeAlgLit, w::GrpFPCoxElt) -> EltIHke
{}
    C := CanonicalBasis(FreeModule(H));
    return H ! _ToBasis(C, L, w);
end intrinsic;
intrinsic _ToBasis(L::IHkeAlgLit, H::IHkeAlgStd, w::GrpFPCoxElt) -> EltIHke
{}
    C := CanonicalBasis(FreeModule(H));
    return ToBasis(L, C, C ! H.w);
end intrinsic;
