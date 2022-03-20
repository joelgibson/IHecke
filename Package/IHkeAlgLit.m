import "EltIHke.m": _AddScaled, _RemoveZeros;

declare type IHkeAlgLit: BasisIHke;
declare attributes IHkeAlgLit:
    TargetBasis,    // "Standard" or "Canonical".
    IntoTarget,     // Assoc, changes from the literal basis into the target basis.
    OutOfTarget;    // Assoc, changes from the target basis into the literal basis.

// TODO: Clean this up
// Plan:
//   - Literal basis should be triangular to either the standard or canonical basis (user may choose).
//   - Needs to work for spherical and antispherical modules.
//   - Needs to be saveable (save must encode group, free module, and coefficients).
//     - Have a function for serialising a basis (returns a serialisable Magma object) and restoring. Let someone else
//       deal with actually saving. (Put a code snippet in the README).

intrinsic CreateLiteralBasis(FMod::FModIHke, TargetBasis::MonStgElt, symbol::MonStgElt, name::MonStgElt) -> IHkeAlgLit
{A basis which is defined "literally", i.e. by a look-up table in terms of another basis. The TargetBasis should be
 either "Standard" or "Canonical", which may affect the performance of the basis, but is otherwise invisible.}
    require TargetBasis in ["Standard", "Canonical"]:
        "The target basis should be either \"Standard\" or \"Canonical\"";

    basis := New(IHkeAlgLit);
    _BasisIHkeInit(~basis, FMod, symbol, name);
    basis`TargetBasis := TargetBasis;
    basis`IntoTarget := AssociativeArray(CoxeterGroup(FMod));
    basis`OutOfTarget := AssociativeArray(CoxeterGroup(FMod));
    return basis;
end intrinsic;

intrinsic LiteralBasisInCan(HAlg::IHkeAlg, symbol::MonStgElt, name::MonStgElt) -> IHkeAlgLit
{A basis of the Hecke algebra which is declared in terms of another basis. This basis may
 also be partial (not fully defined), but we require that if w is defined, then all
 elements less than w are defined. A basis element may never be re-defined. This basis
 stores its elements in terms of the canonical basis, and must be unitriangular to the
 canonical basis: L.w = C.w + (linear combination of C.x for x < w).}
    basis := New(IHkeAlgLit);
    _BasisIHkeInit(~basis, HAlg, symbol, name);
    basis`TargetBasis := "Canonical";
    basis`IntoTarget := AssociativeArray(CoxeterGroup(HAlg));
    basis`OutOfTarget := AssociativeArray(CoxeterGroup(HAlg));
    return basis;
end intrinsic;

// Check if a Laurent polynomial is a unit, i.e. of the form ±v^k for some k. We cannot use the built-in IsUnit function
// in Magma, since we are technically using Laurent series, and there are many more units there.
function IsLaurentPolyUnit(poly)
    intcoeffs := Coefficients(poly);
    return #intcoeffs eq 1 and Abs(intcoeffs[1]) eq 1;
end function;

intrinsic SetBasisElement(L::IHkeAlgLit, w::GrpFPCoxElt, Lw::EltIHke)
{Set L.w = Lw. Throws an error if L.w is already defined.}
    W := CoxeterGroup(L);
    require Parent(w) eq W: "Coxeter group has incorrect parent";
    require not IsDefined(L`IntoTarget, w): "Basis is already defined at", w;

    // Check unit condition on top term.
    require IsLaurentPolyUnit(Coefficient(Lw, w)): "The coefficient", Coefficient(Lw, w), "of", w, "is not a unit.";


    B := L`TargetBasis eq "Standard"
        select StandardBasis(FreeModule(L))
          else CanonicalBasis(FreeModule(L));


    // TODO: Should we check that it is defined on all descents or something?

    // Set the basis element, and the inverse map.
    L`IntoTarget[w] := B ! Lw;

    terms := AssociativeArray(W);
    terms[w] := Coefficient(L`IntoTarget[w], w)^-1;
    for u -> coeff in L`IntoTarget[w]`Terms do
        if u eq w then
            continue;
        end if;
        error if not IsDefined(L`OutOfTarget, u),
            "Could not invert basis transformation (missing lower terms)";

        _AddScaled(~terms, L`OutOfTarget[u]`Terms, -coeff);
    end for;
    _RemoveZeros(~terms);
    L`OutOfTarget[w] := EltIHkeConstruct(L, terms);
end intrinsic;

intrinsic IsDefined(L::IHkeAlgLit, w::GrpFPCoxElt) -> BoolElt
{Check if the literal basis is defined at w, i.e. L.w can be converted into other bases.}
    require Parent(w) eq CoxeterGroup(L): w, "is not a member of", CoxeterGroup(L);
    return result where result := IsDefined(L`IntoTarget, w);
end intrinsic;

intrinsic Keys(L::IHkeAlgLit) -> SetEnum[GrpFPCoxElt]
{Return the set of keys w for which IsDefined(L, w) returns true.}
    return Keys(L`IntoTarget);
end intrinsic;

intrinsic _ToBasis(B::BasisIHke, L::IHkeAlgLit, w::GrpFPCoxElt) -> EltIHke
{Convert from the literal basis into some other basis.}
    error if not IsDefined(L`IntoTarget, w),
        "The basis", BasisName(L), "is not yet defined at", w;

    if L`TargetBasis eq "Standard" then
        return ISA(Type(B), IHkeAlgBaseStd)
            select L`IntoTarget[w]
              else B ! L`IntoTarget[w];
    else
        return ISA(Type(B), IHkeAlgBaseCan)
            select L`IntoTarget[w]
              else B ! L`IntoTarget[w];
    end if;
end intrinsic;

intrinsic _ToBasis(L::IHkeAlgLit, B::BasisIHke, eltB::EltIHke) -> EltIHke
{Convert from some other basis into the literal basis.}
    TargetBasis := L`TargetBasis eq "Standard"
        select StandardBasis(FreeModule(L))
          else CanonicalBasis(FreeModule(L));
    eltT := TargetBasis ! eltB;
    error if not forall(w){w : w -> _ in eltT`Terms | IsDefined(L`OutOfTarget, w)},
        "The basis", BasisName(L), "is not yet defined at", w;

    terms := AssociativeArray(CoxeterGroup(L));
    for w -> coeff in eltT`Terms do
        _AddScaled(~terms, L`OutOfTarget[w]`Terms, coeff);
    end for;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(L, terms);
end intrinsic;

/////////////////////////////////////
// Serialisation of the Literal basis

// The serialisation format should work at Magma print level, and also using the WriteObject API. We return a record
// containing enough data to identify the correct Hecke algebra or module. A bijection is set up between the relevant
// Coxeter group elements and the integers by recording reduced expressions in a list, and then the transition matrices
// are recorded as sparse matrices under this bijection. (Empty rows mean the basis is not defined at that element).
SerialFmt := recformat<
    CoxeterMatrix : Mtrx,       // The Coxeter matrix underlying the Hecke algebra or module.
    FreeModuleType : MonStgElt, // "Hecke algebra", "Right antispherical module", or "Right spherical module".
    Parabolic : SeqEnum,        // For the spherical or antispherical modules, the parabolic subset used.
    TargetBasis : MonStgElt,    // "Standard" or "Canonical".
    BasisName : MonStgElt,      // The user-assigned name for this basis.
    BasisSymbol : MonStgElt,    // A short string naming the basis.
    Rexes : SeqEnum,            // A list of reduced expressions, establishing a bijection with integers.
    IntoTarget : MtrxSprs,      // The row indexed by w gives the coefficients of the literal basis in the target basis.
    OutOfTarget : MtrxSprs      // The inverse of IntoTarget.
>;

intrinsic SerialiseBasis(L::IHkeAlgLit) -> Rec
{Return the basis as a serialisable Magma object, for saving to a file. Restore using DeserialiseBasis.}
    FMod := FreeModule(L);
    W := CoxeterGroup(L);
    LPoly := BaseRing(FMod);

    // Record all Coxeter group elements used into an indexed set.
    grpelts := {@ W | @};
    for w -> elt in L`IntoTarget do
        Include(~grpelts, w);
        for x -> _ in elt`Terms do
            Include(~grpelts, x);
        end for;
    end for;
    for w -> elt in L`OutOfTarget do
        Include(~grpelts, w);
        for x -> _ in elt`Terms do
            Include(~grpelts, x);
        end for;
    end for;

    // Prepare the IntoTarget and OutOfTarget matrices. There may be empty rows.
    IntoTarget := SparseMatrix(LPoly, #grpelts, #grpelts);
    for w -> elt in L`IntoTarget do
        for x -> coeff in elt`Terms do
            SetEntry(~IntoTarget, Index(grpelts, w), Index(grpelts, x), coeff);
        end for;
    end for;

    OutOfTarget := SparseMatrix(LPoly, #grpelts, #grpelts);
    for w -> elt in L`OutOfTarget do
        for x -> coeff in elt`Terms do
            SetEntry(~OutOfTarget, Index(grpelts, w), Index(grpelts, x), coeff);
        end for;
    end for;

    return rec<SerialFmt |
        CoxeterMatrix := CoxeterMatrix(W),
        FreeModuleType := FreeModuleType(FMod),
        Parabolic := Parabolic(FMod),
        TargetBasis := "Canonical",
        BasisName := BasisName(L),
        BasisSymbol := BasisSymbol(L),
        Rexes := [Eltseq(w) : w in grpelts],
        IntoTarget := IntoTarget,
        OutOfTarget := OutOfTarget
    >;
end intrinsic;

intrinsic DeserialiseBasis(FMod::FModIHke, ser::Rec) -> IHkeAlgLit
{Restore a basis which was serialised using SerialiseBasis.}
    W := CoxeterGroup(FMod);

    require CoxeterMatrix(W) eq ser`CoxeterMatrix:
        "The serialised Coxeter matrix is incompatible.";
    require ser`FreeModuleType eq FreeModuleType(FMod):
        Sprintf("Basis to be deserialised is for \"%o\", not for \"%o\"", FreeModuleType(FMod), ser`FreeModuleType);
    require ser`Parabolic eq Parabolic(FMod):
        Sprintf("Basis to be deserialised has incompatible parabolic subset %o", ser`Parabolic);

    // TODO: Something with target basis, once we allow a different target basis.

    grpelts := {@ W | W ! rex : rex in ser`Rexes @};
    L := LiteralBasisInCan(FMod, ser`BasisSymbol, ser`BasisName);

    procedure unpack(~assoc, basis, matrix)
        M := ChangeRing(matrix, BaseRing(FMod));
        for w -> weight in RowWeights(M) do
            if weight eq 0 then
                continue;
            end if;

            assert IsLaurentPolyUnit(M[w, w]);
            terms := AssociativeArray(W);
            for x in Support(M, w) do
                terms[grpelts[x]] := M[w, x];
            end for;
            assoc[grpelts[w]] := EltIHkeConstruct(basis, terms);
        end for;
    end procedure;

    basis := L`TargetBasis eq "Standard" select StandardBasis(FMod) else CanonicalBasis(FMod);
    unpack(~L`IntoTarget, basis, ser`IntoTarget);
    unpack(~L`OutOfTarget, L, ser`OutOfTarget);

    return L;
end intrinsic;
