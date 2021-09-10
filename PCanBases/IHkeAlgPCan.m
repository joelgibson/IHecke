import "../Package/EltIHke.m":
    _AddScaled,
    _RemoveZeros,
    _AddScaledTerm,
    _IsUniTriangular;
import "PCanDB.m": PCanDB;


// Our "implementation" of the p-canonical basis is really just reading pre-calculated bases from a
// file: it is a "literal" basis in the sense that the basis is just given by a lookup table.
declare type IHkeAlgPCan[EltIHke]: BasisIHke;
declare attributes IHkeAlgPCan:
    // The prime and cartan type used.
    Prime,
    CartanName,

    // An associative array, mapping group elements to their p-canonical basis elements, expressed
    // in the canonical basis.
    InCan,

    // The inverse mapping, which we calculate on-the-fly and cache the results here.
    FromCanCache;


////////////////
// Construction

// Sometimes, the p-canonical basis is the same as the canonical basis. Times when this is
// definitely true are encoded in the pCanIsDefinitelyCan function:
//
// Type G2:
//  - [JW2016]: For primes larger than 3, p-can = can.
//  - [JW2016]: Since Cartan mat is symmetric mod 2, the 2-canonical basis is symmetric under swapping s, t.
//
// Types B3, C3:
//  - [JW2016]: The only interesting case is p = 2.
//
// Type D4:
//  - [JW2016]: The p-canonical basis coincides with the canonical basis except at p = 2.
//
// Type An:
//  - [WB12]: The p-canonical and canonical bases coincide for all n <= 6.
//  - [JW2016]: For A7, the p-canonical and canonical bases coincide except at p = 2.
//
// [JW2016]: "The p-canonical basis for Hecke algebras", Jensen and Williamson.
// [WB12] "Modular intersection cohomology complexes on flag varieties", Williamson and Braden.
function pCanIsDefinitelyCan(type, rank, prime)
    if type eq "A" then
        if rank le 6 then return true; end if;
        return false;
    elif type eq "B" then
        if rank in [2, 3] then return prime ne 2; end if;
        return false;
    elif type eq "C" then
        if rank in [2, 3] then return prime ne 2; end if;
        return false;
    elif type eq "D" then
        if rank in [4] then return prime ne 2; end if;
    elif type eq "F" and rank eq 4 then
        return false;
    elif type eq "G" and rank eq 2 then
        return prime ge 5;
    end if;
end function;

intrinsic IHeckeAlgebraPCan(HAlg::IHkeAlg, cartanName::MonStgElt, prime::RngIntElt: quiet := false) -> BasisIHke
{Load a p-canonical basis from a file, or check against a list of rules to check if this p-canonical
 basis is equal to the canonical basis. If neither of these can be determined, throw an error.

 This function will print a message when the p-canonical basis has been successfully loaded from a
 file, or it has instead returned a canonical basis (due to being told that they were equal). This
 message can be suppressed by passing quiet := true.

 Note that the cartanName needs to be specified, and be compatible with the Coxeter group defining
 the Hecke algebra A. It cannot be inferred from A, since the p-canonical basis depends on the root
 system rather than just the Coxeter group (B3 and C3 have different 2-canonical bases, for example).}

    error if not IsPrime(prime), "The number", prime, "provided was not prime.";

    W := CoxeterGroup(HAlg);
    error if CoxeterMatrix(W) ne CoxeterMatrix(cartanName),
        "The given Coxeter group is incompatible with the cartan type", cartanName;

    // Check if we are certain that the p-canonical and canonical bases coincide.
    type := cartanName[1];
    rank := StringToInteger(cartanName[2..#cartanName]);
    if pCanIsDefinitelyCan(type, rank, prime) then
        if not quiet then
            printf "The %o-canonical basis for type %o coincides with the canonical basis.\n",
                prime, cartanName;
        end if;
        return CanonicalBasis(HAlg);
    end if;

    // Check if we have the basis available in the database.
    code := Sprintf("Type-%o-P-%o", cartanName, prime);
    if not IsDefined(PCanDB, code) then
        error Sprintf("The %o-canonical basis for %o is unavailable", prime, cartanName);
    end if;

    // Read the file
    C := CanonicalBasis(HAlg);
    inCan := AssociativeArray(W);
    for i -> line in PCanDB[code] do
        // Look for something like "212: (1)C(212) + (1)C(2)".
        ok, _, matches := Regexp("^(.+).*:[ ]*(.+)\s*$", line);
        error if not ok, Sprintf("Could not understand line %o:\n%o", i, line);

        // Parse out the left and right sides.
        w := matches[1] eq "id"
            select W.0
            else W ! [StringToInteger(matches[1][i]) : i in [1..#matches[1]]];
        inCan[w] := ReadEltIHke(C, matches[2]);
    end for;

    basis := New(IHkeAlgPCan);
    _BasisIHkeInit(~basis, HAlg, Sprintf("p%oC", prime), Sprintf("%o-canonical basis", prime));
    basis`Prime := prime;
    basis`CartanName := cartanName;
    basis`InCan := inCan;
    basis`FromCanCache := AssociativeArray(W);

    if not quiet then
        printf "The %o-canonical basis for type %o was loaded from the database.\n",
            prime, cartanName;
    end if;

    return basis;
end intrinsic;


////////////
// Overrides

intrinsic 'eq'(C1::IHkeAlgPCan, C2::IHkeAlgPCan) -> BoolElt
{}
    return FreeModule(C1) eq FreeModule(C2) and C1`Prime eq C2`Prime and C1`CartanName eq C2`CartanName;
end intrinsic;


//////////////////////////
// Basis change (protocol)
//
// We'll define conversions in and out of the canonical basis, rather than the standard basis. We
// still need to provide an explicit conversion to the standard basis, since it's the default basis
// of the Hecke algebra.

intrinsic _ToBasis(C::IHkeAlgCan, pC::IHkeAlgPCan, w::GrpFPCoxElt) -> EltIHke
{}
    return pC`InCan[w];
end intrinsic;

// For the reverse, use unitriangularity.
intrinsic _ToBasis(pC::IHkeAlgPCan, C::IHkeAlgCan, w::GrpFPCoxElt) -> EltIHke
{}
    if IsDefined(pC`FromCanCache, w) then
        return pC`FromCanCache[w];
    end if;

    if #w eq 0 then
        return pC.0;
    end if;

    W := CoxeterGroup(pC);
    terms := AssociativeArray(W);
    terms[w] := BaseRing(pC) ! 1;
    for u -> coeff in _ToBasis(C, pC, w)`Terms do
        if u eq w then
            continue;
        end if;
        _AddScaled(~terms, _ToBasis(pC, C, u)`Terms, -coeff);
    end for;
    _RemoveZeros(~terms);
    result := EltIHkeConstruct(pC, terms);

    pC`FromCanCache[w] := result;
    return result;
end intrinsic;

intrinsic _ToBasis(H::IHkeAlgStd, pC::IHkeAlgPCan, w::GrpFPCoxElt) -> EltIHke
{}
    C := CanonicalBasis(FreeModule(pC));
    return _ToBasis(H, C, _ToBasis(C, pC, w));
end intrinsic;

intrinsic _ToBasis(pC::IHkeAlgPCan, H::IHkeAlgStd, w::GrpFPCoxElt) -> EltIHke
{}
    C := CanonicalBasis(FreeModule(pC));
    return _ToBasis(pC, C, _ToBasis(C, H, w));
end intrinsic;

/////////////////////
// Input

intrinsic ReadEltIHke(B::BasisIHke, eltStr::MonStgElt) -> EltIHke
{Read a printed Hecke element. The symbol used for printing the basis ("H", "C", etc) are ignored,
 and instead the given parent is used as the basis in which to interpret the string.}
    W := CoxeterGroup(B);
    v := BaseRing(B).1; // Needed for evalling the coefficient.
    terms := AssociativeArray(W);
    line := eltStr;
    while true do
        ok, _, matches := Regexp("^\\(([^)]+)\\)[^(]*\\(([^)]+)\\)( [+] (.*))?$", eltStr);
        if not ok then
            error if not Regexp("^[ ]*$", eltStr), "Could not parse the element", line;
        end if;

        w := matches[2] eq "id"
            select W.0
            else W ! [StringToInteger(matches[2][i]) : i in [1..#matches[2]]];

        _AddScaledTerm(~terms, w, eval matches[1]);
        if #matches eq 2 then
            break;
        end if;
        eltStr := matches[4];
    end while;
    _RemoveZeros(~terms);
    return EltIHkeConstruct(B, terms);
end intrinsic;