// Sometimes mu-coefficients are loaded from a table. This test double-checks that the canonical-to-standard basis
// conversion looks the same for the precalculated table mus and the on-the-fly calculated mus.
//
// This is a long-running test and should be run manually - one small case gets run in the ./run-tests.sh script to
// make sure that the test still executes (eg no functions have been renamed).

AttachSpec("IHecke.spec");
SetColumns(0);
SetQuitOnError(true);
SetVerbose("IHecke", 3);

procedure Usage(missingArg)
    print "";
    printf "Error: missing argument %o\n", missingArg;
    print "Usage example: magma -b type:=B4 Tests/LongTestMusCorrect.m";
    print "";
    print "Arguments:";
    print "  type    (required) Cartan type, eg A2, B3, F4, ...";
end procedure;

if not assigned type then
    Usage("type");
    quit;
end if;

W := CoxeterGroup(GrpFPCox, type);

// Make two copies of the Hecke algebra, the first with mu-coefficients loaded from a table, and the second without.
HAlg, H, C := ShortcutIHeckeAlgebra(W);
HAlg2 := IHeckeAlgebra(W);
H2 := StandardBasis(HAlg2);
C2 := CanonicalBasis(HAlg2 : UseTable := false);

for w in EnumerateCoxeterGroup(W) do
    elts1, coeffs1 := Support(H ! C.w);
    elts2, coeffs2 := Support(H2 ! C2.w);
    assert elts1 eq elts2;
    assert coeffs1 eq coeffs2;
end for;

printf "LongTestMusCorrect in type %o passed.\n", type;
quit;
