// Examples/BenchmarkCanInStd.m
//
// Generates the whole canonical basis in the standard basis, then exits. Can be used to benchmark how fast canonical
// bases can be calculated:
//  $ time magma -b type:=B4 Examples/BenchmarkCanInStd.m

AttachSpec("IHecke.spec");
SetColumns(0);
SetQuitOnError(true);
SetVerbose("IHecke", 3);

procedure Usage(missingArg)
    print "";
    printf "Error: missing argument %o\n", missingArg;
    print "Usage example: magma -b type:=B4 Examples/BenchmarkCanInStd.m";
    print "";
    print "Arguments:";
    print "  type    (required) Cartan type, eg A2, B3, F4, ...";
    print "  notable (optional) If set, do not load the mu-coefficients from a table.";
end procedure;

if not assigned type then
    Usage("type");
    quit;
end if;

W := CoxeterGroup(GrpFPCox, type);
HAlg := IHeckeAlgebra(W);
H := StandardBasis(HAlg);
C := CanonicalBasis(HAlg : UseTable := not assigned notable);

Welts := EnumerateCoxeterGroup(W);
printf "Generating the canonical basis for %o elements...\n", #Welts;
time for w in Welts do
    _ := H ! C.w;
end for;

quit;
