// Examples/BenchmarkCanCanMultiplication.m
//
// Right-multiply every canonical basis element C(w) by the element C(1), expressing the result back in C.
//  $ time magma -b type:=B4 Examples/BenchmarkCanCanMultiplication.m

AttachSpec("IHecke.spec");
SetColumns(0);
SetQuitOnError(true);
SetVerbose("IHecke", 3);

procedure Usage(missingArg)
    print "";
    printf "Error: missing argument %o\n", missingArg;
    print "Usage example: magma -b type:=B4 Examples/BenchmarkCanCanMultiplication.m";
    print "";
    print "Arguments:";
    print "  type    (required) Cartan type, eg A2, B3, F4, ...";
    print "  notable (optional) If set, do not load the mu-coefficients from a table.";
    print "  profile (optional) Profile filename.";
end procedure;

if not assigned type then
    Usage("type");
    quit;
end if;

if assigned profile then
    SetProfile(true);
end if;

W := CoxeterGroup(GrpFPCox, type);
HAlg := IHeckeAlgebra(W);
H := StandardBasis(HAlg);
C := CanonicalBasis(HAlg : UseTable := not assigned notable);

Welts := EnumerateCoxeterGroup(W);
printf "Right multiplying all %o canonical basis elements by C(s)...\n", #Welts;
time for w in Welts do
    _ := C.w * C.1;
end for;

if assigned profile then
    G := ProfileGraph();
    fname := profile cat ".html";
    ProfileHTMLOutput(profile);
    printf "Profile written to %o.html\n", fname;
end if;


quit;
