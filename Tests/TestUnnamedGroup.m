// Test creation of the Hecke algebra for a non-standard group.

SetQuitOnError(true);
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

// The Coxeter graph for a triangle with all sides labelled by 5.
W := CoxeterGroup(GrpFPCox, Matrix(Integers(), [
    [1, 5, 5],
    [5, 1, 5],
    [5, 5, 1]
]));
HAlg := IHeckeAlgebra(W);
assert Sprint(HAlg) eq "Iwahori-Hecke algebra of type X3";

H := StandardBasis(HAlg);
C := CanonicalBasis(HAlg);
x := H ! C.[1,2,3,2,1];

print "TestUnnamedGroup passed.";
quit;
