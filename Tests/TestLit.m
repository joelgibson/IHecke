SetQuitOnError(true);
SetColumns(0);
SetAssertions(3);
AttachSpec("IHecke.spec");

W := CoxeterGroup(GrpFPCox, "C2");
HAlg := IHeckeAlgebra(W);
LPoly<v> := BaseRing(HAlg);
H := StandardBasis(HAlg);
C := CanonicalBasis(HAlg);
p2C := LiteralBasisInCan(HAlg, "p2C", "2-canonical basis");

// Canonical basis accessors
assert Sprint(p2C) eq "2-canonical basis of Iwahori-Hecke algebra of type B2, symbol p2C";
assert FreeModule(p2C) eq HAlg;
assert CoxeterGroup(p2C) eq W;
assert BasisSymbol(p2C) eq "p2C";
assert BasisName(p2C) eq "2-canonical basis";
assert BaseRing(p2C) eq BaseRing(HAlg);
assert p2C eq p2C;

// Set a partial canonical basis (as if we're calculating it using some algorithm).
SetBasisElement(p2C, W.0, C.0);
assert C ! p2C.0 eq C.0;

SetBasisElement(p2C, W.1, C.1);
SetBasisElement(p2C, W.2, C.2);
SetBasisElement(p2C, W![1,2], C.[1,2]);
SetBasisElement(p2C, W![2,1], C.[2,1]);
SetBasisElement(p2C, W![1,2,1], C.[1,2,1] + C.1);

assert C ! p2C.[1,2,1] eq C.[1,2,1] + C.1;
assert p2C ! C.[1,2,1] eq p2C.[1,2,1] - p2C.1;

print "TestLit passed.";
quit;
